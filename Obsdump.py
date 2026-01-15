#!/awips2/python/bin/python
"""
obsdump - part of MatchObsAll routines - version 0.215

Searches through METAR and mesonet netcdf files, and creates ASCII
comma-delimited files of Temp, Dewpoint, Wind Direction, Wind speed,
Wind Gust. If there are many obs in a single hour, the one nearest
the "top of the hour" is saved.

Author: Tim Barker - SOO BOI
  2025/01/15 - version 0.215. KCA. Refactored for clarity and maintainability.
  2021/04/21 - version 0.207. MBC. Added support for WaveHeight.
  ... (previous history)
"""
import datetime
import os
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import dateutil.tz

import SmartScript
from ufpy.dataaccess import DataAccessLayer as DAL
from shapely.geometry import MultiPoint
import AbsTime
import TimeRange

import MatchObsAllConfigCR as MOAC

# =============================================================================
# Constants
# =============================================================================

GMT_ZONE = dateutil.tz.gettz('GMT')
MenuItems = None

# Default config values
CONFIG_DEFAULTS = {
    "DEBUG": 0,
    "PREVHOURS": 8,
    "OBSDIR": "/localapps/runtime/MatchObsAll/data/XXX",
    "sortareas": [],
    "DUMPAREA": None,
}

# Apply defaults
for key, value in CONFIG_DEFAULTS.items():
    if key not in MOAC.Config:
        MOAC.Config[key] = value

# Observation source configurations
OB_SOURCES = {
    "metar": {
        "request_type": "obs",
        "parameters": [
            "elevation", "tempFromTenths", "temperature",
            "dpFromTenths", "dewpoint", "windSpeed",
            "windDir", "windGust", "pkwndSpeed", "visibility",
            "skyCover", "skyLayerBase"
        ],
        "temp_unit": "C",
        "wind_unit": "kt",
    },
    "mesonet": {
        "request_type": "ldadmesonet",
        "parameters": [
            "elevation", "temperature", "dewpoint",
            "windDir", "windSpeed", "windGust"
        ],
        "temp_unit": "K",
        "wind_unit": "m/s",
    },
    "maritime": {
        "request_type": "sfcobs",
        "parameters": [
            "elevation", "temperature", "dewpoint",
            "windDir", "windSpeed", "windGust", "waveHeight"
        ],
        "temp_unit": "K",
        "wind_unit": "m/s",
    },
}

SKY_COVER_VALUES = {"VV": 100, "OVC": 95, "BKN": 75, "SCT": 50, "FEW": 25, "CLR": 5}
CEILING_COVERS = ("BKN", "OVC", "VV")
NON_CEILING_COVERS = ("SCT", "FEW")


# =============================================================================
# Logging
# =============================================================================

def logtime(message: str, priority: int = 0) -> None:
    """Write a timestamped log message."""
    if priority > MOAC.Config["DEBUG"]:
        return
    ts = datetime.datetime.now(GMT_ZONE).strftime("%Y/%m/%d %H:%M:%S")
    print(f"{ts}| {message}")
    sys.stdout.flush()


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class StationData:
    """Station metadata."""
    sid: str
    name: str
    lat: float
    lon: float
    _elev: float = 0

    @property
    def elev(self) -> int:
        return self._elev

    @elev.setter
    def elev(self, value: float) -> None:
        if value > 20000:
            self._elev = 20000
        elif value < -1000:
            logtime(f"Warning: elevation {value} out of range, clamping to -1000")
            self._elev = -1000
        else:
            self._elev = int(value)


@dataclass
class Observation:
    """Single observation record."""
    stnid: str
    obtime: datetime.datetime
    offset: float = float("inf")

    temp: Optional[int] = None
    dpt: Optional[int] = None
    wdir: Optional[int] = None
    wspd: Optional[int] = None
    pkwd: Optional[int] = None
    gust: Optional[int] = None
    vsby: Optional[float] = None
    wvehght: Optional[int] = None

    sky_cover: List[str] = field(default_factory=list)
    sky_layer_base: List[int] = field(default_factory=list)

    def reset_sky(self) -> None:
        """Clear sky cover data."""
        self.sky_cover = []
        self.sky_layer_base = []

    def add_sky_layer(self, cover: str, base: float) -> None:
        """Add a sky layer observation."""
        if cover in SKY_COVER_VALUES:
            self.sky_cover.append(cover)
            if cover == "CLR":
                self.sky_layer_base.append(25000)
            elif base > -9000:
                self.sky_layer_base.append(int(base))
        elif cover:
            logtime(f"Unknown sky cover: {cover}")

    @property
    def skyamt(self) -> Optional[int]:
        """Get sky amount percentage."""
        for cover, value in SKY_COVER_VALUES.items():
            if cover in self.sky_cover:
                return value
        return None

    @property
    def cig(self) -> Optional[int]:
        """Get ceiling height in hundreds of feet."""
        if not self.sky_cover or not self.sky_layer_base:
            return None

        # Sort layers by height
        try:
            layers = sorted(zip(self.sky_layer_base, self.sky_cover))
        except (TypeError, ValueError):
            return None

        for height, cover in layers:
            if cover in CEILING_COVERS:
                return height // 100
            if cover in NON_CEILING_COVERS:
                return height // 100

        return 250  # Clear sky default

    @property
    def effective_gust(self) -> Optional[int]:
        """Get gust, considering peak wind if higher."""
        gust = self.gust if self.gust and self.gust > 0 else None
        pkwd = self.pkwd

        if pkwd and gust:
            return max(pkwd, gust)
        return pkwd or gust

    def format_value(self, value: Optional[float], width: int = 3) -> str:
        """Format a value for output."""
        if value is None:
            return " " * width
        return f"{int(round(value)):{width}d}"

    def format_vis(self) -> str:
        """Format visibility (in 100ths of miles)."""
        if self.vsby is None:
            return "    "
        return f"{int(round(self.vsby * 100)):4d}"


# =============================================================================
# Unit Conversions
# =============================================================================

def c_to_f(t: float) -> float:
    return (t * 9.0 / 5.0) + 32.0

def k_to_f(t: float) -> float:
    return (t - 273.15) * 9.0 / 5.0 + 32.0

def ms_to_kt(speed: float) -> float:
    return speed * 1.943

def m_to_ft(height: float) -> float:
    return height * 3.281


def convert_temp(value: float, unit: str) -> Optional[int]:
    """Convert temperature to Fahrenheit."""
    if value is None or value <= -9000:
        return None

    if unit == 'C':
        result = c_to_f(value)
    elif unit == 'K':
        result = k_to_f(value)
    elif unit == 'F':
        result = value
    else:
        logtime(f"Unknown temperature unit: {unit}")
        return None

    result = int(round(result))
    if -80 <= result <= 130:
        return result

    logtime(f"Temperature out of range: {result}")
    return None


def convert_wind(value: float, unit: str) -> Optional[int]:
    """Convert wind speed to knots."""
    if value is None or value <= -9000:
        return None

    if unit == 'm/s':
        result = ms_to_kt(value)
    elif unit == 'kt':
        result = value
    else:
        logtime(f"Unknown wind unit: {unit}")
        return None

    result = int(round(result))
    if 0 <= result <= 150:
        return result

    logtime(f"Wind speed out of range: {result}")
    return None


# =============================================================================
# Main Procedure
# =============================================================================

class Procedure(SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)
        self.station_names: Dict[str, str] = {}
        self.stations: Dict[str, StationData] = {}
        self.obs: Dict[datetime.datetime, Dict[str, Observation]] = {}
        self.dt_list: List[datetime.datetime] = []
        self.sort_masks: List[np.ndarray] = []
        self.sort_names: List[str] = []
        self.dump_mask: Optional[np.ndarray] = None
        self.lat_min = self.lat_max = self.lon_min = self.lon_max = 0.0

    def execute(self):
        """Generate files with obs data nearest top of hour."""
        t0 = time.time()
        logtime("Starting obsdump")

        self._setup_sort_areas()
        self._setup_dump_area()
        self._setup_grid_limits()
        self._setup_time_list()
        self._load_station_names()

        # Fetch observations from all sources
        for source_name, source_config in OB_SOURCES.items():
            self._fetch_observations(source_name, source_config)
            logtime(f"Execution time after {source_name} = {time.time() - t0:.1f}s")

        self._write_output()
        logtime(f"Total execution time = {time.time() - t0:.1f}s")

    def _setup_sort_areas(self) -> None:
        """Initialize sorting edit areas."""
        for area_name in MOAC.Config["sortareas"]:
            logtime(f"Processing sort area: {area_name}")
            mask = self.encodeEditArea(area_name)
            if mask is not None:
                self.sort_masks.append(mask)
                self.sort_names.append(area_name)
            else:
                logtime(f"Edit area {area_name} is not a valid polygon")

    def _setup_dump_area(self) -> None:
        """Initialize the dump area mask."""
        if MOAC.Config["DUMPAREA"] is None:
            self.dump_mask = self.newGrid(1)
        else:
            mask = self.encodeEditArea(MOAC.Config["DUMPAREA"])
            self.dump_mask = mask.copy() if mask is not None else self.newGrid(1)

    def _setup_grid_limits(self) -> None:
        """Calculate lat/lon bounds from grid corners."""
        ny, nx = self.getGridShape()
        corners = [
            self.getLatLon(0, 0),
            self.getLatLon(0, ny - 1),
            self.getLatLon(nx - 1, 0),
            self.getLatLon(nx - 1, ny - 1),
        ]
        lats, lons = zip(*corners)
        self.lat_min, self.lat_max = float(min(lats)), float(max(lats))
        self.lon_min, self.lon_max = float(min(lons)), float(max(lons))

        logtime(f"Latitude limits: {self.lat_min:.3f} -- {self.lat_max:.3f}", 5)
        logtime(f"Longitude limits: {self.lon_min:.3f} -- {self.lon_max:.3f}", 5)

    def _setup_time_list(self) -> None:
        """Build list of hourly times to process."""
        now = datetime.datetime.now(GMT_ZONE).replace(minute=0, second=0, microsecond=0)
        delta_hour = datetime.timedelta(hours=1)

        for offset in range(MOAC.Config["PREVHOURS"], -1, -1):
            dt = now - (delta_hour * offset)
            self.obs[dt] = defaultdict(dict)
            self.dt_list.append(dt)

    def _load_station_names(self) -> None:
        """Load station names from database."""
        self.station_names = {}
        req = DAL.newDataRequest("common_obs_spatial")
        req.setParameters("stationid", "name")
        envelope = MultiPoint((
            (self.lon_max, self.lat_min),
            (self.lon_min, self.lat_max)
        ))
        req.setEnvelope(envelope)

        for geom in DAL.getGeometryData(req):
            sid = geom.getString("stationid")
            if sid not in self.station_names:
                self.station_names[sid] = geom.getString("name")

        logtime(f"Loaded {len(self.station_names)} station names")

    def _fetch_observations(self, source_name: str, config: dict) -> None:
        """Fetch observations from a single source type."""
        req = DAL.newDataRequest(config["request_type"])
        req.setParameters(*config["parameters"])
        envelope = MultiPoint((
            (self.lon_min, self.lat_min),
            (self.lon_max, self.lat_max)
        ))
        req.setEnvelope(envelope)

        delta_30min = datetime.timedelta(minutes=30)

        for dt in self.dt_list:
            begin_dt = dt - delta_30min
            end_dt = dt + delta_30min
            start = AbsTime.AbsTime(self._get_epoch(begin_dt))
            end = AbsTime.AbsTime(self._get_epoch(end_dt))
            tr = TimeRange.TimeRange(start, end)

            observations = DAL.getGeometryData(req, tr)
            if not observations:
                logtime(f"WARNING: No {source_name} obs retrieved for {dt}")
                continue

            logtime(f"{len(observations)} {source_name} obs retrieved for {dt}", 10)

            for ob in observations:
                self._process_observation(ob, dt, config)

    def _process_observation(self, ob, dt: datetime.datetime, config: dict) -> None:
        """Process a single observation record."""
        sid = ob.getLocationName()

        # Create station data if needed
        if sid not in self.stations:
            name = self.station_names.get(sid, sid)
            elev = ob.getNumber("elevation") * 3.2808 if "elevation" in ob.getParameters() else 0
            self.stations[sid] = StationData(
                sid=sid,
                name=name,
                lat=ob.getGeometry().y,
                lon=ob.getGeometry().x,
                _elev=int(min(20000, max(-1000, elev)))
            )

        ob_time = AbsTime.AbsTime(ob.getDataTime().getRefTime())
        this_offset = abs((dt - ob_time).total_seconds())

        # Check if we should update existing observation
        if self.obs[dt][sid]:
            if this_offset > self.obs[dt][sid].offset:
                return  # Existing ob is closer to top of hour
            if this_offset < self.obs[dt][sid].offset:
                self.obs[dt][sid].reset_sky()
        else:
            self.obs[dt][sid] = Observation(stnid=sid, obtime=ob_time)

        self.obs[dt][sid].offset = this_offset
        self.obs[dt][sid].obtime = ob_time

        # Update observation data
        self._set_ob_data(ob, dt, sid, config)

    def _set_ob_data(self, ob, dt: datetime.datetime, sid: str, config: dict) -> None:
        """Set observation data fields."""
        params = ob.getParameters()
        obs_record = self.obs[dt][sid]
        temp_unit = config["temp_unit"]
        wind_unit = config["wind_unit"]

        # Sky cover
        if 'skyCover' in params:
            cover = ob.getString('skyCover')
            base = ob.getNumber('skyLayerBase') if 'skyLayerBase' in params else -9999
            obs_record.add_sky_layer(cover, base)
            return  # Sky obs are separate records

        # Temperature (prefer tenths if available)
        if 'tempFromTenths' in params and ob.getNumber('tempFromTenths') > -900:
            obs_record.temp = convert_temp(ob.getNumber('tempFromTenths'), temp_unit)
        elif 'temperature' in params:
            obs_record.temp = convert_temp(ob.getNumber('temperature'), temp_unit)

        # Dewpoint (prefer tenths if available)
        if 'dpFromTenths' in params and ob.getNumber('dpFromTenths') > -900:
            obs_record.dpt = convert_temp(ob.getNumber('dpFromTenths'), temp_unit)
        elif 'dewpoint' in params:
            obs_record.dpt = convert_temp(ob.getNumber('dewpoint'), temp_unit)

        # Wind
        if 'windDir' in params:
            wdir = ob.getNumber('windDir')
            if wdir is not None and 0 <= wdir <= 360:
                obs_record.wdir = int(wdir)

        if 'windSpeed' in params:
            obs_record.wspd = convert_wind(ob.getNumber('windSpeed'), wind_unit)

        if 'windGust' in params:
            obs_record.gust = convert_wind(ob.getNumber('windGust'), wind_unit)

        if 'pkwndSpeed' in params:
            obs_record.pkwd = convert_wind(ob.getNumber('pkwndSpeed'), wind_unit)

        # Visibility
        if 'visibility' in params:
            vis = ob.getNumber('visibility')
            if vis is not None and 0 <= vis <= 10:
                obs_record.vsby = vis

        # Wave height
        if 'waveHeight' in params:
            wh = ob.getNumber('waveHeight')
            if wh is not None and wh > -9000:
                obs_record.wvehght = int(m_to_ft(wh))

    def _sort_observations(self, dt: datetime.datetime) -> List[Tuple[int, int, str]]:
        """Sort observations by area and elevation."""
        order = []

        for sid in self.obs[dt]:
            station = self.stations[sid]
            grid_cell = self.getGridCell(station.lat, station.lon)

            if grid_cell[0] is None or grid_cell[1] is None:
                continue

            x, y = int(grid_cell[0]), int(grid_cell[1])

            if self.dump_mask[y, x] < 0.5:
                continue

            # Find which sort area contains this station
            sort_val = 9999
            for i, mask in enumerate(self.sort_masks):
                if mask[y, x] > 0:
                    sort_val = i
                    break

            order.append((sort_val, 20000 - station.elev, sid))

        order.sort()
        return order

    def _write_output(self) -> None:
        """Write observation files."""
        os.makedirs(MOAC.Config["OBSDIR"], exist_ok=True)

        for dt in self.dt_list:
            sorted_stations = self._sort_observations(dt)
            if not sorted_stations:
                continue

            lines = []
            last_area = None

            for area_num, _, sid in sorted_stations:
                # Write area header if changed
                if area_num != last_area:
                    lines.append("#\n")
                    if area_num == 9999:
                        lines.append("#   Misc stations\n")
                        site_name = "MSC"
                    else:
                        lines.append(f"#  {self.sort_names[area_num]}\n")
                        site_name = self.sort_names[area_num]
                    lines.append("#\n")
                    last_area = area_num

                # Write observation line
                stn = self.stations[sid]
                obs = self.obs[dt][sid]

                try:
                    line = (
                        f"{sid:>5},"
                        f"{stn.name:40},"
                        f"{stn.lat:7.4f},"
                        f"{stn.lon:9.4f},"
                        f"{stn.elev:5d},"
                        f"{obs.format_value(obs.temp)},"
                        f"{obs.format_value(obs.dpt)},"
                        f"{obs.format_value(obs.wdir)},"
                        f"{obs.format_value(obs.wspd)},"
                        f"{obs.format_value(obs.effective_gust)},"
                        f"{obs.format_value(obs.skyamt)},"
                        f"{obs.format_value(obs.cig)},"
                        f"{obs.format_vis()},"
                        f"{obs.format_value(obs.wvehght)},"
                        f"{site_name}\n"
                    )
                    lines.append(line)
                except Exception as e:
                    logtime(f"Error formatting {sid}: {e}")

            # Write file
            outfile = os.path.join(MOAC.Config["OBSDIR"], f"{dt.strftime('%Y%m%d%H')}00.dat")
            with open(outfile, 'w') as f:
                f.writelines(lines)

            logtime(f"Wrote {len(lines)} lines to {outfile}")

    def _get_epoch(self, dt: datetime.datetime) -> float:
        """Convert datetime to Unix epoch."""
        dt0 = datetime.datetime(1970, 1, 1, tzinfo=GMT_ZONE)
        return (dt - dt0).total_seconds()
