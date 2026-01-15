#!/awips2/python/bin/python
# ----------------------------------------------------------------------------
# This software is in the public domain, furnished "as is", without technical
# support, and with no warranty, express or implied, as to its usefulness for
# any purpose.
#
# obsdump - part of MatchObsAll routines - version 0.215
#
#    Searches through METAR and mesonet netcdf files, and creates ASCII
#    comma-delimited files of Temp, Dewpoint, Wind Direction, Wind
#    speed, Wind Gust. If there are many obs in a single hour, the one
#    nearest the "top of the hour" is saved.
#
# Author: Tim Barker - SOO BOI
#   2025/01/15 - version 0.215. KA. Refactored for clarity and maintainability.
#                               Fixed station name lookup to fall back to file.
#   2021/04/21 - version 0.207. MBC. Added support for WaveHeight.
#   2020/12/18 - version 0.210. JRW. Updated to Python3 using updatePython
#   2018/06/04 - version 0.206. JRW. Added DUMPAREA to remove obs outside a
#                                    defined edit area.
#   2016/11/04 - version 0.205. JRW. Fixed Low Vis problems
#   2016/08/10 - version 0.204. JRW. Fixed WindGust/pkwnd computation
#   2016/08/10 - version 0.203. JRW. Removed need to import ObsdumpConfig.
#   2016/07/01 - version 0.202. JRW. Fixed bugs related to cloud heights.
#   2016/03/10 - version 0.117. JRW. Modified to run using the DAF.
#   2016/01/16 - version 0.116. JRW. Fixed Sky cover obs to use lowest coverage.
#   2015/12/22 - version 0.115. JRW. Modified Visibility to output in 100's.
#   2015/07/15 - version 0.114. JRW. Fixed Ceiling problem with missing heights.
#   2015/05/05 - version 0.112. JRW. Modified to use new Uengine Structure.
#   2014/12/05 - version 0.111. JRW. Added Aviation elements.
#   2012/01/16 - version 0.110. Change mesonet access for HDF5 files.
#   2011/07/20 - version 0.108. AWIPS-2 port.
#   2008/06/10 - version 0.105. Changed into its own procedure.
#   2008/05/27 - version 0.103. Fix issue with moving ships.
#   2008/04/11 - version 0.100. Added type checks for getGridCell.
#   2006/10/10 - version 0.99. Better handle character type of netCDF data.
#   2005/02/14 - version 0.98. Update version number.
#   2003/05/27 - version 0.95. Fix typo in getRoughLimits.
#   2003/04/08 - version 0.93. Non-polygon edit areas handled.
#   2003/02/26 - version 0.92. Version update only.
#   2003/02/07 - Add maritime data. Better log messages.
#   2003/02/06 - Remove commas from site ID or Name.
#   2003/01/20 - First version.
# ----------------------------------------------------------------------------
#
#  Normally don't want Obsdump to show up in menus - it gets run via cron.
#
MenuItems = None

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

# Default config values
CONFIG_DEFAULTS = {
    "DEBUG": 0,
    "PREVHOURS": 8,
    "OBSDIR": "/localapps/runtime/MatchObsAll/data/XXX",
    "STATIONFILE": "/awips/fxa/ldad/data/mesonetStation.txt",
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
    elev: int = 0

    def __post_init__(self):
        # Clamp elevation to valid range
        if self.elev > 20000:
            self.elev = 20000
        elif self.elev < -1000:
            logtime(f"Warning: elevation {self.elev} out of range, clamping to -1000")
            self.elev = -1000
        else:
            self.elev = int(self.elev)


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
            elif base is not None and base > -9000:
                self.sky_layer_base.append(int(base))
        elif cover:
            logtime(f"Unknown sky cover: {cover}", 10)

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

        for height, cover in layers:
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
    """Convert Celsius to Fahrenheit."""
    return (t * 9.0 / 5.0) + 32.0


def k_to_f(t: float) -> float:
    """Convert Kelvin to Fahrenheit."""
    return (t - 273.15) * 9.0 / 5.0 + 32.0


def ms_to_kt(speed: float) -> float:
    """Convert meters per second to knots."""
    return speed * 1.943


def m_to_ft(height: float) -> float:
    """Convert meters to feet."""
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

    logtime(f"Temperature out of range: {result}", 5)
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

    logtime(f"Wind speed out of range: {result}", 5)
    return None


# =============================================================================
# Main Procedure
# =============================================================================

class Procedure(SmartScript.SmartScript):
    def __init__(self, dbss):
        SmartScript.SmartScript.__init__(self, dbss)
        self.station_names: Dict[str, str] = {}
        self.file_station_names: Dict[str, str] = {}
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

    # -------------------------------------------------------------------------
    # Setup Methods
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # Station Name Loading
    # -------------------------------------------------------------------------

    def _load_station_names(self) -> None:
        """Load station names from database and station file."""
        self.station_names = {}
        self.file_station_names = {}

        # First, load from station file as fallback
        self._load_station_file()

        # Then load from DAL (preferred source)
        try:
            req = DAL.newDataRequest("common_obs_spatial")
            req.setParameters("stationid", "name")
            envelope = MultiPoint((
                (self.lon_max, self.lat_min),
                (self.lon_min, self.lat_max)
            ))
            req.setEnvelope(envelope)

            for geom in DAL.getGeometryData(req):
                sid = geom.getString("stationid")
                name = geom.getString("name")
                if sid not in self.station_names:
                    self.station_names[sid] = name

            logtime(f"Loaded {len(self.station_names)} station names from DAL")
        except Exception as e:
            logtime(f"Error loading station names from DAL: {e}")

        logtime(f"Loaded {len(self.file_station_names)} station names from file")

    def _load_station_file(self) -> None:
        """Load station names from mesonet station file."""
        stnfile = MOAC.Config.get("STATIONFILE", "/awips/fxa/ldad/data/mesonetStation.txt")

        if not os.path.exists(stnfile):
            logtime(f"Station file not found: {stnfile}")
            return

        try:
            with open(stnfile) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    data = line.split("|")
                    if len(data) < 3:
                        continue
                    sid = data[0].strip()
                    name = data[2].strip().replace(",", " ")
                    if sid and name:
                        self.file_station_names[sid] = name
        except IOError as e:
            logtime(f"Error reading station file: {e}")

    def _get_station_name(self, sid: str) -> str:
        """Get station name, falling back to file if DAL returned ID as name."""
        name = self.station_names.get(sid)

        # If DAL returned the ID as the name or no name found, try the file
        if name is None or name == sid:
            name = self.file_station_names.get(sid, sid)

        return name

    # -------------------------------------------------------------------------
    # Observation Fetching
    # -------------------------------------------------------------------------

    def _fetch_observations(self, source_name: str, config: dict) -> None:
        """Fetch observations from a single source type."""
        try:
            req = DAL.newDataRequest(config["request_type"])
            req.setParameters(*config["parameters"])
            envelope = MultiPoint((
                (self.lon_min, self.lat_min),
                (self.lon_max, self.lat_max)
            ))
            req.setEnvelope(envelope)
        except Exception as e:
            logtime(f"Error creating request for {source_name}: {e}")
            return

        delta_30min = datetime.timedelta(minutes=30)

        for dt in self.dt_list:
            begin_dt = dt - delta_30min
            end_dt = dt + delta_30min
            start = AbsTime.AbsTime(self._get_epoch(begin_dt))
            end = AbsTime.AbsTime(self._get_epoch(end_dt))
            tr = TimeRange.TimeRange(start, end)

            try:
                observations = DAL.getGeometryData(req, tr)
            except Exception as e:
                logtime(f"Error fetching {source_name} for {dt}: {e}")
                continue

            if not observations:
                logtime(f"WARNING: No {source_name} obs retrieved for {dt}")
                continue

            logtime(f"{len(observations)} {source_name} obs retrieved for {dt}", 10)

            for ob in observations:
                try:
                    self._process_observation(ob, dt, config)
                except Exception as e:
                    logtime(f"Error processing {source_name} ob: {e}", 5)

    def _process_observation(self, ob, dt: datetime.datetime, config: dict) -> None:
        """Process a single observation record."""
        sid = ob.getLocationName()

        # Create station data if needed
        if sid not in self.stations:
            name = self._get_station_name(sid)
            elev = 0
            if "elevation" in ob.getParameters():
                elev_val = ob.getNumber("elevation")
                if elev_val is not None and elev_val > -9000:
                    elev = elev_val * 3.2808

            self.stations[sid] = StationData(
                sid=sid,
                name=name,
                lat=ob.getGeometry().y,
                lon=ob.getGeometry().x,
                elev=int(min(20000, max(-1000, elev)))
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

        # Sky cover (these are separate records in the data)
        if 'skyCover' in params:
            cover = ob.getString('skyCover')
            base = ob.getNumber('skyLayerBase') if 'skyLayerBase' in params else -9999
            obs_record.add_sky_layer(cover, base)
            return  # Sky obs are separate records, don't process other fields

        # Temperature (prefer tenths if available)
        if 'tempFromTenths' in params:
            temp_val = ob.getNumber('tempFromTenths')
            if temp_val is not None and temp_val > -900:
                obs_record.temp = convert_temp(temp_val, temp_unit)
            elif 'temperature' in params:
                obs_record.temp = convert_temp(ob.getNumber('temperature'), temp_unit)
        elif 'temperature' in params:
            obs_record.temp = convert_temp(ob.getNumber('temperature'), temp_unit)

        # Dewpoint (prefer tenths if available)
        if 'dpFromTenths' in params:
            dp_val = ob.getNumber('dpFromTenths')
            if dp_val is not None and dp_val > -900:
                obs_record.dpt = convert_temp(dp_val, temp_unit)
            elif 'dewpoint' in params:
                obs_record.dpt = convert_temp(ob.getNumber('dewpoint'), temp_unit)
        elif 'dewpoint' in params:
            obs_record.dpt = convert_temp(ob.getNumber('dewpoint'), temp_unit)

        # Wind direction
        if 'windDir' in params:
            wdir = ob.getNumber('windDir')
            if wdir is not None and wdir > -9000 and 0 <= wdir <= 360:
                obs_record.wdir = int(wdir)

        # Wind speed
        if 'windSpeed' in params:
            obs_record.wspd = convert_wind(ob.getNumber('windSpeed'), wind_unit)

        # Wind gust
        if 'windGust' in params:
            obs_record.gust = convert_wind(ob.getNumber('windGust'), wind_unit)

        # Peak wind
        if 'pkwndSpeed' in params:
            obs_record.pkwd = convert_wind(ob.getNumber('pkwndSpeed'), wind_unit)

        # Visibility
        if 'visibility' in params:
            vis = ob.getNumber('visibility')
            if vis is not None and vis > -9000 and 0 <= vis <= 10:
                obs_record.vsby = vis

        # Wave height
        if 'waveHeight' in params:
            wh = ob.getNumber('waveHeight')
            if wh is not None and wh > -9000:
                obs_record.wvehght = int(m_to_ft(wh))

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------

    def _sort_observations(self, dt: datetime.datetime) -> List[Tuple[int, int, str]]:
        """Sort observations by area and elevation."""
        order = []

        for sid in self.obs[dt]:
            if sid not in self.stations:
                continue

            station = self.stations[sid]
            grid_cell = self.getGridCell(station.lat, station.lon)

            if grid_cell[0] is None or grid_cell[1] is None:
                continue

            x, y = int(grid_cell[0]), int(grid_cell[1])

            try:
                if self.dump_mask[y, x] < 0.5:
                    continue
            except IndexError:
                continue

            # Find which sort area contains this station
            sort_val = 9999
            for i, mask in enumerate(self.sort_masks):
                try:
                    if mask[y, x] > 0:
                        sort_val = i
                        break
                except IndexError:
                    continue

            order.append((sort_val, 20000 - station.elev, sid))

        order.sort()
        return order

    def _write_output(self) -> None:
        """Write observation files."""
        try:
            os.makedirs(MOAC.Config["OBSDIR"], exist_ok=True)
        except OSError as e:
            logtime(f"Error creating output directory: {e}")
            return

        for dt in self.dt_list:
            sorted_stations = self._sort_observations(dt)
            if not sorted_stations:
                logtime(f"No observations to write for {dt}")
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
                obs_rec = self.obs[dt][sid]

                try:
                    line = (
                        f"{sid:>5},"
                        f"{stn.name:40},"
                        f"{stn.lat:7.4f},"
                        f"{stn.lon:9.4f},"
                        f"{stn.elev:5d},"
                        f"{obs_rec.format_value(obs_rec.temp)},"
                        f"{obs_rec.format_value(obs_rec.dpt)},"
                        f"{obs_rec.format_value(obs_rec.wdir)},"
                        f"{obs_rec.format_value(obs_rec.wspd)},"
                        f"{obs_rec.format_value(obs_rec.effective_gust)},"
                        f"{obs_rec.format_value(obs_rec.skyamt)},"
                        f"{obs_rec.format_value(obs_rec.cig)},"
                        f"{obs_rec.format_vis()},"
                        f"{obs_rec.format_value(obs_rec.wvehght)},"
                        f"{site_name}\n"
                    )
                    lines.append(line)
                except Exception as e:
                    logtime(f"Error formatting {sid}: {e}")

            # Write file
            outfile = os.path.join(
                MOAC.Config["OBSDIR"],
                f"{dt.strftime('%Y%m%d%H')}00.dat"
            )

            try:
                with open(outfile, 'w') as f:
                    f.writelines(lines)
                logtime(f"Wrote {len(lines)} lines to {outfile}")
            except IOError as e:
                logtime(f"Error writing {outfile}: {e}")

    # -------------------------------------------------------------------------
    # Utility Methods
    # -------------------------------------------------------------------------

    def _get_epoch(self, dt: datetime.datetime) -> float:
        """Convert datetime to Unix epoch."""
        dt0 = datetime.datetime(1970, 1, 1, tzinfo=GMT_ZONE)
        return (dt - dt0).total_seconds()
