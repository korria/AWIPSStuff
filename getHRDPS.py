#!/awips2/python/bin/python
"""
getHRDPS.py

Download the Canadian HRDPS 2.5 km continental run from the MSC Datamart,
regrid it to the BOI polar-stereographic grid, and put one GRIB2 file per
forecast hour in /data/Incoming for LDAD-GRIB ingest.

Run from cron; it is safe to run as often as every 15 minutes. Each run picks
the newest cycle whose last hour is fully posted and delivers the hours not yet
delivered. A failed hour stops the run, which exits 1, and the next run resumes
from that hour.

HISTORY:
2026/11/07 - version 2.0. KA. Rewrite. Yesterday's 18Z is read from its own
             date directory (under /today/ it 404'd every night 00-04Z). Pick
             the newest fully posted run instead of guessing from the clock.
             PTYPE, WEonG precip type and WDIR are regridded nearest-neighbour
             (bilinear made false freezing rain on rain/snow edges and
             averaged 350/10 deg to 180). Hours are recorded as delivered, so a
             run is not downloaded and ingested again; incomplete hours are
             not ingested. Retries on 5xx, timeouts and truncated downloads.
             Atomic move into /data/Incoming, lock against overlapping runs,
             4 parallel downloads. Terrain (_HGT_Sfc) at f000 only. Dropped
             unused pygrib/numpy/xarray and requests. CRAIN/CSNOW/CFRZR/CICEP
             made from PTYPE, as HRRR has them in AWIPS (for ptyp4). Fog
             visibility (lower of WEonG ice/liquid fog) added as Vis SFC.
             Dropped
             WEonG PCPNTYPE (same PTyp slot as PTYPE) and WEonG
             DMNTPCPNTYPE/CONDAPCPN (AWIPS doesn't decode them).
2024/05/24 - version 1.0. KA. Initial version.
"""
import os
import re
import sys
import time
import fcntl
import random
import logging
import datetime
import subprocess
import http.client
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

# --- CONFIGURATION ---
DIR_URL_TEMPLATE = ('https://dd.weather.gc.ca/{date}/WXO-DD/model_hrdps/continental/'
                    '2.5km/{cycle}/{fhour:03d}/')
SRC_FILENAME_TEMPLATE = '{date}T{cycle}Z_MSC_HRDPS{var}_RLatLon0.0225_PT{fhour:03d}H.grib2'
DEST_FILENAME_TEMPLATE = 'LDAD-GRIB-BOI-HRDPS_{date}{cycle}_P{fhour:03d}.grib2'
STAGING_DIR = '/ldad/localapps/getHRDPS/hrdps_data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

CYCLES = ['18', '12', '06', '00']
F_HOURS = range(0, 49)
DOWNLOAD_WORKERS = 4        # the Datamart is a shared server; keep this modest
STATE_MAX_AGE_DAYS = 3

# BOI 2.5 km polar-stereographic grid (wgrib2 -new_grid nps ...) and the
# process id AWIPS uses to recognise the model. Ingest depends on both.
OUT_GRID = ['nps:-117.0:60.0', '-130.0:750:2500', '39.0:600:2500']
PROCESS_ID = '50'

VARIABLES = [
 "_TMP_AGL-2m",
 #"_DPT_AGL-2m",
 #"_DEPR_AGL-2m",
 #"_RH_AGL-2m", 
 "_HPBL_Sfc",
 "_HGT_Sfc",
 "_HGT_ISBL_0500",
 "_HGT_ISBL_0700",
 "_HGT_ISBL_0750",
 "_HGT_ISBL_0800",
 "_HGT_ISBL_0850",
 "_HGT_ISBL_0875",
 "_HGT_ISBL_0900",
 "_HGT_ISBL_0925",
 #"_HGT_ISBL_0950",
 "_HGT_ISBL_1000",
 "_TMP_ISBL_0500",
 "_TMP_ISBL_0700",
 "_TMP_ISBL_0750",
 "_TMP_ISBL_0800",
 "_TMP_ISBL_0850",
 "_TMP_ISBL_0875",
 "_TMP_ISBL_0900",
 "_TMP_ISBL_0925",
 #"_TMP_ISBL_0950",
 "_TMP_ISBL_1000",
 "_RH_ISBL_0500",
 "_RH_ISBL_0700",
 "_RH_ISBL_0750",
 "_RH_ISBL_0800",
 "_RH_ISBL_0850",
 "_RH_ISBL_0875",
 "_RH_ISBL_0900",
 "_RH_ISBL_0925",
 #"_RH_ISBL_0950",
 "_RH_ISBL_1000",
 #"_HGT_ISBY-1000-500",
 #"_GUST_AGL-10m",
 #"_GUST-Max_AGL-10m",
 #"_VGRD_AGL-10m",
 #"_UGRD_AGL-10m",
 "_WIND_AGL-10m",
 "_WDIR_AGL-10m",
 "_APCP-Accum1h_Sfc",
 "_APCP_Sfc",
 "_TCDC_Sfc",
 "_SNOD_Sfc",
 "_PRMSL_MSL",
 "_PRATE_Sfc",
 "_PTYPE_Sfc",
 #"_WEAFR_Sfc",
 #"_WEAPE_Sfc",
 #"_WEARN_Sfc",
 #"_WEASN_Sfc",
 "-WEonG_DPT_Sfc",
 "-WEonG_GUST_Sfc",
 #"-WEonG_HGTSNLVL_Sfc",
 #"-WEonG_CONDAPCPN_Sfc",      # AWIPS has no table entry (0-1-159): decodes as "Missing"
 #"-WEonG_PCPNTYPE_Sfc",       # same code as _PTYPE_Sfc (PTyp SFC): AWIPS keeps only one
 #"-WEonG_DMNTPCPNTYPE_Sfc",   # AWIPS has no table entry (0-1-123): decodes as "Missing"
 "-WEonG_PROBPCPN_Sfc",
 "-WEonG_VISIFG_Sfc",          # ice and liquid fog visibility are delivered
 "-WEonG_VISLFG_Sfc",          # as one Vis SFC field (see FOG_VISIBILITY)
 #"-WEonG_PROBLPCPN_Sfc",
 #"-WEonG_PROBFZPCPN_Sfc",
 #"-WEonG_PROBSN_Sfc",
 #"-WEonG_PROBSNSQ_Sfc",
 #"-WEonG_PROBFZRA_Sfc",
 #"-WEonG_PROBRA_Sfc",
 #"-WEonG_PROBTSOCRNC_Sfc",
 #"-WEonG_PROBPL_Sfc",
 #"-WEonG_PROBDZ_Sfc",
 #"-WEonG_PROBFZDZ_Sfc",
 #"_LAND_Sfc",
 
 ]

# Category codes and directions can't be averaged (rain 1 next to snow 5
# gives 3, freezing rain; 350 and 10 deg give 180), so these are regridded
# nearest-neighbour. Everything else is bilinear.
NEAREST_NEIGHBOUR = {
    '_PTYPE_Sfc',
    '-WEonG_PCPNTYPE_Sfc',
    '-WEonG_DMNTPCPNTYPE_Sfc',
    '_WDIR_AGL-10m',
}

# Fields that don't change with forecast hour; downloaded at f000 only.
F000_ONLY = {'_HGT_Sfc'}   # terrain

# HRRR has no PTYPE in AWIPS, only 0/1 CRAIN/CSNOW/CFRZR/CICEP, which the
# ptyp4 displays use. These are made from PTYPE (WMO code table 4.201) the
# same way getECIFS does: 1 rain, 3 freezing rain, 5 snow, 6 wet snow,
# 7 rain/snow mix, 8 ice pellets, 11 drizzle, 12 freezing drizzle. A mix sets
# CRAIN and CSNOW. HRDPS sets PTYPE everywhere, raining or not, so these are
# the type precipitation would be; ptyp4 multiplies them by TP1hr.
CATEGORICAL_SOURCE = '_PTYPE_Sfc'
PTYPE_TO_CATEGORICAL = [
    ('CRAIN', (1, 7, 11)),
    ('CSNOW', (5, 6, 7)),
    ('CFRZR', (3, 12)),
    ('CICEP', (8,)),
]

# HRDPS has no total visibility, only WEonG fog visibility (m, capped at
# 10 km; precipitation doesn't lower it). AWIPS doesn't decode those codes
# (0-19-34 ice fog, 0-19-33 liquid fog), so the lower of the two is
# delivered as standard visibility at the surface (Vis SFC), as HRRR's is.
FOG_VISIBILITY = ('-WEonG_VISIFG_Sfc', '-WEonG_VISLFG_Sfc')   # ice, liquid

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger(__name__)


def src_filename(date, cycle, fhour, var):
    return SRC_FILENAME_TEMPLATE.format(date=date, cycle=cycle, fhour=fhour, var=var)


def hour_variables(fhour):
    """Fields to download for this forecast hour."""
    return [v for v in VARIABLES if fhour == 0 or v not in F000_ONLY]


def count_grib2(data):
    """Count the GRIB2 messages in data; raise IOError if it is not whole GRIB2."""
    n = pos = 0
    while pos < len(data):
        size = int.from_bytes(data[pos + 8:pos + 16], 'big')
        if size < 16 or data[pos:pos + 4] != b'GRIB' or data[pos + size - 4:pos + size] != b'7777':
            raise IOError('incomplete or corrupt GRIB2 data')
        pos += size
        n += 1
    return n


def fetch(url, max_retries=5):
    """
    GET url. Returns None on 404 (not posted). Other HTTP errors, timeouts and
    truncated or corrupt GRIB2 are retried with backoff, then raise IOError.
    """
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = resp.read()
            if url.endswith('.grib2'):
                count_grib2(data)
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            err = f'HTTP {e.code}'
        except (OSError, http.client.HTTPException) as e:
            err = e
        if attempt < max_retries:
            wait = 5 * 2 ** (attempt - 1) + random.uniform(0, 1)
            log.warning(f'{os.path.basename(url)}: {err}; retrying in {wait:.0f}s')
            time.sleep(wait)
    raise IOError(f'{os.path.basename(url)}: {err} after {max_retries} attempts')


def find_latest_complete_run():
    """Newest cycle (today, then yesterday) whose last hour has every wanted field."""
    now = datetime.datetime.now(datetime.timezone.utc)
    last = F_HOURS[-1]
    for day in (now, now - datetime.timedelta(days=1)):
        date = day.strftime('%Y%m%d')
        for cycle in CYCLES:
            listing = fetch(DIR_URL_TEMPLATE.format(date=date, cycle=cycle, fhour=last))
            if listing is None:
                continue
            posted = set(re.findall(r'href="([^"]+\.grib2)"', listing.decode('utf-8', 'replace')))
            missing = [v for v in hour_variables(last) if src_filename(date, cycle, last, v) not in posted]
            if not missing:
                return date, cycle
            log.info(f'{date} {cycle}Z not complete yet ({len(missing)} fields missing at f{last:03d})')
    return None, None


def delivered_path(date, cycle):
    return os.path.join(STAGING_DIR, f'delivered_{date}_{cycle}z')


def read_delivered(date, cycle):
    """Forecast hours of this run already moved into FINAL_DIR."""
    try:
        with open(delivered_path(date, cycle)) as f:
            return {int(line) for line in f if line.strip()}
    except FileNotFoundError:
        return set()


def wgrib2(data, args):
    """Run wgrib2 on GRIB2 bytes; 'OUT' in args is the output file. Returns its contents."""
    in_path = os.path.join(STAGING_DIR, 'wgrib2_in.grib2')
    out_path = os.path.join(STAGING_DIR, 'wgrib2_out.grib2')
    try:
        with open(in_path, 'wb') as f:
            f.write(data)
        subprocess.check_output([WGRIB2_BIN, in_path] + [out_path if a == 'OUT' else a for a in args],
                                stderr=subprocess.STDOUT)
        with open(out_path, 'rb') as f:
            return f.read()
    finally:
        for path in (in_path, out_path):
            if os.path.exists(path):
                os.remove(path)


def regrid(data, interpolation):
    """Regrid GRIB2 messages to OUT_GRID with wgrib2; return the result."""
    # Levels are left as MSC sends them (gust at 10 m, cloud cover at the
    # surface), which is what the BOI displays are set up for.
    out = wgrib2(data, ['-set', 'analysis_or_forecast_process_id', PROCESS_ID,
                        '-new_grid_winds', 'earth',
                        '-new_grid_interpolation', interpolation,
                        '-new_grid'] + OUT_GRID + ['OUT'])
    if count_grib2(out) != count_grib2(data):
        raise IOError(f'wgrib2 wrote {count_grib2(out)} of {count_grib2(data)} fields ({interpolation})')
    return out


def categorical(ptype):
    """CRAIN/CSNOW/CFRZR/CICEP (0/1) from one PTYPE message, on the same grid."""
    args = ['-set_grib_type', 'c3', '-rpn', 'sto_1']
    for var, codes in PTYPE_TO_CATEGORICAL:
        test = ':'.join(f'rcl_1:{c}:==' + (':+' if i else '') for i, c in enumerate(codes))
        args += ['-rpn', test, '-set_var', var, '-grib_out', 'OUT']
    out = wgrib2(ptype, args)
    if count_grib2(out) != len(PTYPE_TO_CATEGORICAL):
        raise IOError(f'wgrib2 wrote {count_grib2(out)} of {len(PTYPE_TO_CATEGORICAL)} categorical fields')
    return out


def fog_visibility(ice, liquid):
    """The lower of ice- and liquid-fog visibility, as VIS at the surface."""
    # Records are picked by position, not name: older wgrib2 builds don't
    # know these codes. The result keeps the liquid-fog record's metadata.
    out = wgrib2(ice + liquid, ['-if', '^1:', '-rpn', 'sto_1', '-endif',
                                '-if', '^2:', '-rpn', 'rcl_1:min', '-set_var', 'VIS',
                                '-set_lev', 'surface', '-set_grib_type', 'c3',
                                '-grib_out', 'OUT', '-endif'])
    if count_grib2(out) != 1:
        raise IOError(f'wgrib2 wrote {count_grib2(out)} of 1 fog visibility field')
    return out


def process_hour(pool, date, cycle, fhour):
    """Download, regrid and deliver one forecast hour. Returns the field count."""
    dir_url = DIR_URL_TEMPLATE.format(date=date, cycle=cycle, fhour=fhour)
    variables = hour_variables(fhour)
    gribs = pool.map(lambda v: fetch(dir_url + src_filename(date, cycle, fhour, v)), variables)

    groups = {'bilinear': [], 'neighbor': []}
    missing = []
    fog = {}
    for var, grib in zip(variables, gribs):
        if grib is None:
            missing.append(var)
        elif var in FOG_VISIBILITY:
            fog[var] = grib
        else:
            groups['neighbor' if var in NEAREST_NEIGHBOUR else 'bilinear'].append(grib)
        if grib and var == CATEGORICAL_SOURCE:
            groups['neighbor'].append(categorical(grib))
    if missing and fhour > 0:
        raise IOError(f'not posted: {" ".join(missing)}')
    if missing:  # f000 has no accumulations, PTYPE or WEonG fields
        log.info(f'f000 has no {" ".join(missing)} (normal)')
    if len(fog) == len(FOG_VISIBILITY):
        groups['bilinear'].append(fog_visibility(*(fog[v] for v in FOG_VISIBILITY)))

    dest = DEST_FILENAME_TEMPLATE.format(date=date, cycle=cycle, fhour=fhour)
    tmp_path = os.path.join(FINAL_DIR, '.tmp_' + dest)
    data = b''.join(regrid(b''.join(messages), interpolation)
                    for interpolation, messages in groups.items() if messages)
    try:
        with open(tmp_path, 'wb') as out:
            out.write(data)
        os.rename(tmp_path, os.path.join(FINAL_DIR, dest))  # atomic: ingest never sees a partial file
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return count_grib2(data)


def process_run(date, cycle):
    """Deliver every hour not yet delivered. Stops at the first failure; True if all are in."""
    delivered = read_delivered(date, cycle)
    with ThreadPoolExecutor(DOWNLOAD_WORKERS) as pool:
        for fhour in F_HOURS:
            if fhour in delivered:
                continue
            start = time.time()
            try:
                n = process_hour(pool, date, cycle, fhour)
            except subprocess.CalledProcessError as e:
                log.error(f'f{fhour:03d}: wgrib2 failed: {e.output.decode("utf-8", "replace")}')
                return False
            except Exception as e:
                log.error(f'f{fhour:03d}: {e}')
                return False
            with open(delivered_path(date, cycle), 'a') as f:
                f.write(f'{fhour}\n')
            log.info(f'f{fhour:03d}: {n} fields delivered in {time.time() - start:.0f}s')
    return True


def clean_old_state():
    """Remove delivered_* records older than STATE_MAX_AGE_DAYS."""
    cutoff = time.time() - STATE_MAX_AGE_DAYS * 86400
    for name in os.listdir(STAGING_DIR):
        path = os.path.join(STAGING_DIR, name)
        if name.startswith('delivered_') and os.path.getmtime(path) < cutoff:
            os.remove(path)


def main():
    os.makedirs(STAGING_DIR, exist_ok=True)
    lock = open(os.path.join(STAGING_DIR, '.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        log.info('Previous run still active. Exiting.')
        return 0
    clean_old_state()

    try:
        date, cycle = find_latest_complete_run()
    except IOError as e:
        log.error(f'Could not check the Datamart: {e}')
        return 1
    if not date:
        log.error('No complete HRDPS run posted today or yesterday.')
        return 1
    if read_delivered(date, cycle) >= set(F_HOURS):
        log.info(f'Run {date} {cycle}Z already delivered. Exiting.')
        return 0

    log.info(f'>>> PROCESSING: {date} {cycle}Z <<<')
    if not process_run(date, cycle):
        log.error(f'Run {date} {cycle}Z incomplete; the next run resumes it.')
        return 1
    log.info(f'Success! All {len(F_HOURS)} hours of {date} {cycle}Z are in {FINAL_DIR}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
