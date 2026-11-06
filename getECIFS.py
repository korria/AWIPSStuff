#!/awips2/python/bin/python
#
# Get ECMWF-IFS 0.1 deg open data and put into LDAD-GRIB
#
# Modeled on getECAI.py (ECMWF-AIFS). Pulls only the wanted fields out of the
# global 0.1 deg files using the .index byte offsets (~470 MB/hour otherwise),
# then regrids to a CONUS + NE Pacific 0.1 deg lat-lon box, south-to-north.
#
# v2026.11.05 - Korri Anderson - Initial version
# v2026.11.05 - Pull from the Google Cloud mirror first (AWS S3 bucket was returning
#               SlowDown), then S3, then data.ecmwf.int; back off when throttled;
#               fewer, larger range requests. No longer needs boto3.
# v2026.11.05 - 6-hourly 0-240 h (was 3-hourly to 144 h).
# v2026.11.05 - tp -> APCP (x1000, m -> mm), same as getECAI. AWIPS did not
#               decode the raw ECMWF local code 0-1-193.
# v2026.11.06 - Add surface z (f000 only) as HGT (m) terrain, same as getECAI.
#
import os
import json
import time
import random
import shutil
import logging
import datetime
import subprocess
import urllib.error
import urllib.request

# --- CONFIGURATION ---
# ECMWF open-data mirrors, tried in order. The AWS bucket throttles heavily
# (SlowDown), so it is a backup; ECMWF's own server is last, as ECMWF asks.
MIRRORS = [
    'https://storage.googleapis.com/ecmwf-open-data',
    'https://ecmwf-forecasts.s3.amazonaws.com',
    'https://data.ecmwf.int/forecasts',
]
STAGING_DIR = '/ldad/localapps/getECIFS/data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

MAX_LOOKBACK_DAYS = 3
CYCLES = ['12', '00']
# 6-hourly to 240 h, like AIFS (open data also has 3-hourly to 144 h).
F_HOURS = list(range(0, 241, 6))

# Example path: 20261105/00z/ifs/0p1/oper/20261105000000-24h-oper-fc-sfc.grib2
SRC_KEY_TEMPLATE = "{date}/{cycle}z/ifs/0p1/oper/{date}{cycle}0000-{fhour}h-oper-fc-{kind}"
DEST_FILENAME_TEMPLATE = "LDAD-GRIB-ECMWF-IFS_{}_t{}z_f{:03d}.grib2"

# Fields to keep. Surface/soil: the AIFS set plus IFS-only extras
# (10fg gust, mucape, ptype, 3 h and 6 h max/min temp; names a given hour
# doesn't have are simply skipped). Pressure levels: all levels.
SFC_PARAMS = {
    '100u', '100v', '10u', '10v', '2d', '2t', 'cp', 'hcc', 'lcc', 'mcc',
    'msl', 'sf', 'skt', 'sot', 'sp', 'ssrd', 'strd', 'tcc', 'tcw', 'tp',
    'vsw', 'z', '10fg', 'mucape', 'ptype', 'mx2t3', 'mn2t3', 'mx2t6', 'mn2t6',
}
PL_PARAMS = {'gh', 'q', 't', 'u', 'v', 'w'}

# ECMWF precip type codes -> NCEP-style categorical yes/no fields.
# 1 rain, 3 freezing rain, 5 snow, 6 wet snow, 7 rain/snow mix,
# 8 ice pellets, 11 drizzle, 12 freezing drizzle. A mix sets CRAIN and CSNOW.
PTYPE_TO_CATEGORICAL = [
    ('CRAIN', (1, 7, 11)),
    ('CSNOW', (5, 6, 7)),
    ('CFRZR', (3, 12)),
    ('CICEP', (8,)),
]

# Output grid (wgrib2 -new_grid latlon lon0:nlon:dlon lat0:nlat:dlat)
# 160W-60W, 20N-60N at 0.1 deg = 1001 x 401. AWIPS needs a matching grid.
OUT_LON = '200:1001:0.1'
OUT_LAT = '20:401:0.1'

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger(__name__)


def src_key(date_str, cycle, f_hour, kind):
    return SRC_KEY_TEMPLATE.format(date=date_str, cycle=cycle, fhour=f_hour, kind=kind)


def find_latest_complete_run():
    """Newest run whose last wanted hour is posted (both sfc and pl)."""
    now = datetime.datetime.now(datetime.timezone.utc)
    for i in range(MAX_LOOKBACK_DAYS):
        date_str = (now - datetime.timedelta(days=i)).strftime("%Y%m%d")
        for cycle in CYCLES:
            try:
                for kind in ('sfc.index', 'pl.index'):
                    fetch(src_key(date_str, cycle, F_HOURS[-1], kind))
                return date_str, cycle
            except Exception:
                log.info(f"{date_str} {cycle}z not complete yet")
    return None, None


def marker_path(date_str, cycle):
    return os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")


def wanted_ranges(index_text, params, pressure_levels):
    """
    Byte ranges (start, end inclusive) of the wanted fields from an ECMWF
    .index file (one JSON record per line), merged where adjacent.
    """
    ranges = []
    for line in index_text.splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec['param'] not in params or (rec['levtype'] == 'pl') != pressure_levels:
            continue
        start, end = rec['_offset'], rec['_offset'] + rec['_length'] - 1
        if ranges and ranges[-1][1] + 1 == start:
            ranges[-1] = (ranges[-1][0], end)
        else:
            ranges.append((start, end))
    return sorted(ranges)


def group_ranges(ranges, max_gap=5 * 1024 * 1024):
    """
    Group wanted ranges into fewer requests: ranges less than max_gap
    apart share one request. Returns [(span_start, span_end, [ranges])].
    """
    groups = []
    for start, end in ranges:
        if groups and start - groups[-1][1] - 1 <= max_gap:
            groups[-1] = (groups[-1][0], end, groups[-1][2] + [(start, end)])
        else:
            groups.append((start, end, [(start, end)]))
    return groups


def fetch(path, byte_range=None, max_retries=6):
    """
    GET a file (or inclusive byte range) from the open-data mirrors, in order.
    A mirror that throttles or errors passes to the next. A 404 means the
    file isn't posted yet. If every mirror fails, back off and try again.
    """
    headers = {'User-Agent': 'Mozilla/5.0'}
    if byte_range:
        headers['Range'] = 'bytes=%d-%d' % byte_range
    name = os.path.basename(path)
    for attempt in range(max_retries):
        codes = []
        for base in MIRRORS:
            try:
                req = urllib.request.Request(f'{base}/{path}', headers=headers)
                with urllib.request.urlopen(req, timeout=120) as resp:
                    data = resp.read()
                if byte_range and len(data) != byte_range[1] - byte_range[0] + 1:
                    raise IOError(f'got {len(data)} bytes for range {byte_range}')
                return data
            except urllib.error.HTTPError as e:
                codes.append(e.code)
                if e.code != 404:
                    log.warning(f"HTTP {e.code} from {base} for {name}")
            except (urllib.error.URLError, OSError) as e:
                codes.append(None)
                log.warning(f"{base} failed for {name}: {e}")
        if 404 in codes:  # a mirror says not posted; next cron run retries
            raise IOError(f"{name} not posted yet")
        wait = min(5 * 2 ** attempt, 120) + random.uniform(0, 1)
        log.warning(f"All mirrors failed for {name}; waiting {wait:.0f}s "
                    f"(attempt {attempt + 1}/{max_retries})")
        time.sleep(wait)
    raise IOError(f"Could not download {name} after {max_retries} attempts")


def download_fields(date_str, cycle, f_hour, out_path):
    """Write the wanted GRIB messages from the sfc and pl files to out_path."""
    with open(out_path, 'wb') as out:
        for kind, params, is_pl in (('sfc', SFC_PARAMS, False), ('pl', PL_PARAMS, True)):
            base = src_key(date_str, cycle, f_hour, kind)
            index = fetch(base + '.index').decode()
            for span_start, span_end, pieces in group_ranges(wanted_ranges(index, params, is_pl)):
                data = fetch(base + '.grib2', (span_start, span_end))
                for start, end in pieces:
                    out.write(data[start - span_start:end - span_start + 1])


def sanitize_grib(input_path, output_path):
    """
    Step 1: regrid to the output box (also flips ECMWF north-to-south order
            to south-to-north). Nearest neighbour: the box sits on the source
            0.1 deg points, so values are copied, and ptype stays categorical.
    Step 2: tp (ECMWF local 0-1-193, metres) -> APCP (kg m-2 = mm), as in
            getECAI; surface geopotential z (f000 only) -> HGT terrain (m);
            append CRAIN/CSNOW/CFRZR/CICEP derived from PTYPE.
    """
    temp_path = input_path + ".box"
    cmd_step1 = [
        WGRIB2_BIN, input_path,
        '-set_grib_type', 'c3',
        '-new_grid_winds', 'earth',
        '-new_grid_interpolation', 'neighbor',
        # Source winds are already earth-relative on the same points, and
        # ECMWF does not keep U and V adjacent (wgrib2 drops unpaired winds).
        '-new_grid_vectors', 'none',
        '-new_grid', 'latlon', OUT_LON, OUT_LAT,
        temp_path,
    ]
    cmd_step2 = [
        WGRIB2_BIN, temp_path,
        # wgrib2 names 0-1-193 tprate (newer) or var0_1_193 (older builds).
        '-if', ':(tprate|var0_1_193|var0_1_1_98_1_193):',
            '-rpn', '1000:*',
            '-set_var', 'APCP',
        '-endif',
        '-if', ':GP:',
            '-rpn', '9.80665:/',
            '-set_var', 'HGT',
        '-endif',
        '-set_grib_type', 'c3',
        '-grib_out', output_path,
    ]
    # Keep PTYPE (written above) and append NCEP-style 0/1 categorical fields.
    cmd_step2 += ['-if', ':PTYPE:', '-rpn', 'sto_1']
    for var, codes in PTYPE_TO_CATEGORICAL:
        test = ':'.join(f'rcl_1:{c}:==' + (':+' if i else '') for i, c in enumerate(codes))
        cmd_step2 += ['-rpn', test, '-set_var', var, '-grib_out', output_path]
    cmd_step2 += ['-endif']
    try:
        subprocess.check_output(cmd_step1, stderr=subprocess.STDOUT)
        subprocess.check_output(cmd_step2, stderr=subprocess.STDOUT)
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"wgrib2 failed: {e.output.decode('utf-8', 'replace')}")
        return False
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def process_run(date_str, cycle):
    log.info(f">>> PROCESSING: {date_str} {cycle}z <<<")
    os.makedirs(STAGING_DIR, exist_ok=True)
    os.makedirs(FINAL_DIR, exist_ok=True)

    for f_hour in F_HOURS:
        dest_filename = DEST_FILENAME_TEMPLATE.format(date_str, cycle, f_hour)
        final_path = os.path.join(FINAL_DIR, dest_filename)
        temp_move_path = os.path.join(FINAL_DIR, f".tmp_{dest_filename}")
        raw_path = os.path.join(STAGING_DIR, f"raw_{dest_filename}")
        clean_path = os.path.join(STAGING_DIR, dest_filename)

        if os.path.exists(final_path):
            log.info(f"[EXISTS] {dest_filename}")
            continue

        try:
            download_fields(date_str, cycle, f_hour, raw_path)
            if not sanitize_grib(raw_path, clean_path):
                log.error(f"Conversion failed for f{f_hour:03d}")
                return False
            shutil.move(clean_path, temp_move_path)
            os.rename(temp_move_path, final_path)
            log.info(f"Processed f{f_hour:03d} -> {final_path}")
        except Exception as e:
            log.error(f"Failed on f{f_hour:03d}: {e}")
            return False
        finally:
            for path in (raw_path, clean_path):
                if os.path.exists(path):
                    os.remove(path)

    log.info(f"Success! All files for {date_str} {cycle}z moved to {FINAL_DIR}")
    return True


def main():
    date_str, cycle = find_latest_complete_run()
    if not date_str:
        log.error(f"No complete runs found in the last {MAX_LOOKBACK_DAYS} days.")
        return

    if os.path.exists(marker_path(date_str, cycle)):
        log.info(f"Run {date_str} {cycle}z already finished. Exiting.")
        return

    if process_run(date_str, cycle):
        with open(marker_path(date_str, cycle), 'w') as f:
            f.write(f"Download complete at {datetime.datetime.now(datetime.timezone.utc)}")
        log.info(f"Run marked complete: {marker_path(date_str, cycle)}")


if __name__ == "__main__":
    main()
