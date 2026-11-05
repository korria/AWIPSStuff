#!/awips2/python/bin/python
#
# Get ECMWF-IFS 0.1 deg open data and put into LDAD-GRIB
#
# Modeled on getECAI.py (ECMWF-AIFS). Pulls only the wanted fields out of the
# global 0.1 deg files using the .index byte offsets (~470 MB/hour otherwise),
# then regrids to a CONUS + NE Pacific 0.1 deg lat-lon box, south-to-north.
#
# v2026.11.05 - Korri Anderson - Initial version
#
import os
import sys
sys.path.insert(0, "/ldad/localapps/getRRFS/bin/")
import json
import shutil
import logging
import datetime
import subprocess
import boto3
from botocore import UNSIGNED
from botocore.config import Config

# --- CONFIGURATION ---
BUCKET_NAME = 'ecmwf-forecasts'
STAGING_DIR = '/ldad/localapps/getECIFS/data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

MAX_LOOKBACK_DAYS = 3
CYCLES = ['12', '00']
# Open data steps: 3-hourly to 144 h, then 6-hourly. Capped at 240 h like AIFS.
F_HOURS = list(range(0, 145, 3)) + list(range(150, 241, 6))

# Example key: 20261105/00z/ifs/0p1/oper/20261105000000-24h-oper-fc-sfc.grib2
SRC_KEY_TEMPLATE = "{date}/{cycle}z/ifs/0p1/oper/{date}{cycle}0000-{fhour}h-oper-fc-{kind}"
DEST_FILENAME_TEMPLATE = "LDAD-GRIB-ECMWF-IFS_{}_t{}z_f{:03d}.grib2"

# Fields to keep. Surface/soil: the AIFS set plus IFS-only extras
# (10fg gust, mucape, ptype, 3 h max/min temp). Pressure levels: all levels.
SFC_PARAMS = {
    '100u', '100v', '10u', '10v', '2d', '2t', 'cp', 'hcc', 'lcc', 'mcc',
    'msl', 'sf', 'skt', 'sot', 'sp', 'ssrd', 'strd', 'tcc', 'tcw', 'tp',
    'vsw', '10fg', 'mucape', 'ptype', 'mx2t3', 'mn2t3',
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


def get_s3_client():
    return boto3.client('s3', config=Config(
        signature_version=UNSIGNED,
        connect_timeout=10,
        read_timeout=120,
        retries={'max_attempts': 8, 'mode': 'adaptive'}
    ))


def src_key(date_str, cycle, f_hour, kind):
    return SRC_KEY_TEMPLATE.format(date=date_str, cycle=cycle, fhour=f_hour, kind=kind)


def find_latest_complete_run(client):
    """Newest run whose last wanted hour is posted (both sfc and pl)."""
    now = datetime.datetime.now(datetime.timezone.utc)
    for i in range(MAX_LOOKBACK_DAYS):
        date_str = (now - datetime.timedelta(days=i)).strftime("%Y%m%d")
        for cycle in CYCLES:
            try:
                for kind in ('sfc.index', 'pl.index'):
                    client.head_object(Bucket=BUCKET_NAME,
                                       Key=src_key(date_str, cycle, F_HOURS[-1], kind))
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


def download_fields(client, date_str, cycle, f_hour, out_path):
    """Write the wanted GRIB messages from the sfc and pl files to out_path."""
    with open(out_path, 'wb') as out:
        for kind, params, is_pl in (('sfc', SFC_PARAMS, False), ('pl', PL_PARAMS, True)):
            base = src_key(date_str, cycle, f_hour, kind)
            index = client.get_object(Bucket=BUCKET_NAME, Key=base + '.index')['Body'].read()
            for start, end in wanted_ranges(index.decode(), params, is_pl):
                body = client.get_object(Bucket=BUCKET_NAME, Key=base + '.grib2',
                                         Range=f'bytes={start}-{end}')['Body']
                shutil.copyfileobj(body, out, 1024 * 1024)


def sanitize_grib(input_path, output_path):
    """
    Step 1: regrid to the output box (also flips ECMWF north-to-south order
            to south-to-north). Nearest neighbour: the box sits on the source
            0.1 deg points, so values are copied, and ptype stays categorical.
    Step 2: tp (metres) -> APCP (kg m-2 = mm), as AWIPS expects, and
            CRAIN/CSNOW/CFRZR/CICEP derived from PTYPE.
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
        '-if', ':(tprate|TPRATE|TP|var0_1_228):',
            '-rpn', '1000:*',
            '-set_var', 'APCP',
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


def process_run(client, date_str, cycle):
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
            download_fields(client, date_str, cycle, f_hour, raw_path)
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
    client = get_s3_client()

    date_str, cycle = find_latest_complete_run(client)
    if not date_str:
        log.error(f"No complete runs found in the last {MAX_LOOKBACK_DAYS} days.")
        return

    if os.path.exists(marker_path(date_str, cycle)):
        log.info(f"Run {date_str} {cycle}z already finished. Exiting.")
        return

    if process_run(client, date_str, cycle):
        with open(marker_path(date_str, cycle), 'w') as f:
            f.write(f"Download complete at {datetime.datetime.now(datetime.timezone.utc)}")
        log.info(f"Run marked complete: {marker_path(date_str, cycle)}")


if __name__ == "__main__":
    main()
