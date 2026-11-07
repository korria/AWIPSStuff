#!/awips2/python/bin/python
"""
getRRFS.py

Download the operational RRFS 3 km CONUS 2D-field files from the NOAA open
data bucket on AWS, clip them to the western US with wgrib2, and put one
GRIB2 file per forecast hour in /data/Incoming for LDAD-GRIB ingest. RRFS
smoke (organic matter MASSDEN at 8 m and COLMD) is relabelled so AWIPS
decodes it.

Run from cron; it is safe to run as often as every 15 minutes. Each run picks
the newest 00/06/12/18Z cycle with every hour posted and delivers the hours
not yet delivered. A failed hour stops the run, which exits 1, and the next
run resumes from that hour.

HISTORY:
2026/11/07 - version 2.0. KA. Wait until every hour (f001-f084) is posted;
             the f072 check let runs start while f073-f084 were still posting
             (late hours post out of order). Look only in rrfs. date folders
             (it checked 240 rrfsens. paths first on every run). Delivered
             hours are recorded, so a failure resumes instead of downloading
             and ingesting the run again. Lock against overlapping runs,
             logging, exit 1 on failure; S3 errors no longer reported as "no
             run found". Clipped output written with complex packing (c3), as
             RRFS sends it: same values, about half the size.
"""
import os
import sys
import time
import fcntl
import shutil
import logging
import subprocess
sys.path.insert(0, "/ldad/localapps/getRRFS/bin/")
import boto3
from botocore import UNSIGNED
from botocore.client import Config

# --- CONFIGURATION ---
BUCKET_NAME = 'noaa-rrfs-ops-pds'
STAGING_DIR = '/ldad/localapps/getRRFS/data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

CYCLES = ['18', '12', '06', '00']   # the cycles that run to 84 h, newest first
F_HOURS = range(1, 85)              # hourly f001-f084
MAX_LOOKBACK_DAYS = 2
STATE_MAX_AGE_DAYS = 3              # must be more than MAX_LOOKBACK_DAYS

# Western US Clipping Bounds (Longitude: West:East, Latitude: South:North)
# This covers WA, OR, CA, NV, ID, UT, AZ and parts of MT, WY, CO, NM
WEST_BOUNDS = "-126:-101 31:50"

# Source: rrfs.YYYYMMDD/HH/rrfs.tHHz.2dfld.3km.fFFF.conus.grib2
SRC_FILENAME_TEMPLATE = "rrfs.t{}z.2dfld.3km.f{:03d}.conus.grib2"
DEST_FILENAME_TEMPLATE = "LDAD-GRIB-LocalRRFS_West_{}_t{}z_f{:03d}.grib2"

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger(__name__)


def get_s3_client():
    return boto3.client('s3', config=Config(signature_version=UNSIGNED))


def src_key(date_str, cycle, f_hour):
    return f"rrfs.{date_str}/{cycle}/" + SRC_FILENAME_TEMPLATE.format(cycle, f_hour)


def list_bucket(client, prefix, delimiter=''):
    """All keys (or, with a delimiter, folder prefixes) under prefix."""
    found = set()
    for page in client.get_paginator('list_objects_v2').paginate(
            Bucket=BUCKET_NAME, Prefix=prefix, Delimiter=delimiter):
        found.update(obj['Key'] for obj in page.get('Contents', []))
        found.update(p['Prefix'] for p in page.get('CommonPrefixes', []))
    return found


def find_latest_complete_run(client):
    """Newest (date, cycle) with every hour in F_HOURS posted."""
    # Only rrfs.YYYYMMDD/ folders; the bucket also holds rrfsens., refs. and firewx.
    folders = sorted(list_bucket(client, 'rrfs.', '/'), reverse=True)
    for folder in folders[:MAX_LOOKBACK_DAYS]:
        date_str = folder[len('rrfs.'):-1]
        for cycle in CYCLES:
            posted = list_bucket(client, f"rrfs.{date_str}/{cycle}/rrfs.t{cycle}z.2dfld.3km.f")
            missing = [h for h in F_HOURS if src_key(date_str, cycle, h) not in posted]
            if not missing:
                return date_str, cycle
            if len(missing) < len(F_HOURS):
                log.info(f"{date_str} {cycle}z not complete yet ({len(missing)} hours missing)")
    return None, None


def marker_path(date_str, cycle):
    return os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")


def delivered_path(date_str, cycle):
    return os.path.join(STAGING_DIR, f"delivered_{date_str}_{cycle}z")


def read_delivered(date_str, cycle):
    """Forecast hours of this run already moved into FINAL_DIR."""
    try:
        with open(delivered_path(date_str, cycle)) as f:
            return {int(line) for line in f if line.strip()}
    except FileNotFoundError:
        return set()


def clip_grib_west(input_path, output_path, f_hour):
    """
    Uses wgrib2 -small_grib to clip data to Western US coordinates.
    """
    # -small_grib lon_min:lon_max lat_min:lat_max output_file
    temp_met = input_path + ".met"
    temp_smoke_raw = input_path + ".smoke_raw"
    temp_smoke_final = input_path + ".smoke_final"

    lon_bounds = WEST_BOUNDS.split()[0]
    lat_bounds = WEST_BOUNDS.split()[1]

    cmd_met = [
        WGRIB2_BIN,
        input_path,
        '-not', 'aerosol=',
        '-set_grib_type','c3',   # -small_grib otherwise writes simple packing (~2x larger)
        '-small_grib',lon_bounds,lat_bounds,
        temp_met
    ]

    cmd_smoke_clip = [
        WGRIB2_BIN,
        input_path,
        '-match',':MASSDEN:.*Particulate organic matter',
        '-set_pdt','0',
        '-set_var','MASSDEN',
        '-set_lev','8 m above ground',
        '-set_ftime',f'{f_hour} hour fcst',
        '-set','background_process_id','134',
        '-set','analysis_or_forecast_process_id','134',
        '-set_grib_type','c3',   # -small_grib otherwise writes simple packing (~2x larger)
        '-small_grib',lon_bounds,lat_bounds,
        temp_smoke_final
    ]

    cmd_smoke_pdt = [
        WGRIB2_BIN,
        input_path,
        '-match',':COLMD:.*Particulate organic matter',
        '-set_pdt','0',
        '-set_var','COLMD',
        '-set_lev','entire atmosphere (considered as a single layer)',
        '-set_ftime',f'{f_hour} hour fcst',
        '-set','background_process_id','134',
        '-set','analysis_or_forecast_process_id','134',
        '-set_grib_type','c3',   # -small_grib otherwise writes simple packing (~2x larger)
        '-small_grib',lon_bounds,lat_bounds,
        temp_smoke_raw
    ]

    try:
        subprocess.check_output(cmd_met, stderr=subprocess.STDOUT)
        subprocess.check_output(cmd_smoke_pdt, stderr=subprocess.STDOUT)
        if os.path.exists(temp_smoke_raw):
            subprocess.check_output(cmd_smoke_clip, stderr=subprocess.STDOUT)

        with open(output_path, 'wb') as outfile:
            for f_name in [temp_met, temp_smoke_raw, temp_smoke_final]:
                if os.path.exists(f_name):
                    with open(f_name, 'rb') as inFile:
                        shutil.copyfileobj(inFile, outfile)
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"Clipping failed: {e.output.decode('utf-8', 'replace')}")
        return False
    except FileNotFoundError:
        log.error(f"wgrib2 binary not found at {WGRIB2_BIN}")
        return False
    finally:
        for temp_file in [temp_met, temp_smoke_raw, temp_smoke_final]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


def process_run(client, date_str, cycle):
    """Deliver every hour not yet delivered. Stops at the first failure; True if all are in."""
    delivered = read_delivered(date_str, cycle)
    for f_hour in F_HOURS:
        if f_hour in delivered:
            continue
        src_filename = SRC_FILENAME_TEMPLATE.format(cycle, f_hour)
        staging_raw = os.path.join(STAGING_DIR, f"raw_{src_filename}")   # full domain
        staging_clip = os.path.join(STAGING_DIR, src_filename)           # clipped result
        dest_filename = DEST_FILENAME_TEMPLATE.format(date_str, cycle, f_hour)
        temp_move_path = os.path.join(FINAL_DIR, f".tmp_{dest_filename}")
        start = time.time()
        try:
            client.download_file(BUCKET_NAME, src_key(date_str, cycle, f_hour), staging_raw)
            if not clip_grib_west(staging_raw, staging_clip, f_hour):
                return False
            shutil.move(staging_clip, temp_move_path)
            os.rename(temp_move_path, os.path.join(FINAL_DIR, dest_filename))  # atomic for ingest
        except Exception as e:
            log.error(f"f{f_hour:03d}: {e}")
            return False
        finally:
            for path in (staging_raw, staging_clip, temp_move_path):
                if os.path.exists(path):
                    os.remove(path)
        with open(delivered_path(date_str, cycle), 'a') as f:
            f.write(f"{f_hour}\n")
        log.info(f"f{f_hour:03d} delivered in {time.time() - start:.0f}s")
    return True


def clean_old_state():
    """Remove SUCCESS_* and delivered_* records older than STATE_MAX_AGE_DAYS."""
    cutoff = time.time() - STATE_MAX_AGE_DAYS * 86400
    for name in os.listdir(STAGING_DIR):
        path = os.path.join(STAGING_DIR, name)
        if name.startswith(('SUCCESS_', 'delivered_')) and os.path.getmtime(path) < cutoff:
            os.remove(path)


def main():
    os.makedirs(STAGING_DIR, exist_ok=True)
    lock = open(os.path.join(STAGING_DIR, '.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        log.info("Previous run still active. Exiting.")
        return 0
    clean_old_state()

    client = get_s3_client()
    try:
        date_str, cycle = find_latest_complete_run(client)
    except Exception as e:
        log.error(f"Could not list the RRFS bucket: {e}")
        return 1
    if not date_str:
        log.error(f"No complete RRFS run in the newest {MAX_LOOKBACK_DAYS} days.")
        return 1
    if os.path.exists(marker_path(date_str, cycle)):
        log.info(f"Run {date_str} {cycle}z already finished. Exiting.")
        return 0

    log.info(f">>> PROCESSING: {date_str} {cycle}z <<<")
    if not process_run(client, date_str, cycle):
        log.error(f"Run {date_str} {cycle}z incomplete; the next run resumes it.")
        return 1
    with open(marker_path(date_str, cycle), 'w') as f:
        f.write("Download and processing complete.")
    if os.path.exists(delivered_path(date_str, cycle)):
        os.remove(delivered_path(date_str, cycle))
    log.info(f"Success! All files for {date_str} {cycle}z moved to {FINAL_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
