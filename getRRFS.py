#!/awips2/python/bin/python
import itertools
import os
import shutil
import sys
import subprocess
sys.path.insert(0, "/ldad/localapps/getRRFS/bin/")
import boto3
from botocore import UNSIGNED
from botocore.client import Config

# --- CONFIGURATION ---
BUCKET_NAME = 'noaa-rrfs-ops-pds'
PREFIX_ROOT = ''

# Directories
STAGING_DIR = '/ldad/localapps/getRRFS/data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

# Search Priority
CYCLES = ['18','12','06','00']
MAX_F_HOUR = 72

#Different timestamps 1hr vs. 3hr
r1 = range(1,61,1) # 1-60 hourly
r2 = range(61,85,1) # 63-84 3-hourly

# Western US Clipping Bounds (Longitude: West:East, Latitude: South:North)
# This covers WA, OR, CA, NV, ID, UT, AZ and parts of MT, WY, CO, NM
WEST_BOUNDS = "-126:-101 31:50"

# Source Pattern (S3)
SRC_FILENAME_TEMPLATE = "rrfs.t{}z.2dfld.3km.f{:03d}.conus.grib2"

# Destination Pattern (LDAD) - Added _West tag
DEST_FILENAME_TEMPLATE = "LDAD-GRIB-LocalRRFS_West_{}_t{}z_f{:03d}.grib2"

def get_s3_client():
    return boto3.client('s3', config=Config(signature_version=UNSIGNED))

def check_s3_file_exists(client, bucket, key):
    try:
        client.head_object(Bucket=bucket, Key=key)
        return True
    except:
        return False

def find_latest_complete_run(client):
    """
    Finds the latest run (YYYYMMDD + Cycle) that has hour f084 available.
    """
    print("Scanning S3 for the latest complete run...")

    result = client.list_objects_v2(Bucket=BUCKET_NAME, Prefix=PREFIX_ROOT, Delimiter='/')
    if 'CommonPrefixes' not in result:
        return None, None, None

    # Sort dates newest to oldest
    all_date_folders = [x['Prefix'] for x in result['CommonPrefixes']]
    all_date_folders.sort(reverse=True)

    for date_folder in all_date_folders:
        # Extract "20251128" from "rrfs_a/rrfs.20251128/"
        date_str = date_folder.strip('/').split('.')[-1]

        for cycle in CYCLES:
            # Check for f084
            last_file = SRC_FILENAME_TEMPLATE.format(cycle, MAX_F_HOUR)
            key = f"{date_folder}{cycle}/{last_file}"

            if check_s3_file_exists(client, BUCKET_NAME, key):
                return date_folder, date_str, cycle

    return None, None, None

def is_run_processed(date_str, cycle):
    marker = os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")
    return os.path.exists(marker)

def mark_run_complete(date_str, cycle):
    marker = os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")
    with open(marker, 'w') as f:
        f.write("Download and processing complete.")
    print(f"[INFO] Run marked complete: {marker}")

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
        '-small_grib',lon_bounds,lat_bounds,
        temp_smoke_raw   
    ]


    # print(f"   [WGRIB2] Clipping {os.path.basename(input_path)}...")

    try:
        # Run command, capture output to suppress noise
        subprocess.check_output(cmd_met, stderr=subprocess.STDOUT)
        subprocess.check_output(cmd_smoke_pdt, stderr=subprocess.STDOUT)
        if os.path.exists(temp_smoke_raw):
           subprocess.check_output(cmd_smoke_clip, stderr=subprocess.STDOUT)

        with open(output_path,'wb') as outfile:
           for f_name in [temp_met,temp_smoke_raw,temp_smoke_final]:
              if os.path.exists(f_name):
                 with open(f_name, 'rb') as inFile:
                    shutil.copyfileobj(inFile, outfile)
        for temp_file in [temp_met,temp_smoke_raw,temp_smoke_final]:
           if os.path.exists(temp_file): 
              os.remove(temp_file)
        return True

    except subprocess.CalledProcessError as e:
        print(f"   [ERROR] Clipping failed: {e.output.decode('utf-8')}")
        for temp_file in [temp_met,temp_smoke_raw,temp_smoke_final]:
           if os.path.exists(temp_file): 
              os.remove(temp_file)
        return False

    except FileNotFoundError:
        print(f"   [ERROR] wgrib2 binary not found at {WGRIB2_BIN}")
        return False

def process_run(client, date_prefix, date_str, cycle):
    print(f"\n>>> PROCESSING: {date_str} {cycle}z <<<")

    if not os.path.exists(STAGING_DIR): os.makedirs(STAGING_DIR)
    if not os.path.exists(FINAL_DIR): os.makedirs(FINAL_DIR)

    for f_hour in itertools.chain(r1,r2):
        # 1. Define Names
        src_filename = SRC_FILENAME_TEMPLATE.format(cycle, f_hour)
        s3_key = f"{date_prefix}{cycle}/{src_filename}"

        # Staging Paths
        staging_raw = os.path.join(STAGING_DIR, f"raw_{src_filename}")   # Full Domain
        staging_clip = os.path.join(STAGING_DIR, src_filename)           # Clipped Result

        # Final Path
        dest_filename = DEST_FILENAME_TEMPLATE.format(date_str, cycle, f_hour)
        final_path = os.path.join(FINAL_DIR, dest_filename)
        temp_move_path = os.path.join(FINAL_DIR, f".tmp_{dest_filename}")

        # 2. Check Exists
        if os.path.exists(final_path):
            print(f"   [EXISTS] {dest_filename}")
            continue

        # 3. Download Full Domain
        print(f"   Downloading f{f_hour:03d}...", end="\r")
        try:
            client.download_file(BUCKET_NAME, s3_key, staging_raw)
        except Exception as e:
            print(f"\n   [ERROR] Download failed for {src_filename}: {e}")
            return False

        # 4. Clip to Western US
        if clip_grib_west(staging_raw, staging_clip, f_hour):
            # Success: Move Clipped file to Final
            try:
                shutil.move(staging_clip, temp_move_path)
                os.rename(temp_move_path, final_path)
            except Exception as e:
                print(f"\n   [ERROR] Move failed: {e}")
                return False

            # Cleanup Raw File
            if os.path.exists(staging_raw):
                os.remove(staging_raw)
        else:
            # Failure: Cleanup and abort
            print(f"\n   [ERROR] Clipping failed for {src_filename}")
            if os.path.exists(staging_raw):
                os.remove(staging_raw)
            return False

    print(f"\n   Success! All files for {cycle}z moved to {FINAL_DIR}")
    return True

def main():
    client = get_s3_client()

    # 1. Find Run
    date_prefix, date_str, cycle = find_latest_complete_run(client)
    if not date_prefix:
        print("[ERROR] No complete runs found on S3.")
        return

    print(f"[INFO] Target: {date_str} {cycle}z")

    # 2. Check Marker
    if is_run_processed(date_str, cycle):
        print(f"[INFO] Run {date_str} {cycle}z already finished. Exiting.")
        return

    # 3. Process
    if process_run(client, date_prefix, date_str, cycle):
        mark_run_complete(date_str, cycle)

if __name__ == "__main__":
    main()
