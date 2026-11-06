#!/awips2/python/bin/python
#
# Get ECMWF-AIFS and put into LDAD-GRIB
#
# v2025.12.09 - Korri Anderson
# v2026.11.05 - Download over HTTPS from the ECMWF Google Cloud mirror, falling
#               back to the AWS bucket (throttles with SlowDown) and data.ecmwf.int.
#               No longer needs boto3.
# v2026.11.06 - Drop -set_scaling on HGT; default packing is more accurate.
# v2026.11.06 - -new_grid_vectors none: 14 upper-air U/V records were dropped.
# v2026.11.06 - msl relabelled PRMSL so AWIPS stores it as PMSL, not msl-P.
#
import os
import shutil
import datetime
import time
import random
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
STAGING_DIR = '/ldad/localapps/getECMWF/data'
FINAL_DIR = '/data/Incoming'
WGRIB2_BIN = '/ldad/localapps/getHRDPS/wgrib2'

# Look back window (prevents infinite scanning)
MAX_LOOKBACK_DAYS = 3
CYCLES = ['12', '00']
MAX_F_HOUR = 240
STEP_HOURS = 6

# Filename template
# Example: 20251207000000-6h-oper-fc.grib2
SRC_FILENAME_TEMPLATE = "{date}{cycle}0000-{fhour}h-oper-fc.grib2"
DEST_FILENAME_TEMPLATE = "LDAD-GRIB-ECMWF-AIFS_{}_t{}z_f{:03d}.grib2"

def fetch(path, max_retries=6):
    """
    GET a file from the open-data mirrors, in order.
    A mirror that throttles or errors passes to the next. A 404 means the
    file isn't posted yet. If every mirror fails, back off and try again.
    """
    name = os.path.basename(path)
    for attempt in range(max_retries):
        codes = []
        for base in MIRRORS:
            try:
                req = urllib.request.Request(f'{base}/{path}',
                                             headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=120) as resp:
                    return resp.read()
            except urllib.error.HTTPError as e:
                codes.append(e.code)
                if e.code != 404:
                    print(f"\n [WARN] HTTP {e.code} from {base} for {name}")
            except (urllib.error.URLError, OSError) as e:
                codes.append(None)
                print(f"\n [WARN] {base} failed for {name}: {e}")
        if 404 in codes:  # a mirror says not posted; next cron run retries
            raise IOError(f"{name} not posted yet")
        wait = min(5 * 2 ** attempt, 120) + random.uniform(0, 1)
        print(f" [SLOWDOWN] All mirrors failed for {name}. Waiting {wait:.0f}s "
              f"(attempt {attempt + 1}/{max_retries})")
        time.sleep(wait)
    raise IOError(f"Could not download {name} after {max_retries} attempts")

def file_exists(path):
    try:
        print(f"  [CHECKING] {path}...", end='', flush=True)
        fetch(path.replace('.grib2', '.index'))  # small file, same posting
        print(" FOUND!")
        return True
    except Exception:
        print(" Not found.")
        return False

def find_latest_complete_run():
    """
    Calculates dates for the last few days and checks specific paths directly.
    """
    print("Searching for latest complete run...")

    # Use timezone-aware UTC
    start_date = datetime.datetime.now(datetime.timezone.utc)

    for i in range(MAX_LOOKBACK_DAYS):
        # Calculate date string (e.g., 20251207)
        current_date_obj = start_date - datetime.timedelta(days=i)
        date_str = current_date_obj.strftime("%Y%m%d")

        for cycle in CYCLES:
            # Check for the LAST file (MAX_F_HOUR) to confirm completion
            filename = SRC_FILENAME_TEMPLATE.format(
                date=date_str,
                cycle=cycle,
                fhour=MAX_F_HOUR
            )

            # Path: 20251207/00z/aifs-single/0p25/oper/FILENAME
            key = f"{date_str}/{cycle}z/aifs-single/0p25/oper/{filename}"

            if file_exists(key):
                return date_str, cycle

    return None, None

def is_run_processed(date_str, cycle):
    marker = os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")
    return os.path.exists(marker)

def mark_run_complete(date_str, cycle):
    if not os.path.exists(STAGING_DIR): os.makedirs(STAGING_DIR)
    marker = os.path.join(STAGING_DIR, f"SUCCESS_{date_str}_{cycle}z")
    with open(marker, 'w') as f:
        f.write(f"Download complete at {datetime.datetime.now(datetime.timezone.utc)}")
    print(f"[INFO] Run marked complete: {marker}")

def download(path, dest_path):
    data = fetch(path)
    with open(dest_path, 'wb') as f:
        f.write(data)

def sanitize_grib(input_path, output_path):
    """
    Two-Step Processing to fix ECMWF Grid issues:
    Step 1: Convert N->S Gaussian grid to S->N order (Fixes 'fatal error').
    Step 2: Apply Math/Renaming/Regridding on the clean file.
    """

    # Intermediate temp file for South-to-North data
    temp_sn_path = input_path + ".sn"

    # --- STEP 1: FORCE ORDER SOUTH-TO-NORTH ---
    # We strip all processing and just focus on fixing the grid order.
    cmd_step1 = [
        WGRIB2_BIN,
        input_path,
        # --- REGRIDDING ---
        '-set_grib_type', 'c3',      # Complex Packing (Efficient)
        '-new_grid_winds', 'earth',  # Rotate winds
        # Winds are already earth-relative, and ECMWF does not keep U and V
        # adjacent, so wgrib2 would drop the unpaired ones.
        '-new_grid_vectors', 'none',
        '-new_grid', 'latlon',
        '0:1440:0.25',               # Lon: 0 to 359.75
        '0:361:0.25',                # Lat: 0 to 90 (North)
        temp_sn_path
    ]

    # --- STEP 2: PROCESS & REGRID ---
    # Now we work on the clean temp_sn_path
    cmd_step2 = [
        WGRIB2_BIN,
        temp_sn_path,         # INPUT is now the temp file

        # --- FIX: Geopotential (z) -> Geopotential Height (GH) ---
        '-if', ':(GP):',
            '-rpn','9.80665:/',
            '-set_var', 'HGT',
            # No -set_scaling: 1 0 packed to 8 bits (errors up to 5 m);
            # wgrib2's default keeps heights within 0.5 m.
        '-endif',

        # --- FIX: Total Precipitation (tp) ---
        '-if', ':(TPRATE|TP|var0_1_228):',
             #'-rpn','1000:*',
             '-set_var', 'APCP',
             #'-set_scaling','1','0',
        '-endif',

        # --- FIX: Total Snowfall (sf) ---
        #'-if', ':(TSRWE|SF|var0_1_144):',
        #     '-rpn','1000:*',
        #     '-set_var', 'WEASD'
        #'-endif',

        # --- FIX: CONVECTIVE Precipitation ---
        #'-if', ':(CPRATE|CP|var0_1_143):',
        #     '-rpn','1000:*',
        #     '-set_var', 'ACPCP'
        #'-endif',

        # --- FIX: Volumetric Soil Water (vsw) ---
        #'-if', ':var2_3_192:',
        #    '-set_var', 'SOILW',
        #'-endif',

        # --- FIX: Soil Temperature (sot) ---
        #'-if', ':var2_3_18:',
        #    '-set_var', 'TSTMC',
        #'-endif',

        # --- FIX: Mean Sea Level Pressure (msl) ---
        # ECMWF codes msl as plain pressure (0-3-0) at MSL, which AWIPS stores
        # as msl-P; PRMSL (0-3-1) is what AWIPS stores as PMSL, as for GFS.
        '-if', ':PRES:mean sea level:',
            '-set_var', 'PRMSL',
        '-endif',

        '-set_grib_type','c3',
        '-grib_out', output_path
    ]

    print(f"   [WGRIB2] Step 1: Reordering {os.path.basename(input_path)}...")
    try:
        # Run Step 1
        subprocess.check_output(cmd_step1, stderr=subprocess.STDOUT)

        print(f"   [WGRIB2] Step 2: Remapping & Regridding...")
        # Run Step 2
        subprocess.check_output(cmd_step2, stderr=subprocess.STDOUT)

        # Cleanup Temp File
        if os.path.exists(temp_sn_path):
            os.remove(temp_sn_path)

        return True

    except subprocess.CalledProcessError as e:
        print(f"   [ERROR] wgrib2 failed: {e.output.decode('utf-8')}")
        # Cleanup on failure
        if os.path.exists(temp_sn_path):
            os.remove(temp_sn_path)
        return False
    except FileNotFoundError:
        print(f"   [ERROR] wgrib2 binary not found at: {WGRIB2_BIN}")
        return False

def process_run(date_str, cycle):
    print(f"\n>>> PROCESSING: {date_str} {cycle}z <<<")

    if not os.path.exists(STAGING_DIR): os.makedirs(STAGING_DIR)
    if not os.path.exists(FINAL_DIR): os.makedirs(FINAL_DIR)

    for f_hour in range(0, MAX_F_HOUR + 1, STEP_HOURS):

        # 1. Setup Filenames
        src_filename = SRC_FILENAME_TEMPLATE.format(
            date=date_str,
            cycle=cycle,
            fhour=f_hour
        )

        src_path = f"{date_str}/{cycle}z/aifs-single/0p25/oper/{src_filename}"

        # Staging paths
        staging_raw_path = os.path.join(STAGING_DIR, f"raw_{src_filename}")
        staging_clean_path = os.path.join(STAGING_DIR, src_filename)

        # Final paths
        dest_filename = DEST_FILENAME_TEMPLATE.format(date_str, cycle, f_hour)
        final_path = os.path.join(FINAL_DIR, dest_filename)
        temp_move_path = os.path.join(FINAL_DIR, f".tmp_{dest_filename}")

        # 2. Check Exists
        if os.path.exists(final_path):
            print(f"   [EXISTS] {dest_filename}")
            continue

        # 3. Download & Process
        print(f"   Downloading {src_filename}...", end="\r")
        try:
            # Download Raw
            download(src_path, staging_raw_path)
            time.sleep(1.5)
            # Sanitize (Reproject + Compress + Clip)
            if sanitize_grib(staging_raw_path, staging_clean_path):
                # Move Clean file to Final
                shutil.move(staging_clean_path, temp_move_path)
                os.rename(temp_move_path, final_path)

                # Cleanup Raw
                if os.path.exists(staging_raw_path): os.remove(staging_raw_path)
            else:
                print(f"   [SKIP] Conversion failed for {src_filename}")
                if os.path.exists(staging_raw_path): os.remove(staging_raw_path)
                return False

        except Exception as e:
            print(f"\n   [ERROR] Failed on {src_filename}: {e}")
            return False

    print(f"\n   Success! All files for {cycle}z moved to {FINAL_DIR}")
    return True

def main():
    # 1. Find Run using calculated dates
    date_str, cycle = find_latest_complete_run()

    if not date_str:
        print("[ERROR] No complete runs found in the last 3 days.")
        return

    print(f"[INFO] Target found: {date_str} {cycle}z")

    # 2. Check Marker
    if is_run_processed(date_str, cycle):
        print(f"[INFO] Run {date_str} {cycle}z already finished. Exiting.")
        return

    # 3. Process
    if process_run(date_str, cycle):
        mark_run_complete(date_str, cycle)

if __name__ == "__main__":
    main()
