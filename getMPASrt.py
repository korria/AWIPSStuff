#!/awips2/python/bin/python

"""
getMPASrt2.py - Download MPAS-RT from NSSL, clip to PNW, ingest to AWIPS

HISTORY:
  - 2024-06-24 - Initial Version - Korri Anderson
  - 2024-12-23 - Change MPASRT to MPASHT
  - 2025-01-03 - Fix handling and clear directory when downloads fail.
  - 2025-06-30 - Update to use new NSSL SFTP server
  - 2025-09-07 - Update to MPASrn3
  - 2025-12-09 - Added wgrib2 clipping for Pacific Northwest
  - 2026-01-14 - Refactored: state tracking, resume capability, proper logging
  - 2026-10-26 - Switch from SFTP to the NSSL THREDDS archive (one .tar per
                 cycle over HTTPS); stream the tar and gunzip grib2 members
                 on the fly. Drops the paramiko dependency.
"""
import os
import re
import gzip
import sys
import json
import time
import shutil
import socket
import logging
import tarfile
import subprocess
import urllib.error
import urllib.request
from datetime import datetime, timedelta

# =============================================================================
# Configuration
# =============================================================================

# THREDDS archive: one tar per cycle, e.g. .../2026/26102600_mpasrn.tar
# Catalog: https://data.nssl.noaa.gov/thredds/catalog/FRDD/NSSL-MPAS/2026/catalog.html
TAR_URL = ('https://data.nssl.noaa.gov/thredds/fileServer/FRDD/NSSL-MPAS/'
           '{year}/{yymmddhh}_{run}.tar')
RUN_NAME = 'mpasrn'          # also on the server: mpasht2 (00/06/12/18Z)
USER_AGENT = 'Mozilla/5.0'

# Local Paths
BASE_DIR = '/ldad/localapps/getMPASRT'
LOCAL_DIRECTORY = os.path.join(BASE_DIR, 'downloaded_files')
TARGET_DIRECTORY = '/data/Incoming/'
STATE_FILE = os.path.join(BASE_DIR, 'download_state.json')
LOG_FILE = os.path.join(BASE_DIR, 'getMPASrt2.log')

# Clipping Configuration (Pacific Northwest)
WGRIB2_CMD = '/ldad/localapps/getHRDPS/wgrib2'
MIN_LON, MAX_LON = 1, 750
MIN_LAT, MAX_LAT = 1, 1059

# Processing Configuration
FORECAST_HOURS = 84          # f01 through f84
DAYS_TO_CHECK = 2            # Check today and yesterday
CYCLES = ['00', '12']
MAX_RETRIES = 3
RETRY_DELAY = 60             # seconds
STATE_RETENTION_DAYS = 7     # Prune state older than this

# Forecast hour from a tar member name, e.g.
# "./2026102600/mpasrn_nssl_2026102600f07.grib2.gz" (members are gzipped;
# the empty "..._ready" marker files are skipped)
FHR_RE = re.compile(r'f(\d{2,3})\.grib2\.gz$')

# =============================================================================
# Logging Setup
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)

# =============================================================================
# State Management
# =============================================================================

def load_state():
    """Load download state from disk."""
    if not os.path.exists(STATE_FILE):
        return {}
    try:
        with open(STATE_FILE, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        log.warning(f"Could not load state file: {e}")
        return {}


def save_state(state):
    """Persist download state to disk."""
    try:
        with open(STATE_FILE, 'w') as f:
            json.dump(state, f, indent=2)
    except IOError as e:
        log.error(f"Could not save state file: {e}")


def prune_old_state(state):
    """Remove cycles older than STATE_RETENTION_DAYS."""
    cutoff = datetime.utcnow() - timedelta(days=STATE_RETENTION_DAYS)
    cutoff_str = cutoff.strftime('%Y%m%d%H')

    original_count = len(state)
    pruned = {k: v for k, v in state.items() if k >= cutoff_str}

    if len(pruned) < original_count:
        log.info(f"Pruned {original_count - len(pruned)} old cycle(s) from state")

    return pruned

# =============================================================================
# Directory Management
# =============================================================================

def ensure_directories():
    """Create required directories if they don't exist."""
    os.makedirs(LOCAL_DIRECTORY, exist_ok=True)
    os.makedirs(TARGET_DIRECTORY, exist_ok=True)


def cleanup_local_directory():
    """Remove all files from the local download directory."""
    for fname in os.listdir(LOCAL_DIRECTORY):
        fpath = os.path.join(LOCAL_DIRECTORY, fname)
        try:
            os.remove(fpath)
            log.debug(f"Cleaned up {fname}")
        except OSError as e:
            log.warning(f"Failed to remove {fname}: {e}")

# =============================================================================
# File Processing
# =============================================================================

def clip_to_pnw(input_path, output_path):
    """
    Clip GRIB2 file to Pacific Northwest domain using wgrib2.

    Raises subprocess.CalledProcessError on failure.
    """
    cmd = [
        WGRIB2_CMD,
        input_path,
        '-set_grib_type','same',
        '-ijsmall_grib',
        f'{MIN_LON}:{MAX_LON}',
        f'{MIN_LAT}:{MAX_LAT}',
        output_path
    ]

    result = subprocess.run(cmd, check=True, capture_output=True)
    return result


def process_member(tar, member, local_file):
    """
    Extract and gunzip one grib2 member from the tar stream, clip it, and move it to
    the target directory.

    Returns True on success, False if wgrib2 failed. Network/tar errors
    propagate so the caller can retry the cycle.
    """
    temp_path = os.path.join(LOCAL_DIRECTORY, f"full_{local_file}")
    clipped_path = os.path.join(LOCAL_DIRECTORY, local_file)
    final_path = os.path.join(TARGET_DIRECTORY, local_file)

    try:
        with gzip.GzipFile(fileobj=tar.extractfile(member)) as src, \
                open(temp_path, 'wb') as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)

        clip_to_pnw(temp_path, clipped_path)
        shutil.move(clipped_path, final_path)
        log.info(f"Processed {member.name} -> {final_path}")
        return True

    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode() if e.stderr else 'Unknown error'
        log.error(f"wgrib2 clip failed for {member.name}: {stderr}")
        return False

    finally:
        for path in (temp_path, clipped_path):
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass


def open_cycle_tar(date, hour):
    """
    Open the cycle's tar over HTTPS. Returns the response, or None if the
    cycle isn't on the server yet.
    """
    url = TAR_URL.format(year=date.strftime('%Y'),
                         yymmddhh=date.strftime('%y%m%d') + hour,
                         run=RUN_NAME)
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    try:
        return urllib.request.urlopen(req)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise


def stream_cycle(resp, cycle_id, local_prefix, state):
    """
    Walk the tar stream once, processing every wanted forecast hour not
    already in state. Returns the number of newly processed hours.
    """
    # ponytail: a retry re-streams the whole tar (~8 GB) and skips done
    # members; resume via HTTP Range at member offsets if bandwidth matters.
    done = 0
    with tarfile.open(fileobj=resp, mode='r|') as tar:
        for member in tar:
            match = FHR_RE.search(member.name) if member.isfile() else None
            if not match:
                log.debug(f"Skipping tar member {member.name}")
                continue

            fhr = int(match.group(1))
            if not 1 <= fhr <= FORECAST_HOURS or fhr in state[cycle_id]:
                continue

            local_file = f'LDAD-GRIB-{local_prefix}f{fhr:03d}.grib2'
            if process_member(tar, member, local_file):
                state[cycle_id].append(fhr)
                save_state(state)
                done += 1
    return done


def process_cycle(date, hour, state):
    """
    Process all forecast hours for a single model cycle.

    Returns True if cycle is complete, False otherwise.
    """
    date_str = date.strftime('%Y%m%d')
    cycle_id = f"{date_str}{hour}"
    local_prefix = f'mpasrt4nssl_{date_str}{hour}'

    state.setdefault(cycle_id, [])

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = open_cycle_tar(date, hour)
            if resp is None:
                log.debug(f"Cycle {cycle_id} not available yet")
                return False

            log.info(f"Processing cycle {cycle_id} from {resp.url}")
            with resp:
                stream_cycle(resp, cycle_id, local_prefix, state)
            break

        except (urllib.error.URLError, tarfile.TarError,
                socket.error, IOError) as e:
            log.warning(f"Attempt {attempt}/{MAX_RETRIES} for {cycle_id} failed: {e}")
            if attempt == MAX_RETRIES:
                raise
            time.sleep(RETRY_DELAY)

    count = len(state[cycle_id])
    complete = count == FORECAST_HOURS
    log.info(f"Cycle {cycle_id}: {'complete' if complete else f'{count}/{FORECAST_HOURS}'}")
    return complete

# =============================================================================
# Main
# =============================================================================

def main():
    """Main entry point."""
    log.info("=" * 60)
    log.info("getMPASrt2.py starting")

    if not os.path.exists(WGRIB2_CMD):
        log.error(f"wgrib2 not found at {WGRIB2_CMD}")
        sys.exit(1)

    ensure_directories()
    socket.setdefaulttimeout(300)

    # Load and prune state
    state = load_state()
    state = prune_old_state(state)
    save_state(state)

    cycles_processed = 0
    cycles_completed = 0

    try:
        current_date = datetime.utcnow()

        for days_back in range(DAYS_TO_CHECK):
            date = current_date - timedelta(days=days_back)

            for hour in CYCLES:
                cycle_id = f"{date.strftime('%Y%m%d')}{hour}"

                # Skip if already complete
                if len(state.get(cycle_id, [])) == FORECAST_HOURS:
                    log.debug(f"Skipping {cycle_id} (already complete)")
                    continue

                try:
                    if process_cycle(date, hour, state):
                        cycles_completed += 1
                    cycles_processed += 1

                except Exception as e:
                    log.error(f"Failed to process cycle {cycle_id}: {e}")
                    # Continue to next cycle instead of failing entirely
                    continue

    except Exception as e:
        log.exception(f"Fatal error: {e}")
        cleanup_local_directory()
        sys.exit(1)

    # Summary
    if cycles_processed == 0:
        log.info("No new cycles to process")
    else:
        log.info(f"Processed {cycles_processed} cycle(s), {cycles_completed} complete")

    log.info("getMPASrt2.py finished")


if __name__ == '__main__':
    main()
