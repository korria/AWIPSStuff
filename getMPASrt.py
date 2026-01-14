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
"""
import os
import sys
import json
import time
import shutil
import socket
import logging
import subprocess
from datetime import datetime, timedelta
from contextlib import contextmanager

sys.path.insert(0, '/ldad/localapps/getMPASRT/lib')
import paramiko

# =============================================================================
# Configuration
# =============================================================================

# SFTP Connection
HOST = 'cftp.nssl.noaa.gov'
PORT = 2222
USERNAME = ''
PASSWORD = ''
REMOTE_DIR = '/realtime/'

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
FORECAST_HOURS = 85          # f00 through f84
DAYS_TO_CHECK = 2            # Check today and yesterday
CYCLES = ['00', '12']
MAX_RETRIES = 3
RETRY_DELAY = 5              # seconds
STATE_RETENTION_DAYS = 7     # Prune state older than this

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
# SFTP Connection
# =============================================================================

@contextmanager
def sftp_connection():
    """Context manager for SFTP connections with automatic cleanup."""
    transport = None
    sftp = None
    try:
        log.info(f"Connecting to {HOST}:{PORT}")
        transport = paramiko.Transport((HOST, PORT))
        transport.connect(username=USERNAME, password=PASSWORD)
        sftp = paramiko.SFTPClient.from_transport(transport)
        yield sftp
    finally:
        if sftp:
            try:
                sftp.close()
            except Exception:
                pass
        if transport:
            try:
                transport.close()
            except Exception:
                pass
        log.debug("SFTP connection closed")


def reconnect_sftp():
    """Create a new SFTP connection (for use in retry logic)."""
    log.info(f"Reconnecting to {HOST}:{PORT}")
    transport = paramiko.Transport((HOST, PORT))
    transport.connect(username=USERNAME, password=PASSWORD)
    sftp = paramiko.SFTPClient.from_transport(transport)
    return sftp, transport

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


def download_and_process_file(sftp, remote_file, local_file):
    """
    Download a single file, clip it, and move to target directory.

    Returns True on success, False on failure.
    """
    temp_path = os.path.join(LOCAL_DIRECTORY, f"full_{local_file}")
    clipped_path = os.path.join(LOCAL_DIRECTORY, local_file)
    final_path = os.path.join(TARGET_DIRECTORY, local_file)

    try:
        # Download from SFTP
        sftp.get(os.path.join(REMOTE_DIR, remote_file), temp_path)
        log.debug(f"Downloaded {remote_file}")

        # Clip to PNW region
        clip_to_pnw(temp_path, clipped_path)
        log.debug(f"Clipped {local_file}")

        # Move to AWIPS incoming directory
        shutil.move(clipped_path, final_path)
        log.info(f"Processed {remote_file} -> {final_path}")

        return True

    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode() if e.stderr else 'Unknown error'
        log.error(f"wgrib2 clip failed for {remote_file}: {stderr}")
        return False

    except IOError as e:
        log.error(f"I/O error processing {remote_file}: {e}")
        raise  # Re-raise to trigger reconnection logic

    except Exception as e:
        log.error(f"Unexpected error processing {remote_file}: {e}")
        raise

    finally:
        # Always clean up temp files
        for path in (temp_path, clipped_path):
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass

# =============================================================================
# Cycle Processing
# =============================================================================

def cycle_available(sftp, prefix):
    """Check if a complete cycle is available on the server."""
    test_file = f'{prefix}f{FORECAST_HOURS - 1:02d}.grib2'
    try:
        sftp.stat(os.path.join(REMOTE_DIR, test_file))
        return True
    except FileNotFoundError:
        return False


def process_cycle(sftp_holder, date_str, hour, state):
    """
    Process all forecast hours for a single model cycle.

    sftp_holder is a mutable dict {'sftp': sftp, 'transport': transport}
    so we can reconnect if needed.

    Returns True if cycle is complete, False otherwise.
    """
    cycle_id = f"{date_str}{hour}"
    prefix = f'mpasrn3_nssl_{date_str}{hour}'
    local_prefix = f'mpasrt4nssl_{date_str}{hour}'

    sftp = sftp_holder['sftp']

    # Check if cycle exists on server
    if not cycle_available(sftp, prefix):
        log.debug(f"Cycle {cycle_id} not available yet")
        return False

    log.info(f"Processing cycle {cycle_id}")

    # Initialize state for this cycle if needed
    if cycle_id not in state:
        state[cycle_id] = []

    success_count = len(state[cycle_id])  # Count already-processed files

    for fhr in range(FORECAST_HOURS):
        # Skip if already downloaded
        if fhr in state[cycle_id]:
            log.debug(f"Skipping {cycle_id} f{fhr:02d} (already processed)")
            continue

        remote_file = f'{prefix}f{fhr:02d}.grib2'
        local_file = f'LDAD-GRIB-{local_prefix}f{fhr:03d}.grib2'

        for attempt in range(MAX_RETRIES):
            try:
                if download_and_process_file(sftp_holder['sftp'], remote_file, local_file):
                    success_count += 1
                    state[cycle_id].append(fhr)
                    save_state(state)
                    break
                else:
                    # Clip failed but connection is okay, retry
                    log.warning(f"Attempt {attempt + 1}/{MAX_RETRIES} failed for {remote_file}")
                    time.sleep(RETRY_DELAY)

            except (IOError, paramiko.SSHException, socket.error) as e:
                # Connection issue - try to reconnect
                log.warning(f"Connection error on {remote_file}: {e}")
                time.sleep(RETRY_DELAY)

                try:
                    # Close old connection
                    try:
                        sftp_holder['sftp'].close()
                        sftp_holder['transport'].close()
                    except Exception:
                        pass

                    # Reconnect
                    new_sftp, new_transport = reconnect_sftp()
                    sftp_holder['sftp'] = new_sftp
                    sftp_holder['transport'] = new_transport
                    log.info("Reconnected successfully")

                except Exception as reconnect_error:
                    log.error(f"Reconnection failed: {reconnect_error}")
                    if attempt == MAX_RETRIES - 1:
                        raise

        else:
            log.error(f"Skipping {remote_file} after {MAX_RETRIES} failed attempts")

    complete = success_count == FORECAST_HOURS
    status = "complete" if complete else f"{success_count}/{FORECAST_HOURS}"
    log.info(f"Cycle {cycle_id}: {status}")

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
        with sftp_connection() as sftp:
            # We need a mutable holder so process_cycle can reconnect
            transport = sftp.get_channel().get_transport()
            sftp_holder = {'sftp': sftp, 'transport': transport}

            current_date = datetime.utcnow()

            for days_back in range(DAYS_TO_CHECK):
                date = current_date - timedelta(days=days_back)
                date_str = date.strftime('%Y%m%d')

                for hour in CYCLES:
                    cycle_id = f"{date_str}{hour}"

                    # Skip if already complete
                    if cycle_id in state and len(state[cycle_id]) == FORECAST_HOURS:
                        log.debug(f"Skipping {cycle_id} (already complete)")
                        continue

                    try:
                        if process_cycle(sftp_holder, date_str, hour, state):
                            cycles_completed += 1
                        cycles_processed += 1

                    except Exception as e:
                        log.error(f"Failed to process cycle {cycle_id}: {e}")
                        # Continue to next cycle instead of failing entirely
                        continue

    except paramiko.AuthenticationException:
        log.error("Authentication failed. Check NSSL_USER and NSSL_PASS.")
        sys.exit(1)

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
