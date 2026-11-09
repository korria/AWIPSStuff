#!/awips2/python/bin/python
#
# Get ECMWF-IFS 0.1 deg open data and put into LDAD-GRIB
#
# Modeled on getECAI.py (ECMWF-AIFS). Pulls only the wanted fields out of the
# global 0.1 deg files using the .index byte offsets (~470 MB/hour otherwise),
# then regrids to a North America + NE Pacific 0.1 deg lat-lon box, south-to-north.
#
# v2026.11.05 - Korri Anderson - Initial version
# v2026.11.05 - Pull from the Google Cloud mirror first (AWS S3 bucket was returning
#               SlowDown), then S3, then data.ecmwf.int; back off when throttled;
#               fewer, larger range requests. No longer needs boto3.
# v2026.11.05 - 6-hourly 0-240 h (was 3-hourly to 144 h).
# v2026.11.05 - tp -> APCP (x1000, m -> mm), same as getECAI. AWIPS did not
#               decode the raw ECMWF local code 0-1-193.
# v2026.11.06 - Gust labelled instantaneous so AWIPS shows WGS at every hour;
#               f000 gust (all zeros) dropped.
# v2026.11.06 - msl relabelled PRMSL so AWIPS stores it as PMSL, not msl-P.
# v2026.11.06 - Cloud cover tcc/lcc/mcc/hcc -> TCDC/LCDC/MCDC/HCDC (%), as AIFS.
# v2026.11.06 - cp -> CPRAT and sf -> TSRWE in mm, as AIFS.
# v2026.11.06 - Domain 180-40W, 15-70N (1401 x 551), was 160-60W, 20-60N.
# v2026.11.06 - Add surface z (f000 only) as HGT (m) terrain, same as getECAI.
# v2026.11.09 - GFS-style fields: BL wind (copy of 100 m wind labelled 30-0 mb
#               above ground), PW from tcwv, mucape relabelled 255-0 mb above
#               ground (AWIPS has no level 17), WEASD (mm) and SNOD (m) from
#               sd and snow density rsn.
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
    'tcwv', 'sd',
}
# Snow density, downloaded ahead of SFC_PARAMS so sanitize_grib sees it before
# sd (it is needed to turn sd into snow depth).
FIRST_SFC_PARAMS = {'rsn'}
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
# 180-40W, 15N-70N at 0.1 deg = 1401 x 551.
OUT_LON = '180:1401:0.1'
OUT_LAT = '15:551:0.1'

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
        for kind, params, is_pl in (('sfc', FIRST_SFC_PARAMS, False),
                                    ('sfc', SFC_PARAMS, False), ('pl', PL_PARAMS, True)):
            base = src_key(date_str, cycle, f_hour, kind)
            index = fetch(base + '.index').decode()
            for span_start, span_end, pieces in group_ranges(wanted_ranges(index, params, is_pl)):
                data = fetch(base + '.grib2', (span_start, span_end))
                for start, end in pieces:
                    out.write(data[start - span_start:end - span_start + 1])


def sanitize_grib(input_path, output_path, f_hour):
    """
    Step 1: regrid to the output box (also flips ECMWF north-to-south order
            to south-to-north). Nearest neighbour: the box sits on the source
            0.1 deg points, so values are copied, and ptype stays categorical.
    Step 2: tp (ECMWF local 0-1-193, metres) -> APCP (kg m-2 = mm), as in
            getECAI; surface geopotential z (f000 only) -> HGT terrain (m);
            GFS-style PW, MUCAPE, BL wind, WEASD and SNOD;
            append CRAIN/CSNOW/CFRZR/CICEP derived from PTYPE.
    wgrib2 -if matches the record's original name, even after -set_var.
    """
    temp_path = input_path + ".box"
    # sd is ECMWF local 0-1-254 (m of water); names vary by wgrib2 build, as tp.
    sd = ':(sd|var0_1_254|var0_1_1_98_1_254):'
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
    ] + (['-not', ':GUST:'] if f_hour == 0 else []) + [  # f000 gust is all zeros
        # wgrib2 names 0-1-193 tprate (newer) or var0_1_193 (older builds).
        '-if', ':(tprate|var0_1_193|var0_1_1_98_1_193):',
            '-rpn', '1000:*',
            '-set_var', 'APCP',
        '-endif',
        '-if', ':GP:',
            '-rpn', '9.80665:/',
            '-set_var', 'HGT',
        '-endif',
        # cp and sf are ECMWF local codes in metres; write mm with the WMO
        # codes AIFS uses (CPRAT 0-1-37, TSRWE 0-1-53, run totals).
        '-if', ':(cp|var0_1_195|var0_1_1_98_1_195):',
            '-rpn', '1000:*', '-set_var', 'CPRAT', '-endif',
        '-if', ':(sf|mtsfr|var0_1_198|var0_1_1_98_1_198):',
            '-rpn', '1000:*', '-set_var', 'TSRWE', '-endif',
        # IFS cloud cover uses ECMWF local codes (0-6-192..195) as 0-1
        # fractions; write WMO TCDC/LCDC/MCDC/HCDC in %, same layers as AIFS.
        '-if', ':tcc:', '-rpn', '100:*', '-set_var', 'TCDC', '-endif',
        '-if', ':lcc:', '-rpn', '100:*', '-set_var', 'LCDC',
            '-set_lev', 'surface - 800 mb', '-endif',
        '-if', ':mcc:', '-rpn', '100:*', '-set_var', 'MCDC',
            '-set_lev', '800-450 mb', '-endif',
        '-if', ':hcc:', '-rpn', '100:*', '-set_var', 'HCDC',
            '-set_lev', '450 mb - top of atmosphere', '-endif',
        # ECMWF codes msl as plain pressure (0-3-0) at MSL, which AWIPS stores
        # as msl-P; PRMSL (0-3-1) is what AWIPS stores as PMSL, as for GFS.
        '-if', ':PRES:mean sea level:',
            '-set_var', 'PRMSL',
        '-endif',
        # 10fg is a max over the last hour; AWIPS files that as a separate
        # 1-hour parameter, so label it instantaneous to get WGS every hour.
        '-if', ':GUST:',
            '-set_ftime', f'{f_hour} hour fcst',
        '-endif',
        # tcwv (water vapour only) is GFS's PWAT; label it the same so AWIPS
        # shows PW. tcw also counts cloud, rain and snow water. Level type 10
        # is AWIPS EA, as GFS PWAT; GFS's own type 200 is NCEP-local, which
        # AWIPS only knows for NCEP data.
        '-if', ':TCIWV:',
            '-set_var', 'PWAT',
            '-set_lev', 'entire atmosphere',
        '-endif',
        # mucape comes on level 17 (most unstable parcel), reserved in AWIPS's
        # level table; GFS sends most-unstable CAPE as 255-0 mb above ground.
        '-if', ':CAPE:',
            '-set_lev', '255-0 mb above ground',
        '-endif',
        # Snow density (kg m-3) is stored for snow depth below; sd -> WEASD (mm).
        '-if', ':SDEN:', '-rpn', 'sto_2', '-endif',
        '-if', sd, '-rpn', '1000:*', '-set_var', 'WEASD', '-endif',
        '-set_grib_type', 'c3',
        '-grib_out', output_path,
        # Extra copies, written after each record above.
        # BL wind like GFS (30-0 mb above ground): 100 m is ~10-12 mb above the
        # ground, mid-layer; GFS's own 100 m and 0-30 mb speeds differ ~0.2 m/s.
        '-if', ':(UGRD|VGRD):100 m above ground:',
            '-set_lev', '30-0 mb above ground', '-grib_out', output_path,
        '-endif',
        # Snow depth (m) = water equivalent (mm = kg m-2) / density (kg m-3).
        '-if', sd,
            '-rpn', 'rcl_2:/', '-set_var', 'SNOD', '-grib_out', output_path,
        '-endif',
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
            if not sanitize_grib(raw_path, clean_path, f_hour):
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
