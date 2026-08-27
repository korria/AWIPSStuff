#!/usr/bin/env python3
"""
cheatpw -- build precipitable water (PW) percentile guidance products
for station BOI and publish them to AWIPS textdb.
"""

from __future__ import annotations

import calendar
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path

BASEDIR = Path("/localapps/runtime/cheatpw")
MAXOLDHOURS = 30
MAXRUNS = 3

DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
HOURSECS = 60 * 60
DAYSECS = 24 * HOURSECS


# ==============================================================================
# Logging & Utility
# ==============================================================================

def log_date(message: str) -> None:
    """Print message with 'YYYY/MM/DD HH:MM:SS : ' UTC timestamp."""
    stamp = datetime.now(timezone.utc).strftime("%Y/%m/%d %H:%M:%S")
    print(f"{stamp} : {message}")


def get_percentile(val: float, climovals: list[float]) -> int:
    """Calculate the integer percentile (0-100) of a PW value against climo."""
    num = len(climovals)
    for j, c_val in enumerate(climovals):
        if val >= c_val:
            return int(((num - j) / num * 100.0) + 0.5)
    return -9999


# ==============================================================================
# Climate Data
# ==============================================================================

def get_climo() -> tuple[str, list[float], list[str]]:
    """Fetch sorted list of climo PW values for +/-14 days from current time."""
    days_span = 14
    now = time.time()
    _, _, _, _, _, _, _, nyda, _ = time.gmtime(now - (days_span * DAYSECS))
    _, _, _, _, _, _, _, xyda, _ = time.gmtime(now + (days_span * DAYSECS))

    pw_file = BASEDIR / "data" / "pw.dat"
    if not pw_file.exists():
        log_date(f"Error: {pw_file} does not exist")
        return "", [], []

    with open(pw_file, "r") as f:
        alllines = [line.strip() for line in f if line.strip()]

    if not alllines:
        return "", [], []

    lastline = alllines[-1]
    climo: list[str] = []

    for line in alllines[:-1]:
        if ":" not in line:
            continue
        datestr, valstr = line.split(":", 1)
        try:
            yea = int(datestr[0:4])
            mon = int(datestr[4:6])
            day = int(datestr[6:8])
            hou = int(datestr[8:10])
            testtime = calendar.timegm((yea, mon, day, hou, 0, 0))
            _, _, _, _, _, _, _, tyda, _ = time.gmtime(testtime)
        except (ValueError, IndexError):
            continue

        save = False
        if nyda < xyda:
            if nyda <= tyda <= xyda:
                save = True
        else:
            if tyda >= nyda or tyda <= xyda:
                save = True

        if save:
            val_parts = valstr.split(" ", 1)
            try:
                val = float(val_parts[0])
                if val >= 0.00:
                    climo.append(f"{val:4.2f},{datestr}")
            except ValueError:
                continue

    climo.sort(reverse=True)

    values: list[float] = []
    datestrs: list[str] = []
    for cli in climo:
        v_str, d_str = cli.split(",", 1)
        values.append(float(v_str))
        datestrs.append(d_str)

    return lastline, values, datestrs


# ==============================================================================
# Model File Ingestion
# ==============================================================================

def read_mod(allfiles: list[str], num: int, allvdata: dict, allpdata: dict,
             allkeys: list[str]) -> tuple[dict, dict, list[str]]:
    """Read model run files and parse percentile values into data grids."""
    now = time.time()
    for modfile in allfiles:
        mmod = modfile[0:3]
        try:
            myea = int(modfile[3:7])
            mmon = int(modfile[7:9])
            mday = int(modfile[9:11])
            mhou = int(modfile[11:13])
            mmin = int(modfile[13:15])
            mtime = calendar.timegm((myea, mmon, mday, mhou, mmin, 0))
        except (ValueError, IndexError):
            continue

        diff = now - mtime
        if diff > (MAXOLDHOURS * HOURSECS):
            continue

        count = sum(1 for k in allvdata if k[0:3] == mmod)
        if count >= MAXRUNS:
            continue

        modkey = modfile[0:13]
        filepath = BASEDIR / "data" / modfile
        if not filepath.exists():
            continue

        with open(filepath, "r") as f:
            alllines = f.readlines()

        valdata: dict[str, float] = {}
        pctdata: dict[str, int] = {}
        for line in alllines:
            pieces = line.strip().split(",")
            if len(pieces) != 8:
                continue
            key = pieces[0]
            try:
                if num == 10:
                    v = float(pieces[2])
                    p = int(pieces[3])
                elif num == 50:
                    v = float(pieces[4])
                    p = int(pieces[5])
                elif num == 90:
                    v = float(pieces[6])
                    p = int(pieces[7])
                else:
                    continue
            except ValueError:
                continue

            if key not in allkeys:
                allkeys.append(key)

            # Skip 0 or negative values so missing model data does not show up as 0
            if v > 0:
                valdata[key] = v
            if p > 0:
                pctdata[key] = p

        allvdata[modkey] = valdata
        allpdata[modkey] = pctdata

    return allvdata, allpdata, allkeys


# ==============================================================================
# Trend Product Generator
# ==============================================================================

def get_trends(namfiles: list[str], gfsfiles: list[str], raob_time: int,
               stimes: list[int], etimes: list[int], pct_value: int) -> str:
    allvdata: dict = {}
    allpdata: dict = {}
    allkeys: list[str] = []

    allvdata, allpdata, allkeys = read_mod(namfiles, pct_value, allvdata, allpdata, allkeys)
    allvdata, allpdata, allkeys = read_mod(gfsfiles, pct_value, allvdata, allpdata, allkeys)

    allmods = sorted(list(allvdata.keys()), reverse=True)
    allkeys.sort()

    dattext = f"{pct_value:2d}TH PERCENTILE PW FORECAST IN CWA\n\n"

    # 1. Day of Week Header
    dattext += "           "
    for per in range(16):
        stime = stimes[per]
        _, _, _, shou, _, _, swda, _, _ = time.gmtime(stime)
        if shou == 12:
            dattext += f"|  {DAYS[swda]:3s}  "
        elif per == 0 and shou == 0:
            wda = swda - 1 if swda > 0 else 6
            dattext += f" {DAYS[wda]:3s}"
    dattext += "\n"

    # 2. Day of Month Header
    dattext += "           "
    for per in range(16):
        stime = stimes[per]
        _, _, sday, shou, _, _, _, _, _ = time.gmtime(stime)
        spacer = "|" if shou == 12 else " "
        dattext += f"{spacer} {sday:02d}"
    dattext += "\n"

    # 3. DAY/NGT Header
    dattext += "           "
    for per in range(16):
        stime = stimes[per]
        _, _, _, shou, _, _, _, _, _ = time.gmtime(stime)
        dattext += "|DAY" if shou >= 12 else " NGT"
    dattext += "\n"

    # 4. Divider Header
    dattext += "-----------"
    for per in range(16):
        stime = stimes[per]
        _, _, _, shou, _, _, _, _, _ = time.gmtime(stime)
        dattext += "+---" if shou >= 12 else "----"
    dattext += "\n"

    # Select latest complete model runs
    topmods: list[str] = []
    namdone, gfsdone = 0, 0
    for modkey in allmods:
        if modkey[0:3] == "NAM" and namdone == 0:
            topmods.append(modkey)
            if len(allpdata[modkey]) >= 29:
                namdone = 1
        if modkey[0:3] == "GFS" and gfsdone == 0:
            topmods.append(modkey)
            if len(allpdata[modkey]) >= 35:
                gfsdone = 1

    # Render Max Period Table
    for mod in topmods:
        modname = mod[0:3]
        modday = int(mod[9:11])
        modhou = int(mod[11:13])
        dattext += f"{modname:3s} {modday:02d} {modhou:02d}Z:"
        pdata = allpdata[mod]

        for per in range(16):
            stime = stimes[per]
            etime = etimes[per]
            skey = time.strftime("%Y%m%d%H", time.gmtime(stime))
            ekey = time.strftime("%Y%m%d%H", time.gmtime(etime))
            _, _, _, shou, _, _, _, _, _ = time.gmtime(stime)
            spacer = "|" if shou == 12 else " "

            maxp = -1
            for key, p_val in pdata.items():
                if skey <= key <= ekey:
                    if p_val > maxp and p_val > 0:
                        maxp = p_val

            # If no data or maxp <= 0, display spaces instead of 0
            if maxp > 0:
                dattext += f"{spacer}{maxp:3d}"
            else:
                dattext += f"{spacer}   "
        dattext += "\n"

    dattext += "\n"

    # Render Timestep Detail Grid
    dattext += f"            {' ' * 6}{pct_value:2d} PCNTL PW FORECAST   {' ' * 9}PW PERCENTILE   \n"

    dattext += "            "
    for modkey in allmods:
        dattext += f"  {modkey[0:3]:3s}"
    dattext += "  "
    for modkey in allmods:
        dattext += f" {modkey[0:3]:3s}"
    dattext += "\n"

    dattext += "            "
    for modkey in allmods:
        dattext += f"  {modkey[11:13]:2s}Z"
    dattext += "  "
    for modkey in allmods:
        dattext += f" {modkey[11:13]:2s}Z"
    dattext += "\n"

    for key in allkeys:
        kyea, kmon, kday, khou = int(key[0:4]), int(key[4:6]), int(key[6:8]), int(key[8:10])
        ktime = calendar.timegm((kyea, kmon, kday, khou, 0, 0))

        if ktime < raob_time:
            continue

        _, _, fday, fhou, _, _, fwda, _, _ = time.gmtime(ktime)
        dattext += f"{DAYS[fwda]} {fday:02d} {fhou:02d}Z :"

        # Render Model Values (spaces if missing or <= 0)
        for modkey in allmods:
            if key in allvdata[modkey] and allvdata[modkey][key] > 0:
                v = allvdata[modkey][key]
                dattext += f"{v:5.2f}"
            else:
                dattext += "     "

        # Render Model Percentiles (spaces if missing or <= 0)
        dattext += "  "
        for modkey in allmods:
            if key in allpdata[modkey] and allpdata[modkey][key] > 0:
                p = allpdata[modkey][key]
                dattext += f" {p:3d}"
            else:
                dattext += "    "
        dattext += "\n"

    return dattext


# ==============================================================================
# Database Persistence
# ==============================================================================

def save_prod(prodtext: str, prodkey: str) -> int:
    """Save product to file and publish to AWIPS textdb if content changed."""
    data_dir = BASEDIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    filename = data_dir / f"{prodkey}.dat"

    lasttext = ""
    if filename.exists():
        with open(filename, "r") as f:
            lasttext = f.read()

    if lasttext != prodtext:
        log_date(f"Updating product {prodkey}")
        with open(filename, "w") as f:
            f.write(prodtext)

        time.sleep(1)
        os.system(f"textdb -w {prodkey} < {filename}")
        return 1
    return 0


# ==============================================================================
# Main
# ==============================================================================

def main() -> None:
    lastpw, climovals, climodates = get_climo()
    num_climo = len(climovals)

    if not num_climo or not lastpw:
        log_date("Error: Climate database empty or unreadable - Exiting")
        return

    datestr, reststr = lastpw.split(":", 1)
    valstr, _ = reststr.split(" ", 1)
    val = float(valstr)
    pct = get_percentile(val, climovals)

    yea, mon, day, hou = int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8]), int(datestr[8:10])
    ftime = calendar.timegm((yea, mon, day, hou, 0, 0))
    _, _, _, _, _, _, gwda, _, _ = time.gmtime(ftime)

    outtext = "\n"
    outtext += f"BOI Raob on {DAYS[gwda]} {mon:02d}/{day:02d} {hou:02d}Z : {val:4.2f} {pct:3d}%\n\n"

    low = climovals[-1]
    lowdatestr = climodates[-1]
    lowyea, lowmon, lowday, lowhou = int(lowdatestr[0:4]), int(lowdatestr[4:6]), int(lowdatestr[6:8]), int(lowdatestr[8:10])
    lowdate = f"{lowmon:02d}/{lowday:02d}/{lowyea:04d} {lowhou:02d}Z"

    p99 = climovals[int(num_climo * 0.01)]
    p95 = climovals[int(num_climo * 0.05)]
    p90 = climovals[int(num_climo * 0.10)]
    p75 = climovals[int(num_climo * 0.25)]
    p50 = climovals[int(num_climo * 0.50)]
    p25 = climovals[int(num_climo * 0.75)]
    p10 = climovals[int(num_climo * 0.90)]
    p05 = climovals[int(num_climo * 0.95)]
    p01 = climovals[int(num_climo * 0.99)]

    high = climovals[0]
    highdatestr = climodates[0]
    highyea, highmon, highday, highhou = int(highdatestr[0:4]), int(highdatestr[4:6]), int(highdatestr[6:8]), int(highdatestr[8:10])
    highdate = f"{highmon:02d}/{highday:02d}/{highyea:04d} {highhou:02d}Z"

    tot = sum(v for v in climovals if v >= 0.0)
    cnt = sum(1 for v in climovals if v >= 0.0)
    avg = tot / float(cnt) if cnt else 0.0

    outtext += f"BOI Raob PW Climo for 30 days centered on {mon:02d}/{day:02d}\n"
    outtext += f"        high: {high:4.2f} : {highdate}\n"
    outtext += f"  99th pcntl: {p99:4.2f}\n"
    outtext += f"  95th pcntl: {p95:4.2f}\n"
    outtext += f"  90th pcntl: {p90:4.2f}\n"
    outtext += f"  75th pcntl: {p75:4.2f}\n"
    outtext += f"  50th pcntl: {p50:4.2f} : average:{avg:4.2f}\n"
    outtext += f"  25th pcntl: {p25:4.2f}\n"
    outtext += f"  10th pcntl: {p10:4.2f}\n"
    outtext += f"   5th pcntl: {p05:4.2f}\n"
    outtext += f"   1st pcntl: {p01:4.2f}\n"
    outtext += f"         low: {low:4.2f} : {lowdate}\n\n"

    # Find model data files
    nampat = re.compile(r"NAM\d{12}\.dat")
    gfspat = re.compile(r"GFS\d{12}\.dat")
    namfiles: list[str] = []
    gfsfiles: list[str] = []

    data_dir = BASEDIR / "data"
    if data_dir.exists():
        for filename in os.listdir(data_dir):
            if nampat.match(filename):
                namfiles.append(filename)
            elif gfspat.match(filename):
                gfsfiles.append(filename)

    namfiles.sort(reverse=True)
    gfsfiles.sort(reverse=True)

    # Establish 12-hour forecast display periods
    gyea, gmon, gday, ghou, _, _, _, _, _ = time.gmtime()
    bhou = 12 if ghou > 11 else 0
    ftime_b = calendar.timegm((gyea, gmon, gday, bhou, 0, 0))

    stimes: list[int] = []
    etimes: list[int] = []
    for per in range(16):
        stime = ftime_b + (per * (12 * HOURSECS))
        etime = stime + (12 * HOURSECS) - 1
        stimes.append(stime)
        etimes.append(etime)

    # Generate 10th, 50th, and 90th percentile products
    lowtext = get_trends(namfiles, gfsfiles, ftime_b, stimes, etimes, 10)
    medtext = get_trends(namfiles, gfsfiles, ftime_b, stimes, etimes, 50)
    higtext = get_trends(namfiles, gfsfiles, ftime_b, stimes, etimes, 90)

    lowall = outtext + lowtext
    medall = outtext + medtext
    higall = outtext + higtext

    savel = save_prod(lowall, "BOICHTPWL")
    savem = save_prod(medall, "BOICHTPW")
    saveh = save_prod(higall, "BOICHTPWH")

    if (savel + savem + saveh) == 0:
        log_date("No updates to products needed")


if __name__ == "__main__":
    main()