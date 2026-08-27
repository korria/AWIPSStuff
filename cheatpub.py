#!/usr/bin/env python3
"""
cheatpub -- build a guidance cheat sheet for a station and store it in the
AWIPS text database.

Usage:
    cheatpub.py BOI                       normal cron invocation
    cheatpub.py BOI --dry-run             build the file, skip the textdb write
    cheatpub.py BOI --stdout              print it
    cheatpub.py BOI --json                decoded grid as JSON
    cheatpub.py BOI --dump -v             log every row label in every bulletin
    cheatpub.py BOI --now 2026-08-22T18:00:00      pin the clock
    cheatpub.py BOI --from-files a.txt b.txt       no AWIPS needed
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass, field
from datetime import date as Date
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterator

try:
    from zoneinfo import ZoneInfo
    LOCAL_TZ = ZoneInfo("America/Boise")
except Exception:
    ZoneInfo = None
    LOCAL_TZ = None

__version__ = "2.7.0"

log = logging.getLogger("cheatpub")

# ---------------------------------------------------------------------------
# Configuration & Pre-Compiled Regexes
# ---------------------------------------------------------------------------

CCC = "BOI"                                   # node the finished sheet lives on
DEFAULT_DIR = Path("/localapps/runtime/cheatpub")
TEXTDB_CMD = "textdb"
A2DBAUTH = "/awips2/fxa/bin/a2dbauth"
PGUSER = "pguser"
PGDATABASE = "climate"

LABEL_WIDTH = 12
COL_WIDTH = 3
RULE = "-" * 76
MIN_HOUR = 6                                  # slot hour for minimum / night (06Z/12Z)
MAX_HOUR = 18                                 # slot hour for maximum / day (18Z/00Z)
DAY_ABBR = ("SUN", "MON", "TUE", "WED", "THU", "FRI", "SAT")

TEMP_LOW, TEMP_HIGH = -99, 140
POP_LOW, POP_HIGH = 0, 100
SLP_LOW, SLP_HIGH = 0, 999                    # MSLP prints without the thousands digit
CLIMATE_MISSING = 9999
CLIMATE_DAYS = range(-1, 8)                   # day-1 through day+7

# A 12-hour LAMP window is 13 hourly samples; require most of them before
# publishing a rolled-up extreme.
LAMP_MIN_SAMPLES = 9

# Hour of the day, in the issuing office's local time, at which the afternoon
# package starts.  CCF switches over at 18Z, SFT at 20Z (matches the Perl).
CCF_EVENING_CUTOFF = 18
SFT_EVENING_CUTOFF = 20

STATION_NAMES = {
    "BOI": "Boise",
    "PIH": "Pocatello",
    "MSO": "Missoula",
    "EKO": "Elko",
    "PDT": "Pendleton",
    "SJT": "San Angelo",
}

_HEADER_RE = re.compile(
    r"^\s*(?P<station>\w{3,5})\s+(?P<source>.+?)\s+GUIDANCE\s+"
    r"(?P<month>\d{1,2})/(?P<day>\d{1,2})/(?P<year>\d{2,4})\s+"
    r"(?P<cycle>\d{3,4})\s*UTC",
    re.MULTILINE,
)
_NBM_SOURCE_RE = re.compile(r"NBM\s+V(?P<version>\S+)\s+(?P<product>NB[SEXHP])")
_TOKEN_RE = re.compile(r"[^\s|]+")
_MAJOR_SEC_RE = re.compile(r"\s*(X/N|N/X|P12|POP12|MN/MX)(?=\s|\||$)")
_ENSEMBLE_ROW_RE = re.compile(r"\s*(MEAN|AVG|HI|LO|LOW)(?=\s|\||$)")
_WMO_RE = re.compile(r"FPUS\d+\s+K\w{3}\s+(\d{2})(\d{2})(\d{2})?")

# NOTE: the site is a *named* group, so m.groups()[0] is the site and every
# field index below is shifted by one relative to the Perl.
#   groups[0]      site
#   groups[1:6]    line 1 temperatures      (5)
#   groups[6:9]    line 1 PoP characters    (3)
#   groups[9:18]   line 2 temperatures      (9)
#   groups[18:29]  line 2 PoP characters    (11)
_CCF_RE = re.compile(
    r"^(?P<site>[A-Z0-9]{3,5})\s+..\s+(\d{3})/(\d{3})\s+(\d{3})/(\d{3})\s+(\d{3})\s+..(\S)(\S)(\S)\s*\n"
    r"\s*.....\s+(\d{3})/(\d{3})\s+(\d{3})/(\d{3})\s+(\d{3})/(\d{3})\s+(\d{3})/(\d{3})\s+(\d{3})\s+(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)(\S)",
    re.MULTILINE
)
_CCF_TEMPS = slice(1, 6), slice(9, 18)
_CCF_POPS = slice(6, 9), slice(18, 29)

SOURCE_CODES = {
    "GFS MOS": "MAV",
    "GFSX MOS": "MEX",
    "NAM MOS": "MET",
    "ECM MOS": "ECS",
    "ECMWF MOS": "ECS",
    "ECMX MOS": "ECX",
    "ECME MOS": "ECM",
    "ECM ENSEMBLE MOS": "ECM",
    "ECMWF ENSEMBLE MOS": "ECM",
    "ENSEMBLE MOS": "MEN",
}

INDEX_ROWS = ("FHR", "UTC", "HR", "DAY", "VD")


# ===========================================================================
# Small shared helpers
# ===========================================================================

def normalise_cycle(raw: str | int) -> int:
    """'1200' -> 12, '0600' -> 6, '12' -> 12."""
    cycle = int(raw)
    hour = cycle // 100 if cycle > 23 else cycle
    if not 0 <= hour <= 23:
        raise ValueError(f"implausible cycle {raw!r}")
    return hour


def row_values(line: str) -> list[int]:
    """Every integer on a data line, with the row label's own digits removed.

    'P06  10 20' -> [10, 20]   (the 06 in the label is dropped)
    'TMP  33 35' -> [33, 35]   (no digits in the label, nothing dropped)
    """
    fields = line.split(maxsplit=1)
    label = fields[0] if fields else ""
    tokens = re.findall(r"-?\d+", line)
    if any(character.isdigit() for character in label):
        tokens = tokens[1:]
    return [int(token) for token in tokens]


def utc_shift(tz, when: datetime) -> int:
    """Hours the local clock runs *behind* UTC (6 for MDT, 7 for MST)."""
    if tz is None:
        return MAX_HOUR - 12
    offset = when.replace(tzinfo=timezone.utc).astimezone(tz).utcoffset()
    return -int(offset.total_seconds() // 3600)


def resolve_package_base(day_of_month: int, issue_hour: int, now: datetime,
                         evening_cutoff: int) -> datetime | None:
    """Turn a WMO DDHH stamp into the base time of a public forecast package.

    The first 12-hour period of the package is always base + 18h, so:
        morning package  -> base 00Z, first slot 18Z (a maximum)
        evening package  -> base 12Z, first slot 06Z next day (a minimum)
    An evening package that has already crossed 00Z carries the *next* day's
    DD, so it is walked back a day first.
    """
    day = now.date()
    for _ in range(7):
        if day.day == day_of_month:
            break
        day -= timedelta(days=1)
    else:
        log.warning("cannot match day %02d to a date near %s", day_of_month, now.date())
        return None

    if issue_hour < 7:                       # evening package, past 00Z
        return datetime(day.year, day.month, day.day, 12) - timedelta(days=1)
    if issue_hour < evening_cutoff:          # morning / midday package
        return datetime(day.year, day.month, day.day, 0)
    return datetime(day.year, day.month, day.day, 12)     # evening package


# ===========================================================================
# Bulletin decoding
# ===========================================================================

class ParseError(ValueError):
    """The text does not look like a bulletin we can decode."""


@dataclass(frozen=True, slots=True)
class Column:
    """One forecast period: where it sits on the line and when it is valid."""

    index: int
    start: int
    end: int
    valid: datetime


@dataclass(frozen=True, slots=True)
class Bulletin:
    station: str
    source: str
    code: str
    version: str | None
    issued: datetime
    columns: tuple[Column, ...]
    text: str
    lines: tuple[str, ...]

    def labels(self) -> list[str]:
        """Every row label in the bulletin, in order. Used by --dump."""
        out: list[str] = []
        for line in self.lines:
            if (match := re.match(r"\s*([^\s|]+)", line)) and match.group(1) not in out:
                out.append(match.group(1))
        return out

    def raw_row(self, label: str, section: str | None = None) -> dict[datetime, str] | None:
        """Raw field text per valid time, optionally restricted within a major section header."""
        line = self._find_line(label, section=section)
        if line is None:
            return None
        out: dict[datetime, str] = {}
        cut = line.index(label) + len(label)
        for column in self.columns:
            if cut >= len(line):
                break
            field_text = line[cut:column.end].replace("|", " ").strip()
            if field_text:
                out[column.valid] = field_text
            cut = column.end
        return out

    def row(self, label: str, section: str | None = None, *, low: float = float("-inf"),
            high: float = float("inf")) -> dict[datetime, int]:
        raw = self.raw_row(label, section=section)
        if raw is None:
            return {}
        out: dict[datetime, int] = {}
        for valid, text in raw.items():
            if re.fullmatch(r"-?\d+", text) and low <= (value := int(text)) <= high:
                out[valid] = value
        return out

    def has(self, label: str, section: str | None = None) -> bool:
        return self._find_line(label, section=section) is not None

    def _find_line(self, label: str, section: str | None = None) -> str | None:
        pattern = re.compile(r"\s*" + re.escape(label) + r"(?=\s|\||$)")
        if not section:
            return next((line for line in self.lines if pattern.match(line)), None)

        sec_pattern = re.compile(r"\s*" + re.escape(section) + r"(?=\s|\||$)")
        in_section = False

        for line in self.lines:
            if sec_pattern.match(line):
                in_section = True
                continue
            if in_section:
                if _MAJOR_SEC_RE.match(line) and not pattern.match(line):
                    break
                if pattern.match(line):
                    return line
        return None


def fix_year(year: int) -> int:
    if year < 50:
        return year + 2000
    if year < 100:
        return year + 1900
    return year


def _index_row(lines: tuple[str, ...]) -> tuple[str, str]:
    for label in INDEX_ROWS:
        pattern = re.compile(r"\s*" + label + r"(?=\s|\||$)")
        for line in lines:
            if pattern.match(line):
                return label, line
    raise ParseError(f"no index row (looked for {', '.join(INDEX_ROWS)})")


def _walk_hours(issued: datetime, hours: list[int]) -> Iterator[datetime]:
    date = issued.date()
    previous = issued.hour
    first = True
    for hour in hours:
        if (hour < previous) if first else (hour <= previous):
            date += timedelta(days=1)
        yield datetime(date.year, date.month, date.day, hour)
        previous = hour
        first = False


def parse_bulletin(text: str) -> Bulletin:
    text = text.replace("\r", "").replace("\0", "")
    header = _HEADER_RE.search(text)
    if header is None:
        raise ParseError("no '<MODEL> GUIDANCE m/d/y HHMM UTC' header line")

    try:
        hour = normalise_cycle(header["cycle"])
    except ValueError as exc:
        raise ParseError(str(exc)) from None
    issued = datetime(fix_year(int(header["year"])), int(header["month"]),
                      int(header["day"]), hour)

    source = header["source"].strip()
    version: str | None = None
    body = text[header.start():]
    lines = tuple(body.split("\n"))

    if nbm := _NBM_SOURCE_RE.search(source):
        code, version = nbm["product"], nbm["version"]
    elif "ECME" in source or "ECM ENSEMBLE" in source or "ECMWF ENSEMBLE" in source:
        code = "ECM"
    elif source in SOURCE_CODES:
        code = SOURCE_CODES[source]
        # An ECS keyname carrying HI/MEAN/LO *rows* is really the ensemble
        # product.  Match on line starts, not on a substring of the whole text.
        if code == "ECS" and any(_ENSEMBLE_ROW_RE.match(line) for line in lines):
            code = "ECM"
    else:
        raise ParseError(f"unrecognised guidance source {source!r}")

    label, index_line = _index_row(lines)
    tokens = list(_TOKEN_RE.finditer(index_line))[1:]
    if not tokens:
        raise ParseError(f"{label} row has no entries")

    if label == "FHR":
        valids = [issued + timedelta(hours=int(t.group())) for t in tokens]
    else:
        valids = list(_walk_hours(issued, [int(t.group()) for t in tokens]))

    columns = tuple(Column(i, t.start(), t.end(), v)
                    for i, (t, v) in enumerate(zip(tokens, valids)))
    return Bulletin(station=header["station"], source=source, code=code,
                    version=version, issued=issued, columns=columns,
                    text=body, lines=lines)


# ===========================================================================
# Non-Column Bulletin Parsers (SFT, CCF, LAV, FLP, FRH)
# ===========================================================================

@dataclass(frozen=True, slots=True)
class ParsedProduct:
    rows: list[Row] = field(default_factory=list)
    raw_text: str = ""


def _package_slots(base: datetime, count: int) -> list[datetime]:
    """The 12-hourly slot keys of a public package, starting at base + 18h."""
    first = base + timedelta(hours=18)
    return [first + timedelta(hours=12 * i) for i in range(count)]


MONTHS = {name.lower(): number for number, name in enumerate(
    ("Jan", "Feb", "Mar", "Apr", "May", "Jun",
     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), start=1)}

_SFT_DATE_FIELD_RE = re.compile(r"([A-Za-z]{3})\s+(\d{1,2})(?!\d)")
_SFT_MIN_COLUMNS = 3
_SFT_STAMP_RE = re.compile(
    r"^\s*(\d{3,4})\s+(AM|PM)\s+[A-Z]{2,4}\s+[A-Za-z]{3}\s+([A-Za-z]{3})\s+(\d{1,2})\s+(\d{4})",
    re.IGNORECASE)


def _sft_stamp(line: str) -> datetime | None:
    """'230 AM MDT Tue Aug 25 2026' -> local issuance time, for ordering bodies."""
    match = _SFT_STAMP_RE.match(line)
    if not match:
        return None
    month = MONTHS.get(match.group(3).lower())
    if month is None:
        return None
    clock = int(match.group(1))
    hour, minute = clock // 100, clock % 100
    if not (0 <= minute <= 59 and 1 <= hour <= 12):
        return None
    hour = hour % 12 + (12 if match.group(2).upper() == "PM" else 0)
    try:
        return datetime(int(match.group(5)), month, int(match.group(4)), hour, minute)
    except ValueError:
        return None


def _wmo_datetime(text: str, now: datetime) -> datetime | None:
    """Resolve the DDHHMM stamp of a WMO heading to a real time near `now`."""
    match = _WMO_RE.search(text)
    if not match:
        return None
    day_of_month, hour = int(match.group(1)), int(match.group(2))
    if not 0 <= hour <= 23:
        return None
    day = now.date()
    for _ in range(7):
        if day.day == day_of_month:
            return datetime(day.year, day.month, day.day, hour)
        day -= timedelta(days=1)
    return None


def _resolve_month_day(month: int, day: int, now: datetime) -> datetime | None:
    """Attach a year to a bare 'Aug 26', choosing the one nearest to `now`."""
    best = None
    for year in (now.year - 1, now.year, now.year + 1):
        try:
            candidate = datetime(year, month, day)
        except ValueError:
            continue
        if best is None or abs((candidate - now).days) < abs((best - now).days):
            best = candidate
    return best if best is not None and abs((best - now).days) <= 200 else None


def _sft_columns(lines: list[str], now: datetime) -> list[tuple[float, datetime]]:
    """Centre position and date of each forecast column, from the 'Aug 25 Aug 26 ...' row."""
    for line in lines:
        matches = [m for m in _SFT_DATE_FIELD_RE.finditer(line)
                   if m.group(1).lower() in MONTHS]
        if len(matches) < _SFT_MIN_COLUMNS:
            continue
        columns: list[tuple[float, datetime]] = []
        for match in matches:
            when = _resolve_month_day(MONTHS[match.group(1).lower()],
                                      int(match.group(2)), now)
            if when is None:
                return []
            columns.append(((match.start() + match.end()) / 2, when))
        return columns
    return []


def _sft_pair(token: str) -> tuple[int | None, int | None] | None:
    """'69/97' -> (69, 97); '/96' -> (None, 96); 'MM/10' -> (None, 10).

    Returns None when the token is not a numeric pair at all, which is how the
    weather row ('Sunny', 'Rain/Snow') is told apart from the data rows.
    """
    if "/" not in token:
        return None
    left, _, right = token.partition("/")
    out: list[int | None] = []
    for piece in (left, right):
        piece = piece.strip()
        if not piece or piece.upper() == "MM":
            out.append(None)
        elif re.fullmatch(r"-?\d+", piece):
            out.append(int(piece))
        else:
            return None
    return out[0], out[1]


def _sft_data_row(line: str, columns: list[tuple[float, datetime]]
                  ) -> dict[datetime, tuple[int | None, int | None]] | None:
    """Map a low/high or night/day row onto the dated columns.

    Fields are matched to the column whose centre is nearest, because the values
    do not line up character-for-character with the date heading.
    """
    tokens = [(m, _sft_pair(m.group())) for m in re.finditer(r"\S+", line)]
    if not tokens:
        return None
    good = [(m, pair) for m, pair in tokens if pair is not None]
    if len(good) < 2 or len(good) * 2 < len(tokens):
        return None

    out: dict[datetime, tuple[int | None, int | None]] = {}
    for match, pair in good:
        centre = (match.start() + match.end()) / 2
        when = min(columns, key=lambda c: abs(c[0] - centre))[1]
        if when in out:
            log.warning("SFT: two fields fell in the %s column; check alignment",
                        when.strftime("%m/%d"))
        out[when] = pair
    return out


def parse_sft(text: str, site: str, now: datetime, sft_station: str | None = None) -> ParsedProduct:
    """Parse the Tabular State Forecast (SFT).

    Per the product's own legend the rows are:
        temperatures  early morning low / daytime high
        PoP           nighttime 6PM-6AM / daytime 6AM-6PM
    Both halves of a column belong to the same calendar date, so the low and the
    night PoP go to that date's 06Z slot and the high and day PoP to its 18Z
    slot. Column dates are read straight off the 'Aug 25 Aug 26 ...' heading, so
    nothing here depends on guessing the package from the issuance hour.
    """
    text = text.replace("\r", "").replace("\0", "")
    lines = text.splitlines()

    # A single textdb read can carry more than one body (a correction appended
    # to the original, say). Segment on the dated column headings, keep the
    # newest, and confine both the column dates and the station search to it --
    # otherwise one body's dates could be paired with another body's numbers.
    headings = [i for i, line in enumerate(lines) if _sft_columns([line], now)]
    if not headings:
        log.warning("SFT parsing skipped: no dated column heading found")
        return ParsedProduct()

    segments: list[tuple[int, int, datetime | None]] = []
    for n, head in enumerate(headings):
        stop = headings[n + 1] if n + 1 < len(headings) else len(lines)
        stamp = None
        for line in lines[headings[n - 1] if n else 0: head]:
            stamp = _sft_stamp(line) or stamp
        segments.append((head, stop, stamp))

    head, stop, stamp = max(segments, key=lambda s: s[2] or datetime.min)
    if len(segments) > 1:
        log.info("SFT: %d bodies in this product; using the one issued %s",
                 len(segments), stamp.strftime("%m/%d %H%M") if stamp else "first")
    lines = lines[head:stop]
    columns = _sft_columns(lines, now)

    targets: list[str] = []
    if sft_station:
        targets.append(sft_station.strip())
    targets.extend([site.strip(), STATION_NAMES.get(site.strip().upper(), site.strip())])

    unique_targets: list[str] = []
    for t in targets:
        if t and t.lower() not in [u.lower() for u in unique_targets]:
            unique_targets.append(t)

    # The station sits alone on its own indented line above its rows.
    site_idx, matched_target = -1, ""
    for exact in (True, False):
        for t in unique_targets:
            for i, line in enumerate(lines):
                stripped = line.strip()
                hit = (stripped.lower() == t.lower() if exact
                       else re.match(rf"^\s*{re.escape(t)}\b", line, re.IGNORECASE))
                if hit:
                    site_idx, matched_target = i, t
                    break
            if site_idx != -1:
                break
        if site_idx != -1:
            break

    if site_idx == -1:
        log.warning("SFT parsing skipped: station %s (searched %s) not found", site, unique_targets)
        return ParsedProduct()

    # First numeric-pair row under the name is temperatures, the next is PoPs.
    # The weather row in between is rejected by _sft_pair.
    following = [l for l in lines[site_idx + 1: site_idx + 7] if l.strip()]
    rows_found: list[dict[datetime, tuple[int | None, int | None]]] = []
    for line in following:
        parsed = _sft_data_row(line, columns)
        if parsed is not None:
            rows_found.append(parsed)
            log.debug("SFT %s: data row %r", matched_target, line.rstrip())
        elif rows_found:
            break
        if len(rows_found) == 2:
            break

    if not rows_found:
        log.warning("SFT parsing skipped: no low/high row found under %s", matched_target)
        return ParsedProduct()

    temp_vals: dict[datetime, int | str] = {}
    pop_vals: dict[datetime, int | str] = {}
    for target, source in ((temp_vals, rows_found[0]),
                           (pop_vals, rows_found[1] if len(rows_found) > 1 else {})):
        for when, (night, day) in source.items():
            if night is not None:
                target[when.replace(hour=MIN_HOUR)] = night
            if day is not None:
                target[when.replace(hour=MAX_HOUR)] = day

    issued = _wmo_datetime(text, now)
    if issued is None:
        log.warning("SFT: no usable WMO header; labelling with the current time")
        issued = now
    log.debug("SFT %s: %d columns %s..%s, %d temps, %d pops", matched_target,
              len(columns), columns[0][1].strftime("%m/%d"),
              columns[-1][1].strftime("%m/%d"), len(temp_vals), len(pop_vals))

    label = make_label("SFT", site, issued)
    rows = []
    if temp_vals:
        rows.append(Row("SFT", "temps", label, temp_vals))
    if pop_vals:
        rows.append(Row("SFT", "pops", label, pop_vals))
    return ParsedProduct(rows=rows, raw_text="")



def parse_ccf(text: str, site: str, now: datetime) -> ParsedProduct:
    """Parse Coded City Forecast (CCF)."""
    text = text.replace("\r", "").replace("\0", "")
    m = _CCF_RE.search(text)
    if not m or m.group("site") != site:
        log.warning("CCF parsing skipped: site %s pattern not matched", site)
        return ParsedProduct(raw_text=text)

    groups = m.groups()
    temp_fields = list(groups[_CCF_TEMPS[0]]) + list(groups[_CCF_TEMPS[1]])
    pop_fields = list(groups[_CCF_POPS[0]]) + list(groups[_CCF_POPS[1]])

    # ---- anchor to the *product*, not the wall clock -----------------------
    header_m = _WMO_RE.search(text)
    base = None
    if header_m:
        base = resolve_package_base(int(header_m.group(1)), int(header_m.group(2)),
                                    now, CCF_EVENING_CUTOFF)
    if base is None:
        log.warning("CCF: no usable WMO header; falling back to the current date")
        base = now.replace(hour=0, minute=0, second=0, microsecond=0)

    temps: list[int] = []
    for raw in temp_fields:
        if not re.fullmatch(r"\d{3}", raw):
            log.warning("CCF: unexpected temperature field %r; stopping there", raw)
            break
        value = int(raw)
        temps.append(-(value - 900) if value > 900 else value)

    pops: list[int | None] = []
    for character in pop_fields:
        if character == "+":
            pops.append(100)
        elif character == "-":
            pops.append(5)
        elif character.isdigit():
            pops.append(int(character) * 10)
        else:
            pops.append(None)

    temp_vals = dict(zip(_package_slots(base, len(temps)), temps))
    pop_vals = {slot_key: value
                for slot_key, value in zip(_package_slots(base, len(pops)), pops)
                if value is not None}

    label = make_label("CCF", site, base)
    rows = []
    if temp_vals:
        rows.append(Row("CCF", "temps", label, temp_vals))
    if pop_vals:
        rows.append(Row("CCF", "pops", label, pop_vals))
    return ParsedProduct(rows=rows, raw_text=text)


def _combine_pops(p1: int, p2: int) -> int:
    return int(100 - (((100 - p1) / 100.0) * ((100 - p2) / 100.0) * 100.0) + 0.5)


def parse_lamp(text: str, site: str, is_lav: bool = True, tz=LOCAL_TZ) -> ParsedProduct:
    """Parse LAV or FLP hourly LAMP forecast and roll up to 12-hr Min/Max & PoPs.

    Temperature windows follow the *local* day (noon-to-midnight for the max,
    midnight-to-noon for the min) so this is not hard-wired to Mountain time.
    PoP periods stay on the standard 00Z/12Z 12-hour boundaries, which is what
    P12 means everywhere else on the sheet.
    """
    code = "LAV" if is_lav else "FLP"
    text = text.replace("\r", "").replace("\0", "")
    header = _HEADER_RE.search(text)
    if not header:
        log.warning("%s parsing skipped: header missing", code)
        return ParsedProduct(raw_text=text)

    try:
        cycle = normalise_cycle(header["cycle"])
    except ValueError as exc:
        log.warning("%s parsing skipped: %s", code, exc)
        return ParsedProduct(raw_text=text)
    issued = datetime(fix_year(int(header["year"])), int(header["month"]),
                      int(header["day"]), cycle)

    lines = text.splitlines()
    utc_line = next((l for l in lines if re.match(r"\s*UTC\s", l)), None)
    tmp_line = next((l for l in lines if re.match(r"\s*(TMP|TEMP)\s", l)), None)
    pop_line = next((l for l in lines if re.match(r"\s*(P06|POP06)\s", l)), None)

    if not utc_line or not tmp_line:
        log.warning("%s parsing skipped: UTC/TMP lines missing", code)
        return ParsedProduct(raw_text=text)

    utchrs = row_values(utc_line)
    tmpvals = row_values(tmp_line)
    if not utchrs:
        log.warning("%s parsing skipped: empty UTC row", code)
        return ParsedProduct(raw_text=text)

    # One shared time walk, so temps and PoPs can never disagree about the date.
    valids = list(_walk_hours(issued, utchrs))

    temp_vals: dict[datetime, int] = {}
    pop_vals: dict[datetime, int] = {}

    shift = utc_shift(tz, issued)
    maxes: dict[Date, list[int]] = {}
    mins: dict[Date, list[int]] = {}
    for valid, tmp in zip(valids, tmpvals):
        local = valid - timedelta(hours=shift)
        if local.hour >= 12:                                  # local afternoon/evening
            maxes.setdefault(local.date(), []).append(tmp)
        elif local.hour == 0:                                 # closes the previous day
            maxes.setdefault(local.date() - timedelta(days=1), []).append(tmp)
        if local.hour <= 12:                                  # local night/morning
            mins.setdefault(local.date(), []).append(tmp)

    for day, values in maxes.items():
        if len(values) >= LAMP_MIN_SAMPLES:
            temp_vals[datetime(day.year, day.month, day.day, MAX_HOUR)] = max(values)
    for day, values in mins.items():
        if len(values) >= LAMP_MIN_SAMPLES:
            temp_vals[datetime(day.year, day.month, day.day, MIN_HOUR)] = min(values)

    if pop_line:
        pop_tokens = row_values(pop_line)
        first_half: int | None = None
        for valid, pop in zip(valids, pop_tokens):
            if valid.hour in (18, 6):
                first_half = pop
            elif valid.hour in (12, 0) and first_half is not None:
                if (slot_key := slot(valid)) is not None:
                    pop_vals[slot_key] = _combine_pops(first_half, pop)
                first_half = None

    label = make_label(code, site, issued)
    rows = []
    if temp_vals:
        rows.append(Row(code, "temps", label, temp_vals))
    if pop_vals:
        rows.append(Row(code, "pops", label, pop_vals))
    return ParsedProduct(rows=rows, raw_text=text)


def parse_frh(text: str, site: str, max_profile_lines: int = 10) -> ParsedProduct:
    """Extract NAM sounding profile block (FRH)."""
    lines = text.replace("\r", "").replace("\0", "").splitlines()
    out: list[str] = []
    reading = False
    profile_lines = 0
    for line in lines:
        if "OUTPUT FROM" in line:
            out.append(line[:29])
        elif "TTPTTR" in line:
            out.append(line[:33])
        elif re.match(rf"^({site}//|{site}\s+)", line) or re.search(rf"({site}//|{site}\s*$)", line):
            reading = True
            profile_lines = 1
            out.append(line[:33])
        elif reading:
            out.append(line[:33])
            profile_lines += 1
            if profile_lines >= max_profile_lines:
                reading = False
    raw_res = "\n".join(out) + "\n" if out else ""
    return ParsedProduct(rows=[], raw_text=raw_res)


# ===========================================================================
# Sheet assembly and rendering
# ===========================================================================

def slot(valid: datetime) -> datetime | None:
    match valid.hour:
        case 0:
            return (valid - timedelta(days=1)).replace(hour=MAX_HOUR)
        case 18:
            return valid.replace(hour=MAX_HOUR)
        case 12:
            return valid.replace(hour=MIN_HOUR)
        case 6:
            return valid.replace(hour=MIN_HOUR)
        case _:
            return None


@dataclass(frozen=True, slots=True)
class RowSpec:
    code: str
    block: str
    labels: tuple[str, ...]
    low: int = TEMP_LOW
    high: int = TEMP_HIGH
    section: str | None = None


def _mos(code: str) -> tuple[RowSpec, ...]:
    return (RowSpec(code, "temps", ("X/N", "N/X", "MN/MX")),
            RowSpec(code, "pops", ("P12", "POP12"), POP_LOW, POP_HIGH))


# Labels are a *preference list*: the first one present in the bulletin wins.
# NBP labels are the real v5.0 element names (see the NBM station card):
#   TXNP1/P2/P5/P7/P9, G24P1..P9, SLPP1..P9.
ROW_SPECS: dict[str, tuple[RowSpec, ...]] = {
    "MAV": _mos("MAV"),
    "MET": _mos("MET"),
    "MEX": _mos("MEX"),
    "ECS": _mos("ECS"),
    "ECX": _mos("ECX"),
    "FAN": _mos("FAN"),
    "NBE": (RowSpec("NBM", "temps", ("TXN", "X/N", "N/X")),
            RowSpec("NBM", "pops", ("P12",), POP_LOW, POP_HIGH),
            RowSpec("DPT", "expanded", ("DPT",)),
            RowSpec("GST", "expanded", ("GST",), 0, 199)),
    "NBS": (RowSpec("NBS", "temps", ("TXN", "X/N", "N/X", "MN/MX")),
            RowSpec("NBS", "pops", ("P12", "POP12"), POP_LOW, POP_HIGH)),
    "NBX": (RowSpec("NBX", "temps", ("TXN",)),
            RowSpec("NBX", "pops", ("P12",), POP_LOW, POP_HIGH)),
    "NBP": (RowSpec("N90", "temps", ("TXNP9",)),
            RowSpec("N50", "temps", ("TXNP5",)),
            RowSpec("N10", "temps", ("TXNP1",)),
            RowSpec("G90", "expanded", ("G24P9",), 0, 199),
            RowSpec("SLP", "expanded", ("SLPP5",), SLP_LOW, SLP_HIGH)),
    "ECM": (RowSpec("ECH", "temps", ("HI",), section="X/N"),
            RowSpec("ECH", "temps", ("HI",), section="N/X"),
            RowSpec("ECM", "temps", ("MEAN", "AVG"), section="X/N"),
            RowSpec("ECM", "temps", ("MEAN", "AVG"), section="N/X"),
            RowSpec("ECL", "temps", ("LO", "LOW"), section="X/N"),
            RowSpec("ECL", "temps", ("LO", "LOW"), section="N/X"),
            RowSpec("ECH", "pops", ("HI",), POP_LOW, POP_HIGH, section="P12"),
            RowSpec("ECM", "pops", ("MEAN", "AVG"), POP_LOW, POP_HIGH, section="P12"),
            RowSpec("ECL", "pops", ("LO", "LOW"), POP_LOW, POP_HIGH, section="P12")),
    "MEN": (RowSpec("MEH", "temps", ("HI",), section="X/N"),
            RowSpec("MEH", "temps", ("HI",), section="N/X"),
            RowSpec("MEN", "temps", ("MEAN", "AVG"), section="X/N"),
            RowSpec("MEN", "temps", ("MEAN", "AVG"), section="N/X"),
            RowSpec("MEL", "temps", ("LO", "LOW"), section="X/N"),
            RowSpec("MEL", "temps", ("LO", "LOW"), section="N/X"),
            RowSpec("MEH", "pops", ("HI",), POP_LOW, POP_HIGH, section="P12"),
            RowSpec("MEN", "pops", ("MEAN", "AVG"), POP_LOW, POP_HIGH, section="P12"),
            RowSpec("MEL", "pops", ("LO", "LOW"), POP_LOW, POP_HIGH, section="P12")),
}

TEMP_ORDER = ("CCF", "SFT", "NBM", "NBS", "N90", "N50", "N10", "FLP", "LAV", "NBX", "MET",
              "FAN", "MAV", "MEX", "MEH", "MEN", "MEL", "ECS", "ECX", "ECH", "ECM",
              "ECL", "NORMAL", "RECORD", "YEAR")
POP_ORDER = ("CCF", "SFT", "NBM", "NBS", "FLP", "LAV", "NBX", "MET", "FAN", "MAV", "MEX",
             "MEH", "MEN", "MEL", "ECS", "ECX", "ECH", "ECM", "ECL")
EXPANDED_ORDER = ("DPT", "GST", "G90", "SLP")
BLOCK_ORDER = {"temps": TEMP_ORDER, "pops": POP_ORDER, "expanded": EXPANDED_ORDER}
RAW_ORDER = ("CCF", "SFT", "FLP", "LAV", "MET", "FAN", "MAV", "NBE", "NBS", "NBP", "NBX",
             "ECS", "MEX", "MEN", "ECX", "ECM", "FRH")

# Codes that mark where the sheet should start: the official forecast.
ANCHOR_CODES = ("CCF", "SFT")
# Rows that are not model guidance and are exempt from the min/max sanity check.
CLIMATE_CODES = ("NORMAL", "RECORD", "YEAR")


@dataclass(slots=True)
class Row:
    code: str
    block: str
    label: str
    values: dict[datetime, int | str] = field(default_factory=dict)


def make_label(code: str, site: str, issued: datetime | None) -> str:
    if issued is None:
        return code.ljust(LABEL_WIDTH)
    text = f"{code}{site} {issued.day:2d}/{issued.hour:<2d}"
    return text[:LABEL_WIDTH].ljust(LABEL_WIDTH)


def extract(bulletin: Bulletin, site: str,
            specs: tuple[RowSpec, ...] | None = None) -> list[Row]:
    specs = specs if specs is not None else ROW_SPECS.get(bulletin.code, ())
    grouped_values: dict[tuple[str, str], dict[datetime, int | str]] = {}
    spec_order: list[tuple[str, str]] = []

    for spec in specs:
        key = (spec.code, spec.block)
        if key not in grouped_values:
            grouped_values[key] = {}
            spec_order.append(key)

        for label in spec.labels:
            if not bulletin.has(label, section=spec.section):
                continue
            values = bulletin.row(label, section=spec.section,
                                  low=spec.low, high=spec.high)
            dropped = 0
            for valid, value in values.items():
                if (slot_key := slot(valid)) is not None:
                    grouped_values[key][slot_key] = value
                else:
                    dropped += 1
            if dropped:
                log.debug("%s/%s: %d of %d %s columns fell outside 00/06/12/18Z",
                          bulletin.code, label, dropped, len(values), spec.code)
            # First label present wins; do not let a later alias overwrite it.
            break

    rows: list[Row] = []
    for code, block in spec_order:
        values = grouped_values[(code, block)]
        if values:
            rows.append(Row(code, block,
                            make_label(code, site, bulletin.issued),
                            values))
    return rows


def render_row(dates: list[datetime], values: dict[datetime, int | str],
               center: bool = False) -> str:
    parts: list[str] = []
    skip = False
    for i, when in enumerate(dates):
        if skip:
            skip = False
            continue
        prefix = "|" if when.hour == MIN_HOUR else " "
        value = values.get(when)
        if value is None:
            parts.append(prefix + " " * COL_WIDTH)
            continue
        cell = f"{value:>{COL_WIDTH}}"
        if len(cell) > COL_WIDTH:
            cell = " " * COL_WIDTH
        if center and i < len(dates) - 1 and when.day == dates[i + 1].day:
            parts.append(f"{prefix}  {cell}  ")
            skip = True
        else:
            parts.append(prefix + cell)
    return "".join(parts) + "\n"


def render_headers(dates: list[datetime], title: str) -> str:
    days = {d: DAY_ABBR[(d.weekday() + 1) % 7] for d in dates}
    nums: dict[datetime, int | str] = {d: d.day for d in dates}
    match title:
        case "POPs":
            slots = {d: ("MRN" if d.hour == MIN_HOUR else "DAY") for d in dates}
        case "Expanded":
            slots = {d: (" 12" if d.hour == MIN_HOUR else " 00") for d in dates}
        case _:
            slots = {d: ("LOW" if d.hour == MIN_HOUR else " HI") for d in dates}
    return (title.ljust(LABEL_WIDTH) + render_row(dates, days, center=True)
            + " " * LABEL_WIDTH + render_row(dates, nums, center=True)
            + " " * LABEL_WIDTH + render_row(dates, slots))


@dataclass(slots=True)
class Sheet:
    site: str
    now: datetime = field(default_factory=lambda: datetime.now(timezone.utc).replace(tzinfo=None))
    rows: list[Row] = field(default_factory=list)
    raw: list[tuple[str, str]] = field(default_factory=list)

    def _start_cutoff(self) -> datetime:
        """Where the sheet begins: with the official forecast if we have one."""
        # Local midnight of the current day, expressed in UTC.
        if LOCAL_TZ is not None:
            local_now = self.now.replace(tzinfo=timezone.utc).astimezone(LOCAL_TZ)
            local_midnight = local_now.replace(hour=0, minute=0, second=0, microsecond=0)
            floor = local_midnight.astimezone(timezone.utc).replace(tzinfo=None)
        else:
            floor = self.now.replace(hour=0, minute=0, second=0, microsecond=0)

        anchors = [d for row in self.rows if row.code in ANCHOR_CODES for d in row.values]
        if not anchors:
            return floor
        # CCF/SFT trims the stale overnight column, but it is clamped: a stale or
        # mis-decoded package must never be able to hide the current day's
        # guidance from every other row on the sheet.
        return min(max(min(anchors), floor), floor + timedelta(hours=12))

    @property
    def dates(self) -> list[datetime]:
        start_cutoff = self._start_cutoff()
        seen: set[datetime] = set()
        for row in self.rows:
            for d in row.values:
                if d >= start_cutoff:
                    seen.add(d)
        return sorted(seen)

    def block(self, name: str, title: str) -> str:
        rows = [r for r in self.rows if r.block == name and r.values]
        if not rows:
            return ""
        order = BLOCK_ORDER[name]
        rank = {code: i for i, code in enumerate(order)}
        rows.sort(key=lambda r: rank.get(r.code, len(order)))
        dates = self.dates
        body = "".join(r.label + render_row(dates, r.values) for r in rows)
        return render_headers(dates, title) + body

    def render(self) -> str:
        site = self.site
        out = [f"TTAA00 K{site} DDHHMM\n", f"CHT{site}\n\n",
               f"{'-' * 28}CHEAT SHEET FOR {site}{'-' * 29}\n"]
        for name, title in (("temps", "Temps"), ("pops", "POPs"),
                            ("expanded", "Expanded")):
            if text := self.block(name, title):
                out.append(text)
                out.append(RULE + "\n")
        for _, text in self.raw:
            out.append(text if text.endswith("\n") else text + "\n")
            out.append(RULE + "\n")
        return "".join(out)

    def to_json(self) -> dict:
        return {
            "site": self.site,
            "columns": [d.isoformat() + "Z" for d in self.dates],
            "rows": [{"code": r.code, "block": r.block, "label": r.label.strip(),
                      "values": {k.isoformat() + "Z": v
                                 for k, v in sorted(r.values.items())}}
                     for r in self.rows],
        }


def validate_temps(sheet: Sheet) -> int:
    """Warn when a day's low outranks its high -- the signature of a bad decode."""
    problems = 0
    for row in sheet.rows:
        if row.block != "temps" or row.code in CLIMATE_CODES:
            continue
        for when, low in row.values.items():
            if when.hour != MIN_HOUR or not isinstance(low, int):
                continue
            high = row.values.get(when.replace(hour=MAX_HOUR))
            if isinstance(high, int) and low > high:
                log.warning("%s: %s low %d exceeds high %d - check column alignment",
                            row.label.strip(), when.strftime("%m/%d"), low, high)
                problems += 1
    return problems


# ===========================================================================
# Climate database
# ===========================================================================

_RECORDS_SQL = """
SELECT n.day_of_year, n.max_temp_record, n.max_temp_rec_yr1,
       n.min_temp_record, n.min_temp_rec_yr1
  FROM day_climate_norm n
  JOIN cli_sta_setup s ON s.station_id = n.station_id
 WHERE s.station_code = :'code'
   AND n.day_of_year = ANY (string_to_array(:'days', ','));
"""

_NORMALS_SQL = """
SELECT n.day_of_year, n.max_temp_mean, n.min_temp_mean
  FROM day_climate_norm n
  JOIN cli_sta_setup s ON s.station_id = n.station_id
 WHERE s.station_code = :'code'
   AND n.day_of_year = ANY (string_to_array(:'days', ','));
"""


@dataclass(slots=True)
class ClimateRow:
    temps: dict[datetime, int] = field(default_factory=dict)
    years: dict[datetime, str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class PsqlRunner:
    a2dbauth: str = A2DBAUTH
    user: str = PGUSER
    database: str = PGDATABASE
    timeout: int = 120

    def __call__(self, sql: str, variables: dict[str, str]) -> list[str]:
        cmd = [self.a2dbauth, "psql", "-U", self.user, "-d", self.database,
               "-t", "-A", "-F", "|", "--no-psqlrc", "-f", "-"]
        for name, value in variables.items():
            cmd[-2:-2] = ["-v", f"{name}={value}"]
        try:
            proc = subprocess.run(cmd, input=sql, capture_output=True,
                                  text=True, timeout=self.timeout)
        except FileNotFoundError:
            log.error("%s not found; climate rows omitted", self.a2dbauth)
            return []
        except subprocess.TimeoutExpired:
            log.error("climate query timed out after %ds", self.timeout)
            return []
        if proc.returncode != 0:
            log.error("psql exited %d: %s", proc.returncode,
                      proc.stderr.strip() or "(no stderr)")
            return []
        return [line for line in proc.stdout.splitlines() if line.strip()]


def _climate_valid(value: int) -> bool:
    return abs(value) != CLIMATE_MISSING and TEMP_LOW <= value <= TEMP_HIGH


def fetch_climate(kind: str, code: str, now: datetime, *,
                  days: range = CLIMATE_DAYS,
                  runner: PsqlRunner | None = None) -> ClimateRow:
    runner = runner or PsqlRunner()
    targets = {
        (now + timedelta(days=off)).strftime("%m-%d"):
        (now + timedelta(days=off)).replace(hour=0, minute=0, second=0, microsecond=0)
        for off in days
    }
    sql = _RECORDS_SQL if kind == "records" else _NORMALS_SQL
    rows = runner(sql, {"code": code, "days": ",".join(targets)})

    out = ClimateRow()
    for line in rows:
        parts = line.strip().split("|")
        day = targets.get(parts[0])
        if day is None:
            continue
        try:
            numbers = [int(p) for p in parts[1:]]
        except ValueError:
            log.warning("unparseable climate row: %r", line.strip())
            continue

        low_slot = day.replace(hour=MIN_HOUR)
        high_slot = day.replace(hour=MAX_HOUR)
        match kind, numbers:
            case "records", [high, high_year, low, low_year, *_]:
                if _climate_valid(low):
                    out.temps[low_slot] = low
                    if 0 < low_year < CLIMATE_MISSING:
                        out.years[low_slot] = f"{low_year % 1000:03d}"
                if _climate_valid(high):
                    out.temps[high_slot] = high
                    if 0 < high_year < CLIMATE_MISSING:
                        out.years[high_slot] = f"{high_year % 1000:03d}"
            case "normals", [high, low, *_]:
                if _climate_valid(low):
                    out.temps[low_slot] = low
                if _climate_valid(high):
                    out.temps[high_slot] = high
            case _:
                log.warning("unexpected column count in %r", line.strip())
    return out


# ===========================================================================
# Site configuration
# ===========================================================================

EXPECTED = {
    "NBM": "NBE", "NBS": "NBS", "NBP": "NBP", "NBX": "NBX", "MAV": "MAV",
    "MET": "MET", "MEX": "MEX", "ECS": "ECS", "ECX": "ECX",
    "FAN": "FAN", "ECM": "ECM", "MEN": "MEN", "CCF": "CCF",
    "SFT": "SFT", "LAV": "LAV", "FLP": "FLP", "FRH": "FRH"
}

UNSUPPORTED: set[str] = set()


@dataclass(slots=True)
class SiteConfig:
    site: str
    products: dict[str, str] = field(default_factory=dict)
    sft_station: str | None = None
    climate_code: str | None = None
    records: bool = False
    normals: bool = False
    timezone: str | None = None
    source: Path | None = None

    def tzinfo(self):
        if self.timezone and ZoneInfo is not None:
            try:
                return ZoneInfo(self.timezone)
            except Exception:
                log.warning("unknown timezone %r; using the default", self.timezone)
        return LOCAL_TZ


def _derive_nbm_siblings(config: SiteConfig) -> None:
    """NBP/NBS keynames can be inferred from the NBM (NBE) keyname."""
    if "NBM" not in config.products:
        return
    nbm_key = config.products["NBM"]
    for sibling in ("NBP", "NBS"):
        if sibling in config.products:
            continue
        config.products[sibling] = (nbm_key.replace("NBE", sibling)
                                    if "NBE" in nbm_key else f"{sibling}{config.site}")


def load_config(directory: Path, site: str) -> SiteConfig:
    toml_path = directory / f"{site}.toml"
    cfg_path = directory / f"{site}.cfg"
    if toml_path.exists():
        return _load_toml(toml_path, site)
    if cfg_path.exists():
        return _load_legacy(cfg_path, site)
    raise FileNotFoundError(f"no {site}.toml or {site}.cfg in {directory}")


def _load_toml(path: Path, site: str) -> SiteConfig:
    data = tomllib.loads(path.read_text())
    climate = data.get("climate", {})
    config = SiteConfig(
        site=data.get("site", site),
        products={k.upper(): v for k, v in data.get("products", {}).items()},
        sft_station=data.get("sft", {}).get("station_name"),
        climate_code=climate.get("station_code"),
        records=bool(climate.get("records", False)),
        normals=bool(climate.get("normals", False)),
        timezone=data.get("timezone"),
        source=path,
    )
    _derive_nbm_siblings(config)
    _warn_unsupported(config, path)
    return config


def _load_legacy(path: Path, site: str) -> SiteConfig:
    config = SiteConfig(site=site, source=path)
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(":")]
        if len(parts) < 2 or not parts[0]:
            log.warning("%s:%d: ignoring %r", path.name, lineno, raw.strip())
            continue
        code = parts[0].upper()
        match code:
            case "TZ" | "TIMEZONE":
                config.timezone = parts[1]
            case "SFT" if len(parts) >= 3:
                config.products[code] = parts[1]
                config.sft_station = parts[2]
            case "REC" | "NORM" if len(parts) >= 3:
                config.climate_code = parts[2]
                if code == "REC":
                    config.records = True
                else:
                    config.normals = True
            case "REC" | "NORM":
                log.warning("%s:%d: %s needs an id and a station code", path.name, lineno, code)
            case _:
                config.products[code] = parts[1]

    _derive_nbm_siblings(config)
    _warn_unsupported(config, path)
    return config


def _warn_unsupported(config: SiteConfig, path: Path) -> None:
    for code in sorted(set(config.products) & UNSUPPORTED):
        log.warning("%s: %s is not decoded; row omitted", path.name, code)
    for code in sorted(set(config.products) - UNSUPPORTED - set(EXPECTED)):
        log.warning("%s: unknown product %s; row omitted", path.name, code)


# ===========================================================================
# AWIPS text database
# ===========================================================================

@dataclass(frozen=True, slots=True)
class TextDB:
    command: str = TEXTDB_CMD
    timeout: int = 60

    def read(self, key: str, version: int = 0) -> str:
        try:
            proc = subprocess.run([self.command, "-r", f"-{version}:{key}"],
                                  capture_output=True, text=True,
                                  timeout=self.timeout)
        except FileNotFoundError:
            log.error("%s not on PATH", self.command)
            return ""
        except subprocess.TimeoutExpired:
            log.error("timed out reading %s version %d", key, version)
            return ""
        if proc.returncode:
            log.debug("textdb -r -%d:%s exited %d", version, key, proc.returncode)
        return proc.stdout.replace("\r", "").replace("\0", "")

    def read_latest(self, key: str, versions: int = 2) -> str:
        for version in range(versions):
            if text := self.read(key, version):
                if version:
                    log.info("%s: fell back to version %d", key, version)
                return text
        return ""

    def write(self, key: str, text: str) -> bool:
        try:
            proc = subprocess.run([self.command, "-w", key], input=text,
                                  text=True, timeout=self.timeout)
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            log.error("textdb write of %s failed: %s", key, exc)
            return False
        if proc.returncode:
            log.error("textdb -w %s exited %d", key, proc.returncode)
            return False
        return True


# ===========================================================================
# Main collection & Execution
# ===========================================================================

class UtcFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        stamp = datetime.now(timezone.utc).strftime("%Y/%m/%d %H:%M:%S")
        prefix = "" if record.levelno <= logging.INFO else f"{record.levelname}: "
        return f"{stamp} - {prefix}{record.getMessage()}"


def setup_logging(verbose: bool) -> None:
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(UtcFormatter())
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)


def collect(texts: dict[str, str], site: str, now: datetime,
            sft_station: str | None = None, dump: bool = False, tz=LOCAL_TZ) -> Sheet:
    sheet = Sheet(site, now=now)
    decoded: dict[str, Bulletin] = {}
    raw_map: dict[str, str] = {}

    def stash_raw(code: str, text: str) -> None:
        if not text:
            return
        if code in raw_map:
            log.warning("%s: raw text already captured; keeping the first copy", code)
            return
        raw_map[code] = text

    for key, text in texts.items():
        if not text.strip():
            log.warning("%s: empty product", key)
            continue

        res: ParsedProduct | None = None
        match key:
            case "SFT":
                res = parse_sft(text, site, now, sft_station=sft_station)
            case "CCF":
                res = parse_ccf(text, site, now)
            case "LAV":
                res = parse_lamp(text, site, is_lav=True, tz=tz)
            case "FLP":
                res = parse_lamp(text, site, is_lav=False, tz=tz)
            case "FRH":
                res = parse_frh(text, site)

        if res is not None:
            if res.rows:
                sheet.rows.extend(res.rows)
            stash_raw(key, res.raw_text)
            continue

        try:
            bulletin = parse_bulletin(text)
        except ParseError as exc:
            log.warning("%s: %s", key, exc)
            continue

        expected = EXPECTED.get(key)
        if expected and bulletin.code != expected:
            log.warning("%s: expected a %s bulletin, got %s", key, expected, bulletin.code)

        if dump:
            log.info("%s (%s %s): %s", key, bulletin.code,
                     bulletin.issued.strftime("%d/%HZ"), " ".join(bulletin.labels()))

        if bulletin.code in decoded:
            log.warning("%s: a %s bulletin was already decoded; skipping the duplicate",
                        key, bulletin.code)
            continue
        decoded[bulletin.code] = bulletin
        stash_raw(bulletin.code, bulletin.text)
        sheet.rows.extend(extract(bulletin, site))

    ordered = [code for code in RAW_ORDER if code in raw_map]
    ordered += [code for code in raw_map if code not in RAW_ORDER]
    sheet.raw = [(code, raw_map[code]) for code in ordered]
    return sheet


def add_climate(sheet: Sheet, cfg: SiteConfig, now: datetime,
                runner: PsqlRunner | None = None) -> None:
    if not cfg.climate_code:
        return
    if cfg.normals:
        result = fetch_climate("normals", cfg.climate_code, now, runner=runner)
        if result.temps:
            sheet.rows.append(Row("NORMAL", "temps", "NORMAL".ljust(LABEL_WIDTH), dict(result.temps)))
    if cfg.records:
        result = fetch_climate("records", cfg.climate_code, now, runner=runner)
        if result.temps:
            sheet.rows.append(Row("RECORD", "temps", "RECORD".ljust(LABEL_WIDTH), dict(result.temps)))
        if result.years:
            sheet.rows.append(Row("YEAR", "temps", "YEAR".ljust(LABEL_WIDTH), dict(result.years)))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cheatpub.py",
        description="Build a guidance cheat sheet and store it in textdb.")
    parser.add_argument("site", nargs="?", help="station id, e.g. BOI")
    parser.add_argument("--dir", type=Path, default=DEFAULT_DIR,
                        help=f"runtime directory (default {DEFAULT_DIR})")
    parser.add_argument("--out", type=Path, help="write here instead of data/cht<site>.dat")
    parser.add_argument("--stdout", action="store_true", help="print the sheet instead of writing a file")
    parser.add_argument("--json", action="store_true", help="emit the decoded grid as JSON instead of text")
    parser.add_argument("--dry-run", action="store_true", help="build the sheet but do not write it to textdb")
    parser.add_argument("--from-files", type=Path, nargs="+", metavar="FILE",
                        help="decode these bulletin files instead of reading textdb; implies --dry-run")
    parser.add_argument("--no-climate", action="store_true", help="skip the records/normals database lookups")
    parser.add_argument("--dump", action="store_true", help="log the row labels found in each bulletin")
    parser.add_argument("--timezone", help="IANA zone for LAMP day/night windows (default America/Boise)")
    parser.add_argument("--now", help="pin the clock (ISO) for reproducible runs")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)

    if not args.site:
        log.error("No station specified - Exiting")
        return 1
    site = args.site.upper()
    log.info("Starting cheatpub for %s", site)

    now = (datetime.fromisoformat(args.now).replace(tzinfo=None) if args.now
           else datetime.now(timezone.utc).replace(tzinfo=None))

    sft_station = None
    if args.from_files:
        texts = {path.stem.upper(): path.read_text() for path in args.from_files}
        cfg = SiteConfig(site=site)
    else:
        try:
            cfg = load_config(args.dir / "etc", site)
            sft_station = cfg.sft_station
        except FileNotFoundError as exc:
            log.error("%s", exc)
            return 1
        db = TextDB()
        texts = {code: db.read_latest(key) for code, key in cfg.products.items()}

    if args.timezone:
        cfg.timezone = args.timezone

    sheet = collect(texts, site, now, sft_station=sft_station, dump=args.dump,
                    tz=cfg.tzinfo())
    if not args.no_climate and not args.from_files:
        add_climate(sheet, cfg, now)

    if not sheet.rows:
        log.error("no usable guidance decoded - refusing to publish an empty sheet")
        return 1

    validate_temps(sheet)

    payload = (json.dumps(sheet.to_json(), indent=2) + "\n" if args.json else sheet.render())

    if args.stdout or args.from_files:
        sys.stdout.write(payload)
        return 0

    out_path = args.out or args.dir / "data" / f"cht{site.lower()}.dat"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    try:
        tmp_path.write_text(payload)
        tmp_path.replace(out_path)
    except OSError as exc:
        log.error("Cannot write output to %s: %s", out_path, exc)
        return 1

    if args.dry_run or args.json:
        log.info("Dry run - %s written, textdb not updated", out_path)
        return 0

    if not TextDB().write(f"{CCC}CHT{site}".upper(), payload):
        return 1
    log.info("Normal completion of cheatpub for %s", site)
    return 0


if __name__ == "__main__":
    sys.exit(main())
