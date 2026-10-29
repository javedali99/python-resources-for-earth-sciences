#!/usr/bin/env python3
"""
update_readme.py
================

Keeps the README of *Python Resources for Earth Sciences* in sync with its
own library tables. It needs only the Python standard library, so it runs
as-is in GitHub Actions or on any laptop.

Run it from the repository root::

    python scripts/update_readme.py

What it does on every run
-------------------------
1. **Counts** the library rows (lines starting with ``| **[``) inside every
   domain and Hydrology subfield of the README.
2. **Sorts** every library table alphabetically (case-insensitive), so a
   contributor can add a row anywhere in a table.
3. **Rewrites the auto-generated blocks** in the README. Each block sits
   between two HTML comments that GitHub does not display::

       <!-- AUTO:name:START -->  ...generated...  <!-- AUTO:name:END -->

   Blocks: ``stats`` (badges), ``grid`` (domain cards), ``toc``, ``connect`` and one
   ``count:<id>`` badge per domain and subfield.
4. **Redraws the graphics** in ``assets/`` from the same counts:

   * ``icons/<id>.svg``                       icon tiles for every domain,
                                              Hydrology subfield and section
   * ``about-stats.svg``                      About section stat tiles
   * ``collection-map.svg``                   hub-and-spoke map of domains
   * ``connect/<id>.svg``                     contact cards in Connect

   The graphics use one light design (like the banner) with their own
   background, so they look the same in GitHub's light and dark themes.
   A ``<picture>``/``prefers-color-scheme`` switch is avoided on purpose:
   it follows the computer's appearance setting, not the GitHub theme.

Adding a new domain or subfield
-------------------------------
Add it to ``DOMAINS`` below (id, name, colour, icon) and add its anchor
(``<a id="..."></a>``) and table to the README. Counts, cards, charts and
the table of contents then pick it up automatically.
"""

import math
import os
import re
import sys
from html import escape

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
README = os.path.join(ROOT, "README.md")
ASSETS = os.path.join(ROOT, "assets")

# ---------------------------------------------------------------------------
# Brand palette, taken from the banner (navy ink, teal accent).
# ---------------------------------------------------------------------------
THEMES = {
    "light": {
        "panel": "#F5F8FB", "card": "#FFFFFF", "ink": "#0B2545", "muted": "#5B6B82",
        "line": "#DCE4EE", "track": "#E8EEF5", "accent": "#12A39A",
    },
}
FONT = "'Segoe UI', -apple-system, BlinkMacSystemFont, 'Helvetica Neue', Arial, sans-serif"

def rrect(x, y, w, h, r):
    """Rounded rectangle as path data (ICONS entries are stroked paths)."""
    return (f"M{x + r} {y}h{w - 2 * r}a{r} {r} 0 0 1 {r} {r}v{h - 2 * r}a{r} {r} 0 0 1-{r} {r}"
            f"h-{w - 2 * r}a{r} {r} 0 0 1-{r}-{r}v-{h - 2 * r}a{r} {r} 0 0 1 {r}-{r}Z")


# ---------------------------------------------------------------------------
# Line icons on a 24 x 24 grid, drawn in the banner's thin-line style.
# Each value is a list of SVG path "d" strings, stroked (never filled).
# ---------------------------------------------------------------------------
ICONS = {
    "thermometer": [
        "M10 13.6V5a2 2 0 0 1 4 0v8.6a4 4 0 1 1-4 0Z", "M12 9.5v7",
        "M16.5 6h2.5", "M16.5 9h2.5", "M16.5 12h1.5",
    ],
    "globe": [
        "M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18Z", "M12 3c-3 2.6-4.3 5.6-4.3 9s1.3 6.4 4.3 9",
        "M12 3c3 2.6 4.3 5.6 4.3 9s-1.3 6.4-4.3 9", "M3.4 9h17.2", "M3.4 15h17.2",
    ],
    "strata": [
        "M2.5 19.5 9 7.5l3.6 6.2 2.6-3.7 6.3 9.5Z", "M6 16h4.5", "M13.5 17h4",
        "M9 7.5l1.2 3.2",
    ],
    "map": [
        "M3 6.5 9 4l6 2.5L21 4v13.5L15 20l-6-2.5L3 20Z", "M9 4v13.5", "M15 6.5V20",
    ],
    "drop": [
        "M12 3.2c-3.6 4.6-6 8-6 11a6 6 0 0 0 12 0c0-3-2.4-6.4-6-11Z",
        "M9.4 14.6a2.6 2.6 0 0 0 2.6 2.6",
    ],
    "cloud": [
        "M7.5 19h9.8a3.7 3.7 0 0 0 .3-7.4 5.6 5.6 0 0 0-10.8 1.1A3.2 3.2 0 0 0 7.5 19Z",
        "M8.5 5.5v-2", "M4.2 7.3 2.8 5.9", "M12.8 7.3l1.4-1.4",
    ],
    "waves": [
        "M2.5 7.5q2.4-2.6 4.8 0t4.8 0 4.8 0 4.8 0",
        "M2.5 12.5q2.4-2.6 4.8 0t4.8 0 4.8 0 4.8 0",
        "M2.5 17.5q2.4-2.6 4.8 0t4.8 0 4.8 0 4.8 0",
    ],
    "seismo": ["M2 12.5h3.8l2-5 3 10.5 3.2-14 2.8 12 1.8-3.5H22"],
    # --- utility sections -------------------------------------------------
    "info": ["M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18Z", "M12 11.2v5.3", "M12 7.7h.01"],
    "grid": [rrect(4, 4, 6.5, 6.5, 1.6), rrect(13.5, 4, 6.5, 6.5, 1.6),
             rrect(4, 13.5, 6.5, 6.5, 1.6), rrect(13.5, 13.5, 6.5, 6.5, 1.6)],
    "compass": ["M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18Z", "M15.6 8.4 13.5 13.5 8.4 15.6 10.5 10.5Z"],
    "pull": ["M6 3.5a2.2 2.2 0 1 0 0 4.4a2.2 2.2 0 1 0 0-4.4Z", "M6 16.1a2.2 2.2 0 1 0 0 4.4a2.2 2.2 0 1 0 0-4.4Z",
             "M18 16.1a2.2 2.2 0 1 0 0 4.4a2.2 2.2 0 1 0 0-4.4Z", "M6 7.9v8.2",
             "M18 16.1V10a3 3 0 0 0-3-3h-4", "M13 4.8 10.8 7l2.2 2.2"],
    "trend": ["M3.5 20.5h17", "M4.5 16l4.6-4.6 3.6 3 7-7.4", "M15.6 7h4.1v4.1"],
    "chat": ["M5.5 4.5h13a2.5 2.5 0 0 1 2.5 2.5v7.5a2.5 2.5 0 0 1-2.5 2.5h-6.3L7.5 20.5V17h-2A2.5 2.5 0 0 1 3 14.5V7a2.5 2.5 0 0 1 2.5-2.5Z",
             "M8 10.8h.01", "M12 10.8h.01", "M16 10.8h.01"],
    # --- Hydrology subfields ----------------------------------------------
    "database": ["M12 3.5c4.4 0 7.5 1.3 7.5 2.8s-3.1 2.8-7.5 2.8-7.5-1.3-7.5-2.8 3.1-2.8 7.5-2.8Z",
                 "M4.5 6.3v5.7c0 1.5 3.1 2.8 7.5 2.8s7.5-1.3 7.5-2.8V6.3",
                 "M4.5 12v5.7c0 1.5 3.1 2.8 7.5 2.8s7.5-1.3 7.5-2.8V12"],
    "aquifer": ["M3 5.5h18", "M12 8v5.5", "M9.6 11.1 12 13.5l2.4-2.4", "M3 17.5q2.25-2 4.5 0t4.5 0 4.5 0 4.5 0", "M3 21h18"],
    "river": ["M12 21v-6", "M12 15 7.2 10.2", "M12 15l4.8-5.6", "M7.2 10.2 4 5.5", "M7.2 10.2 8.4 4",
              "M16.8 9.4 20 4.5", "M16.8 9.4 14.8 4"],
    "toolbox": [rrect(3, 8.5, 18, 11.5, 1.8), "M8.5 8.5V6.2a1.2 1.2 0 0 1 1.2-1.2h4.6a1.2 1.2 0 0 1 1.2 1.2v2.3",
                "M3 13h18", "M10 13v2.2h4V13"],
    "bell": ["M3 19.5h18", "M4 18.5c3.2 0 4.4-11.5 8-11.5s4.8 11.5 8 11.5", "M12 7v12.5"],
    "series": ["M3.5 4v16.5H20", "M6.5 15.5l4-5 3.2 3 5.8-7"],
    "envelope": [rrect(3, 5, 18, 14, 2.2), "M3.6 6.6 12 13l8.4-6.4"],
}

# ---------------------------------------------------------------------------
# Domains in alphabetical order. "subs" lists Hydrology's subfields as
# (anchor id, heading, short label used in graphics).
# ---------------------------------------------------------------------------
DOMAINS = [
    {"id": "climatology", "name": "Climatology", "short": "Climatology", "color": "#D9822B", "icon": "thermometer"},
    {"id": "earth-system-science", "name": "Earth System Science", "short": "Earth System", "color": "#0E9F8E", "icon": "globe"},
    {"id": "geology", "name": "Geology", "short": "Geology", "color": "#A0603A", "icon": "strata"},
    {"id": "geospatial-analysis-and-mapping", "name": "Geospatial Analysis and Mapping", "short": "Geospatial", "color": "#3F8F4E", "icon": "map"},
    {"id": "hydrology", "name": "Hydrology", "short": "Hydrology", "color": "#1F6FB5", "icon": "drop",
     "subs": [
         ("data-collection", "Data Collection", "Data Collection"),
         ("groundwater-modelling", "Groundwater Modelling", "Groundwater"),
         ("hydrological-modelling", "Hydrological Modelling", "Hydrological Models"),
         ("miscellaneous", "Miscellaneous", "Miscellaneous"),
         ("optimization-uncertainty-statistics", "Optimization, Uncertainty, Statistics", "Optimization & Stats"),
         ("time-series-analysis", "Time Series Analysis", "Time Series"),
     ]},
    {"id": "meteorology", "name": "Meteorology", "short": "Meteorology", "color": "#5B5BD6", "icon": "cloud"},
    {"id": "oceanography", "name": "Oceanography", "short": "Oceanography", "color": "#0B4F8A", "icon": "waves"},
    {"id": "seismology", "name": "Seismology", "short": "Seismology", "color": "#C2413B", "icon": "seismo"},
]

# Icons for the non-domain sections, drawn on a deep navy tile from the banner.
NAVY = "#14365C"
SECTION_ICONS = [
    {"id": "about", "name": "About", "color": NAVY, "icon": "info"},
    {"id": "browse", "name": "Browse by Domain", "color": NAVY, "icon": "grid"},
    {"id": "at-a-glance", "name": "The Collection at a Glance", "color": NAVY, "icon": "compass"},
    {"id": "contributing", "name": "Contributing", "color": NAVY, "icon": "pull"},
    {"id": "star-history", "name": "Star History", "color": NAVY, "icon": "trend"},
    {"id": "connect", "name": "Connect", "color": NAVY, "icon": "chat"},
]
# Contact cards for the Connect section. Brand logos (LinkedIn, X) are from
# simple-icons (CC0); "line" logos reuse the stroked ICONS set.
CONNECT = [
    {"id": "website", "name": "Website", "label": "WEBSITE", "handle": "javedali.net", "url": "https://javedali.net",
     "from": "#0B2545", "to": "#12A39A", "logo": ("line", "globe")},
    {"id": "linkedin", "name": "LinkedIn", "label": "LINKEDIN", "handle": "in/javedali18", "url": "https://www.linkedin.com/in/javedali18",
     "from": "#004182", "to": "#2E8BE6", "logo": ("fill",
        "M20.447 20.452h-3.554v-5.569c0-1.328-.027-3.037-1.852-3.037-1.853 0-2.136 1.445-2.136 2.939v5.667H9.351V9h3.414v1.561h.046c.477-.9 1.637-1.85 3.37-1.85 3.601 0 4.267 2.37 4.267 5.455v6.286zM5.337 7.433c-1.144 0-2.063-.926-2.063-2.065 0-1.138.92-2.063 2.063-2.063 1.14 0 2.064.925 2.064 2.063 0 1.139-.925 2.065-2.064 2.065zm1.782 13.019H3.555V9h3.564v11.452z")},
    {"id": "x", "name": "X", "label": "X (TWITTER)", "handle": "@javedali99", "url": "https://twitter.com/javedali99",
     "from": "#0F1419", "to": "#5B5BD6", "logo": ("fill",
        "M18.901 1.153h3.68l-8.04 9.19L24 22.846h-7.406l-5.8-7.584-6.638 7.584H.474l8.6-9.83L0 1.154h7.594l5.243 6.932ZM17.61 20.644h2.039L6.486 3.24H4.298Z")},
    {"id": "email", "name": "Email", "label": "EMAIL", "handle": "javedali28@gmail.com", "url": "mailto:javedali28@gmail.com",
     "from": "#C5221F", "to": "#F29900", "logo": ("line", "envelope")},
]

# Hydrology subfield icons, keyed by anchor id, drawn on the Hydrology blue.
SUB_ICONS = {
    "data-collection": "database", "groundwater-modelling": "aquifer",
    "hydrological-modelling": "river", "miscellaneous": "toolbox",
    "optimization-uncertainty-statistics": "bell", "time-series-analysis": "series",
}

ROW_RE = re.compile(r"^\| \*\*\[([^\]]+)\]")
ANCHOR_RE = re.compile(r'^<a id="([a-z0-9-]+)"></a>\s*$')


# ===========================================================================
# Counting and sorting
# ===========================================================================
def sort_key(line):
    """Alphabetical key for a table row: the library name, case-insensitive."""
    name = ROW_RE.match(line).group(1)
    return re.sub(r"^[^A-Za-z0-9]+", "", name).lower()


def sort_tables(lines):
    """Sort every contiguous run of library rows A to Z. Returns new lines."""
    out, run = [], []
    for line in lines + [""]:
        if ROW_RE.match(line):
            run.append(line)
            continue
        if run:
            out.extend(sorted(run, key=sort_key))
            run = []
        out.append(line)
    return out[:-1]


def count_rows(lines):
    """Count library rows per anchor inside the LIBRARIES region.

    Rows are attributed to the most recent ``<a id>`` anchor, so a Hydrology
    row lands on its subfield. Domain totals add their subfields together.
    """
    known = {d["id"] for d in DOMAINS} | {s[0] for d in DOMAINS for s in d.get("subs", [])}
    counts, current, inside = {k: 0 for k in known}, None, False
    for line in lines:
        if "<!-- LIBRARIES:START -->" in line:
            inside = True
        elif "<!-- LIBRARIES:END -->" in line:
            inside = False
        if not inside:
            continue
        m = ANCHOR_RE.match(line)
        if m and m.group(1) in known:
            current = m.group(1)
        elif ROW_RE.match(line) and current:
            counts[current] += 1
    for d in DOMAINS:
        d["count"] = counts[d["id"]] + sum(counts[s[0]] for s in d.get("subs", []))
        d["sub_counts"] = [(sid, name, short, counts[sid]) for sid, name, short in d.get("subs", [])]
    return sum(d["count"] for d in DOMAINS)


# ===========================================================================
# README blocks
# ===========================================================================
def pct(n, total):
    """Share as a whole percent; shows '<1%' instead of a misleading 0%."""
    if not total:
        return "0%"
    share = n / total
    return "<1%" if 0 < share < 0.005 else f"{share:.0%}"


def libs(n):
    """'1 library' or 'n libraries'."""
    return f"{n} library" if n == 1 else f"{n} libraries"


def badge(message, color, label="", style="flat-square"):
    """shields.io static badge URL. Dashes and spaces are escaped for shields."""
    esc = lambda s: s.replace("-", "--").replace("_", "__").replace(" ", "%20")
    return f"https://img.shields.io/badge/{esc(label)}-{esc(message)}-{color.lstrip('#')}?style={style}"


def block_stats(total):
    """Second badge row under the banner: totals plus Python and last commit."""
    repo = "javedali99/python-resources-for-earth-sciences"
    return (
        '<p align="center">\n'
        f'  <img src="{badge(str(total), "#0B2545", "libraries")}" alt="{libs(total)}">\n'
        f'  <img src="{badge(str(len(DOMAINS)), "#12A39A", "domains")}" alt="{len(DOMAINS)} domains">\n'
        '  <img src="https://img.shields.io/badge/made%20with-Python-3776AB?style=flat-square&logo=python&logoColor=white" alt="Made with Python">\n'
        f'  <a href="https://github.com/{repo}/commits"><img src="https://img.shields.io/github/last-commit/{repo}?style=flat-square&color=5B6B82" alt="Last commit"></a>\n'
        "</p>"
    )


# Width of the README content area on github.com (repo page, desktop).
GRID_WIDTH = 830
GRID_COLS = 4


def grid_column_widths():
    """Pixel width for each column of the domain grid.

    GitHub renders tables only as wide as their content, so the grid would
    not fill the page on its own. Each column gets the width its longest
    name needs on one line (about 8.6 px per character in GitHub's bold
    16 px font, plus 26 px cell padding). Columns are otherwise equal, and
    together they add up to GRID_WIDTH, so long names like "Geospatial
    Analysis and Mapping" get a wider column instead of wrapping.
    """
    cols = [DOMAINS[c::GRID_COLS] for c in range(GRID_COLS)]
    need = [max(len(d["name"]) * 8.6 + 26 for d in col) for col in cols]
    # "Water-filling": columns are as equal as possible; a column whose
    # name needs more than the equal share keeps its own width, and the
    # remaining space is split evenly between the others.
    fixed = set()
    while True:
        free = [c for c in range(GRID_COLS) if c not in fixed]
        share = (GRID_WIDTH - sum(need[c] for c in fixed)) / max(1, len(free))
        wide = [c for c in free if need[c] > share]
        if not wide:
            break
        fixed.update(wide)
    return [int(need[c]) if c in fixed else int(share) for c in range(GRID_COLS)]


def block_grid():
    """Clickable domain cards: icon tile, name and count badge, 4 per row.

    Names use non-breaking spaces so each stays on one line, and every
    cell carries its column's pixel width so the grid spans the page.
    """
    widths = grid_column_widths()
    rows = ["<table>"]
    for i in range(0, len(DOMAINS), GRID_COLS):
        rows.append("  <tr>")
        for c, d in enumerate(DOMAINS[i:i + GRID_COLS]):
            name = d["name"].replace(" ", "&nbsp;")
            rows.append(
                f'    <td align="center" width="{widths[c]}" valign="top">\n'
                f'      <a href="#{d["id"]}"><img src="assets/icons/{d["id"]}.svg" width="56" height="56" alt="{d["name"]} icon"></a><br>\n'
                f'      <a href="#{d["id"]}"><b>{name}</b></a><br>\n'
                f'      <img src="{badge(libs(d["count"]), d["color"])}" alt="{libs(d["count"])}">\n'
                f"    </td>"
            )
        rows.append("  </tr>")
    rows.append("</table>")
    return "\n".join(rows)


def block_toc():
    """Collapsible full table of contents with counts."""
    lines = ["<details>", "<summary><b>Full table of contents</b></summary>", "<br>", ""]
    for d in DOMAINS:
        lines.append(f'- [{d["name"]}](#{d["id"]}) <sub>{libs(d["count"])}</sub>')
        for sid, name, _, n in d["sub_counts"]:
            lines.append(f"  - [{name}](#{sid}) <sub>{libs(n)}</sub>")
    lines += ["", "</details>"]
    return "\n".join(lines)


def block_count(color, n):
    """Count badge shown under a domain or subfield heading."""
    return f'<img src="{badge(libs(n), color)}" alt="{libs(n)}">'


def replace_block(text, name, content):
    """Replace the content between AUTO:name markers (markers are kept)."""
    pattern = re.compile(rf"(<!-- AUTO:{re.escape(name)}:START -->)(.*?)(<!-- AUTO:{re.escape(name)}:END -->)", re.S)
    if not pattern.search(text):
        print(f"warning: block AUTO:{name} not found", file=sys.stderr)
        return text
    return pattern.sub(lambda m: f"{m.group(1)}\n{content}\n{m.group(3)}", text)


# ===========================================================================
# SVG helpers
# ===========================================================================
def text_width(s, size, bold=False):
    """Rough width of system sans-serif text, used to size chips."""
    return len(s) * size * (0.58 if bold else 0.54)


def icon_group(icon, x, y, size, stroke, width=1.7):
    """Stroked line icon scaled from the 24 x 24 grid to `size` at (x, y)."""
    s = size / 24
    paths = "".join(f'<path d="{d}"/>' for d in ICONS[icon])
    return (f'<g transform="translate({x:.1f} {y:.1f}) scale({s:.4f})" fill="none" stroke="{stroke}" '
            f'stroke-width="{width / s:.2f}" stroke-linecap="round" stroke-linejoin="round">{paths}</g>')


def icon_tile(d, x, y, size, radius=None):
    """Coloured rounded tile with a white line icon (used in every graphic)."""
    r = radius if radius is not None else size * 0.28
    pad = size * 0.2
    return (f'<rect x="{x}" y="{y}" width="{size}" height="{size}" rx="{r:.1f}" fill="{d["color"]}"/>'
            + icon_group(d["icon"], x + pad, y + pad, size - 2 * pad, "#FFFFFF"))


def svg_doc(w, h, body, title):
    """Wrap SVG body markup with a root element, title and default font."""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" '
            f'role="img" aria-label="{escape(title)}" font-family="{FONT}">'
            f"<title>{escape(title)}</title>{body}</svg>\n")


def header(t, eyebrow, headline, x=56, y=64):
    """Eyebrow + headline + short accent rule, echoing the banner layout."""
    return (f'<text x="{x}" y="{y}" fill="{t["accent"]}" font-size="13" font-weight="700" letter-spacing="3">{escape(eyebrow)}</text>'
            f'<text x="{x}" y="{y + 40}" fill="{t["ink"]}" font-size="28" font-weight="700">{escape(headline)}</text>'
            f'<rect x="{x}" y="{y + 58}" width="44" height="3" rx="1.5" fill="{t["accent"]}"/>')


# ===========================================================================
# Graphics
# ===========================================================================
def svg_icon_file(d):
    """Standalone 64 x 64 icon tile for the domain grid and headings."""
    return svg_doc(64, 64, icon_tile(d, 0, 0, 64, 16), f'{d["name"]} icon')


def svg_about_stats(t, total):
    """Four stat tiles for the About section."""
    hydro = next(d for d in DOMAINS if d.get("subs"))
    tiles = [
        (str(total), "LIBRARIES", "Open-source Python packages"),
        (str(len(DOMAINS)), "DOMAINS", "From climate to seismology"),
        (str(len(hydro["subs"])), "HYDROLOGY SUBFIELDS", "Data, models and statistics"),
        ("A–Z", "ORDER", "Every list sorted alphabetically"),
    ]
    W, H, gap = 1200, 196, 24
    tw = (W - 3 * gap) / 4
    body = []
    for i, (num, label, cap) in enumerate(tiles):
        x = i * (tw + gap)
        body.append(
            f'<rect x="{x + 1:.1f}" y="1" width="{tw - 2:.1f}" height="{H - 2}" rx="18" fill="{t["card"]}" stroke="{t["line"]}" stroke-width="1.5"/>'
            f'<rect x="{x + 30:.1f}" y="34" width="36" height="3" rx="1.5" fill="{t["accent"]}"/>'
            f'<text x="{x + 30:.1f}" y="104" fill="{t["ink"]}" font-size="54" font-weight="700" letter-spacing="-1">{escape(num)}</text>'
            f'<text x="{x + 30:.1f}" y="140" fill="{t["accent"]}" font-size="12.5" font-weight="700" letter-spacing="2.4">{label}</text>'
            f'<text x="{x + 30:.1f}" y="166" fill="{t["muted"]}" font-size="15">{escape(cap)}</text>'
        )
    return svg_doc(W, H, "".join(body), f"{total} libraries across {len(DOMAINS)} domains")


def svg_collection_map(t, total):
    """Hub-and-spoke map: globe in the centre, a card per domain on each
    side with icon, count and share bar, Hydrology subfield chips below."""
    W, H = 1200, 690
    cx, cy, R = 600, 368, 112
    cw, ch, top, step = 330, 82, 178, 98
    body = [f'<rect x="1" y="1" width="{W - 2}" height="{H - 2}" rx="24" fill="{t["panel"]}" stroke="{t["line"]}" stroke-width="1.5"/>',
            header(t, "THE COLLECTION AT A GLANCE", f"{total} libraries across {len(DOMAINS)} Earth science domains")]

    # --- Globe, drawn like the banner: navy outline, teal graticule -------
    grat = []
    for k in (0.42, 0.78):  # meridians as ellipses of shrinking width
        grat.append(f'<ellipse cx="{cx}" cy="{cy}" rx="{R * k:.1f}" ry="{R}"/>')
    grat.append(f'<line x1="{cx}" y1="{cy - R}" x2="{cx}" y2="{cy + R}"/>')
    for dy in (-56, 0, 56):  # parallels
        half = math.sqrt(R * R - dy * dy)
        grat.append(f'<line x1="{cx - half:.1f}" y1="{cy + dy}" x2="{cx + half:.1f}" y2="{cy + dy}"/>')
    body.append(
        f'<circle cx="{cx}" cy="{cy}" r="{R + 22}" fill="{t["accent"]}" opacity="0.07"/>'
        f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="{t["card"]}"/>'
        f'<g transform="rotate(-18 {cx} {cy})" fill="none" stroke="{t["accent"]}" stroke-width="1.3" opacity="0.85">{"".join(grat)}</g>'
        f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" stroke="{t["ink"]}" stroke-width="2.2"/>'
        f'<circle cx="{cx + 34}" cy="{cy + 22}" r="6" fill="{t["accent"]}"/>'
        f'<text x="{cx}" y="{cy + R + 46}" text-anchor="middle" fill="{t["ink"]}" font-size="15" font-weight="700">Python for Earth Sciences</text>'
    )

    # --- Domain cards and connectors ------------------------------------------
    left_angles, right_angles = [212, 194, 166, 148], [-32, -14, 14, 32]
    biggest = max(d["count"] for d in DOMAINS)
    for i, d in enumerate(DOMAINS):
        side = 0 if i < 4 else 1
        row = i % 4
        x = 56 if side == 0 else W - 56 - cw
        y = top + row * step
        ang = math.radians((left_angles if side == 0 else right_angles)[row])
        gx, gy = cx + R * math.cos(ang), cy + R * math.sin(ang)
        sx = x + cw if side == 0 else x
        sy = y + ch / 2
        c1 = sx + (70 if side == 0 else -70)
        c2 = gx + (-70 if side == 0 else 70)
        body.append(
            f'<path d="M{sx:.1f} {sy:.1f} C{c1:.1f} {sy:.1f} {c2:.1f} {gy:.1f} {gx:.1f} {gy:.1f}" fill="none" '
            f'stroke="{d["color"]}" stroke-width="1.6" opacity="0.7"/>'
            f'<circle cx="{gx:.1f}" cy="{gy:.1f}" r="4" fill="{d["color"]}"/>'
            f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="3.5" fill="{t["card"]}" stroke="{d["color"]}" stroke-width="1.6"/>'
        )
        share = d["count"] / total if total else 0
        bar_w = cw - 104
        body.append(
            f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}" rx="16" fill="{t["card"]}" stroke="{t["line"]}" stroke-width="1.4"/>'
            + icon_tile(d, x + 18, y + 17, 48, 13)
            + f'<text x="{x + 82}" y="{y + 33}" fill="{t["ink"]}" font-size="15.5" font-weight="700">{escape(d["short"])}</text>'
            f'<text x="{x + cw - 18}" y="{y + 33}" text-anchor="end" fill="{t["ink"]}" font-size="15.5" font-weight="700">{d["count"]}</text>'
            f'<rect x="{x + 82}" y="{y + 46}" width="{bar_w}" height="7" rx="3.5" fill="{t["track"]}"/>'
            f'<rect x="{x + 82}" y="{y + 46}" width="{max(7, bar_w * d["count"] / biggest):.1f}" height="7" rx="3.5" fill="{d["color"]}"/>'
            f'<text x="{x + 82}" y="{y + 70}" fill="{t["muted"]}" font-size="12.5">{escape(pct(d["count"], total))} of the collection</text>'
        )

    # --- Hydrology subfield chips ---------------------------------------------
    hydro = next(d for d in DOMAINS if d.get("subs"))
    y0 = 598
    body.append(f'<text x="56" y="{y0}" fill="{t["accent"]}" font-size="12" font-weight="700" letter-spacing="2.4">HYDROLOGY SUBFIELDS</text>')
    chips, x = [], 56
    for _, _, short, n in hydro["sub_counts"]:
        w = text_width(short, 13.5) + text_width(str(n), 13.5, True) + 46
        chips.append((x, short, n, w))
        x += w + 10
    for cx0, short, n, w in chips:
        num_fill = hydro["color"] if t is THEMES["light"] else t["accent"]
        body.append(
            f'<rect x="{cx0:.1f}" y="{y0 + 16}" width="{w:.1f}" height="36" rx="18" fill="{t["card"]}" stroke="{t["line"]}" stroke-width="1.3"/>'
            f'<circle cx="{cx0 + 18:.1f}" cy="{y0 + 34}" r="4.5" fill="{hydro["color"]}"/>'
            f'<text x="{cx0 + 30:.1f}" y="{y0 + 39}" fill="{t["ink"]}" font-size="13.5">{escape(short)} '
            f'<tspan font-weight="700" fill="{num_fill}">{n}</tspan></text>'
        )
    return svg_doc(W, H, "".join(body), f"Map of {total} libraries across {len(DOMAINS)} domains")


def svg_connect_card(c):
    """Gradient contact card: logo tile, small label, handle and an arrow.

    Drawn at 640 x 176 and shown at 320 x 88 in the README, so it stays
    sharp on high-density screens. Colours work on light and dark themes.
    """
    W, H = 640, 176
    gid = f"g-{c['id']}"
    kind, val = c["logo"]
    tile = 112
    tx, ty = 32, (H - tile) / 2
    if kind == "line":
        logo = icon_group(val, tx + 24, ty + 24, tile - 48, "#FFFFFF", 4.6)
    else:  # filled brand mark on the 24 x 24 grid
        sc = (tile - 56) / 24
        logo = f'<path transform="translate({tx + 28:.1f} {ty + 28:.1f}) scale({sc:.3f})" d="{val}" fill="#FFFFFF"/>'
    size = 36 if len(c["handle"]) <= 14 else 30
    body = (
        f'<defs><linearGradient id="{gid}" x1="0" y1="0" x2="1" y2="1">'
        f'<stop offset="0" stop-color="{c["from"]}"/><stop offset="1" stop-color="{c["to"]}"/></linearGradient>'
        f'<clipPath id="c-{c["id"]}"><rect width="{W}" height="{H}" rx="34"/></clipPath></defs>'
        f'<g clip-path="url(#c-{c["id"]})">'
        f'<rect width="{W}" height="{H}" fill="url(#{gid})"/>'
        # soft decorative rings in the corner, echoing the banner's globe lines
        f'<circle cx="{W - 40}" cy="-10" r="120" fill="none" stroke="#FFFFFF" stroke-opacity="0.14" stroke-width="2"/>'
        f'<circle cx="{W - 40}" cy="-10" r="170" fill="none" stroke="#FFFFFF" stroke-opacity="0.09" stroke-width="2"/>'
        f'<circle cx="{W - 40}" cy="-10" r="220" fill="none" stroke="#FFFFFF" stroke-opacity="0.05" stroke-width="2"/>'
        f'<rect x="{tx}" y="{ty}" width="{tile}" height="{tile}" rx="30" fill="#FFFFFF" fill-opacity="0.18"/>'
        f'<rect x="{tx + 1}" y="{ty + 1}" width="{tile - 2}" height="{tile - 2}" rx="29" fill="none" stroke="#FFFFFF" stroke-opacity="0.35" stroke-width="2"/>'
        + logo +
        f'<text x="176" y="74" fill="#FFFFFF" fill-opacity="0.8" font-size="20" font-weight="700" letter-spacing="4">{escape(c["label"])}</text>'
        f'<text x="176" y="122" fill="#FFFFFF" font-size="{size}" font-weight="700">{escape(c["handle"])}</text>'
        f'<circle cx="{W - 62}" cy="{H / 2}" r="26" fill="#FFFFFF" fill-opacity="0.18"/>'
        f'<path d="M{W - 71} {H / 2 + 9} l18 -18 M{W - 67} {H / 2 - 9} h14 v14" fill="none" stroke="#FFFFFF" stroke-width="3.5" stroke-linecap="round" stroke-linejoin="round"/>'
        "</g>"
    )
    return svg_doc(W, H, body, f'{c["name"]}: {c["handle"]}')


def block_connect():
    """Two rows of two clickable contact cards."""
    def link(c):
        return (f'<a href="{c["url"]}"><img src="assets/connect/{c["id"]}.svg" width="320" '
                f'alt="{escape(c["name"])}: {escape(c["handle"])}"></a>')
    rows = [CONNECT[i:i + 2] for i in range(0, len(CONNECT), 2)]
    inner = "\n  <br>\n".join("\n".join("  " + link(c) for c in row) for row in rows)
    return f'<p align="center">\n{inner}\n</p>'


def write(path, content):
    """Write a file only when its content changed (keeps git diffs clean)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and open(path, encoding="utf-8").read() == content:
        return False
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return True


# ===========================================================================
# Main
# ===========================================================================
def main(readme=README, assets=ASSETS):
    """Sort tables, recount, rewrite README blocks and redraw graphics."""
    with open(readme, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    lines = sort_tables(lines)
    total = count_rows(lines)
    text = "\n".join(lines)

    text = replace_block(text, "stats", block_stats(total))
    text = replace_block(text, "grid", block_grid())
    text = replace_block(text, "toc", block_toc())
    for d in DOMAINS:
        text = replace_block(text, f"count:{d['id']}", block_count(d["color"], d["count"]))
        for sid, _, _, n in d["sub_counts"]:
            text = replace_block(text, f"count:{sid}", block_count(d["color"], n))

    text = replace_block(text, "connect", block_connect())

    changed = [write(readme, text)]
    for d in DOMAINS:
        changed.append(write(os.path.join(assets, "icons", f"{d['id']}.svg"), svg_icon_file(d)))
        for sid, name, _ in d.get("subs", []):
            sub = {"name": name, "color": d["color"], "icon": SUB_ICONS[sid]}
            changed.append(write(os.path.join(assets, "icons", f"{sid}.svg"), svg_icon_file(sub)))
    for c in CONNECT:
        changed.append(write(os.path.join(assets, "connect", f"{c['id']}.svg"), svg_connect_card(c)))
    for sec in SECTION_ICONS:
        changed.append(write(os.path.join(assets, "icons", f"{sec['id']}.svg"), svg_icon_file(sec)))
    t = THEMES["light"]
    changed.append(write(os.path.join(assets, "about-stats.svg"), svg_about_stats(t, total)))
    changed.append(write(os.path.join(assets, "collection-map.svg"), svg_collection_map(t, total)))

    print(f"{total} libraries in {len(DOMAINS)} domains; "
          + ", ".join(f"{d['short']} {d['count']}" for d in DOMAINS))
    print("files updated:", sum(changed))


if __name__ == "__main__":
    main()
