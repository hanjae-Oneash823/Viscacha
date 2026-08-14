"""DOSSIER — interactive homepage for browsing every hit that enters
master_surveyor. Self-contained HTML (embedded data + vanilla JS), same
"no server needed" philosophy as the dossiers themselves -- open the file
directly, or `python -m http.server` from outputs/dossier/ for network access.

Usage:
    python generate_index.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from dossier import palette
from dossier.config import MASTER_SURVEYOR_GROUPS, OUT_DIR
from dossier.data import dossier_filename, load_hits
from dossier import manifest as manifest_mod
from dossier.render import THEME_TOGGLE_HTML, THEME_TOGGLE_JS, _font_face_css
from master_surveyor.m0_select import representative_row as _representative_row

# Plottable fields -- defined once here, used both to build the axis/color
# <select> options (grouped into <optgroup>s by the 3rd element) and (as
# JSON) for the JS-side field-label lookup, so the label text can't drift
# between the two.
NUMERIC_FIELDS = [
    ("delta_pp", "Δ usage Control→AD (pp)", "Usage & significance"),
    ("neg_log10_padj", "-log10(padj)", "Usage & significance"),
    ("padj", "padj (raw)", "Usage & significance"),
    ("control_pct", "Control usage (%)", "Usage & significance"),
    ("ad_pct", "AD usage (%)", "Usage & significance"),
    ("n_alts", "# ranked alts", "Usage & significance"),
    ("protein_length_diff", "Protein length change (aa)", "Protein"),
    ("pct_identity", "Protein identity to canonical", "Protein"),
    ("changed_aa_fraction", "Fraction of protein changed", "Protein"),
    ("n_domains", "Pfam domains (canonical)", "Protein"),
    ("chembl_max_phase", "ChEMBL max phase", "Drug evidence"),
    ("n_drugs", "ChEMBL drugs (n)", "Drug evidence"),
    ("chembl_bioactive_compounds", "ChEMBL bioactive compounds (n)", "Drug evidence"),
    ("chembl_best_pchembl", "ChEMBL best pChEMBL", "Drug evidence"),
    ("ot_max_phase", "Open Targets max phase", "Drug evidence"),
    ("ot_n_drugs", "Open Targets drugs (n)", "Drug evidence"),
    ("ot_trials_total", "Clinical trials on record (n)", "Drug evidence"),
    ("ot_trials_terminated", "Trials stopped early (n)", "Drug evidence"),
    ("pharos_n_ligands", "Pharos ligands (n)", "Drug evidence"),
    ("pharos_n_drugs", "Pharos drugs (n)", "Drug evidence"),
    ("dgidb_interactions", "DGIdb interactions (n)", "Drug evidence"),
]
CATEGORICAL_FIELDS = [
    ("master_group", "Group", "General"),
    ("cell_type", "Cell type", "General"),
    ("protein_change_type", "Protein change type", "Protein"),
    ("change_bucket", "Length change", "Protein"),
    ("has_evidence", "Has drug evidence", "Drug evidence"),
    ("has_failed_trial", "Has a stopped clinical trial", "Drug evidence"),
    ("pharos_tdl", "Pharos target development level", "Drug evidence"),
]

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def build_manifest() -> list[dict]:
    df = load_hits()
    df = df[df["master_group"].isin(MASTER_SURVEYOR_GROUPS)]
    n_alts_map = manifest_mod.n_alts_by_hit(df)

    records = []
    for (gene, ct, hit_enst), g in df.groupby(["gene_name", "cell_type", "hit_ENST_ID"]):
        r = _representative_row(g)
        record = manifest_mod.format_record(
            r, n_alts_map[(gene, ct, hit_enst)], dossier_filename(gene, ct, hit_enst),
        )
        records.append(record)
    return sorted(records, key=lambda r: (r["gene"], r["cell_type"]))


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

_TF_COLOR = "#4a3aa7"   # violet -- deliberately NOT blue/red, which mean
_DR_COLOR = "#1baf7a"   # aqua   -- Control/AD everywhere else in this tool

_CSS = f"""
:root {{ color-scheme: light; }}
.viz-root {{
  /* Pure white page background -- index-specific (dossiers keep the shared
     off-white palette.PAGE_PLANE); --surface-1 stays a hair off-white so
     cards still read as a distinct layer even before the border kicks in. */
  --surface-1: {palette.SURFACE}; --page-plane: #ffffff;
  --text-primary: {palette.TEXT_PRIMARY}; --text-secondary: {palette.TEXT_SECONDARY};
  --muted: {palette.MUTED}; --gridline: {palette.GRIDLINE}; --baseline: {palette.BASELINE};
  --border: {palette.BORDER}; --bar-black: #000000;
  --accent-control: {palette.ACCENT_CONTROL_HEX}; --accent-ad: {palette.ACCENT_AD_HEX};
  --status-good: #0ca30c; --status-critical: #d03b3b; --glow: none;
  --tf: {_TF_COLOR}; --dr: {_DR_COLOR};
}}
/* Dark theme is a deliberate, distinct world (pure black + neon), not an
   inverted tint of the light one -- shared with render.py via palette.py
   so the dossiers and this homepage can't drift into two different "dark
   modes" when navigating between them. */
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{
    {palette.NEON_DARK_CSS_VARS} --tf: {palette.NEON_DARK_TF}; --dr: {palette.NEON_DARK_DR};
  }}
}}
:root[data-theme="dark"] .viz-root {{
  {palette.NEON_DARK_CSS_VARS} --tf: {palette.NEON_DARK_TF}; --dr: {palette.NEON_DARK_DR};
}}
.viz-root {{
  font-family: "Tamzen", ui-monospace, monospace;
  background: var(--page-plane); color: var(--text-primary);
  max-width: 1400px; margin: 0 auto; padding: 24px 20px 64px;
  line-height: 1.5;
}}
.card {{ background: var(--surface-1); border: 1px solid var(--border); border-radius: 0; padding: 20px 24px; margin-bottom: 20px; }}
h1 {{ font-size: 28px; margin: 0 0 4px; letter-spacing: 0.02em; text-shadow: var(--glow); }}
h2 {{ font-size: 15px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--text-secondary); margin: 0 0 14px; font-weight: 700; }}
.subtle {{ color: var(--text-secondary); font-size: 13px; }}
.muted {{ color: var(--muted); }}
svg text {{ font-family: "Tamzen", ui-monospace, monospace; fill: var(--text-secondary); }}
.tip {{
  position: fixed; pointer-events: none; z-index: 50; display: none;
  background: #ffffff; color: #0b0b0b; border: 1px solid rgba(11,11,11,0.15);
  box-shadow: 0 2px 10px rgba(11,11,11,0.18);
  font-size: 11px; padding: 6px 9px; border-radius: 0; line-height: 1.4;
  white-space: pre; max-width: 260px;
}}
.theme-toggle {{
  position: fixed; top: 16px; right: 16px; z-index: 60;
  width: 34px; height: 34px; border-radius: 0;
  border: 1px solid var(--border); background: var(--surface-1); color: var(--text-primary);
  font-size: 16px; cursor: pointer; display: flex; align-items: center; justify-content: center;
}}
.theme-toggle:hover {{ border-color: var(--accent-control); }}
.stat-row {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 8px; }}
.stat-tile {{ background: var(--page-plane); border: 1px solid var(--border); border-radius: 0; padding: 7px 12px; display: flex; align-items: baseline; gap: 7px; }}
.stat-tile .value {{ font-size: 17px; font-weight: 700; text-shadow: var(--glow); }}
.stat-tile .label {{ font-size: 10px; color: var(--text-secondary); text-transform: uppercase; letter-spacing: 0.05em; }}
.charts-row {{ display: grid; grid-template-columns: 2fr 1fr; gap: 20px; align-items: start; }}
@media (max-width: 900px) {{ .charts-row {{ grid-template-columns: 1fr; }} }}
.legend-row {{ display: flex; gap: 16px; font-size: 12px; color: var(--text-secondary); margin-top: 8px; }}
.legend-row .swatch {{ display: inline-block; width: 10px; height: 10px; border-radius: 0; margin-right: 5px; vertical-align: -1px; }}
.controls {{ display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 16px; }}
.controls input[type=text] {{
  font-family: inherit; font-size: 13px; padding: 7px 10px; border: 1px solid var(--border);
  border-radius: 0; background: var(--surface-1); color: var(--text-primary); min-width: 220px;
}}
.controls select {{
  font-family: inherit; font-size: 12px; padding: 6px 8px; border: 1px solid var(--border);
  border-radius: 0; background: var(--surface-1); color: var(--text-primary);
}}
.controls .count {{ font-size: 12px; color: var(--text-secondary); margin-left: auto; }}
.control-panel {{ display: flex; flex-wrap: wrap; gap: 16px; margin-bottom: 16px; }}
.control-group {{ width: 220px; }}
.control-label {{ font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-secondary); font-weight: 700; margin-bottom: 6px; }}
.option-list {{ max-height: 168px; overflow-y: auto; border: 1px solid var(--border); background: var(--surface-1); }}
.option-group-label {{
  padding: 5px 10px 3px; font-size: 10px; text-transform: uppercase; letter-spacing: 0.04em;
  color: var(--muted); position: sticky; top: 0; background: var(--surface-1);
}}
.option-row {{ padding: 5px 10px; font-size: 12px; cursor: pointer; display: flex; align-items: center; gap: 8px; }}
.option-row:hover {{ background: var(--page-plane); }}
.option-row.selected {{ background: var(--accent-control); color: #ffffff; }}
.option-row .swatch-preview {{ width: 30px; height: 9px; flex-shrink: 0; border: 1px solid var(--border); }}
table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
th, td {{ text-align: left; padding: 8px 10px; border-bottom: 1px solid var(--gridline); white-space: nowrap; }}
th {{ color: var(--text-secondary); font-weight: 700; cursor: pointer; user-select: none; position: sticky; top: 0; background: var(--surface-1); }}
th.sorted::after {{ content: " \\25BE"; }}
tbody tr {{ cursor: pointer; }}
tbody tr:hover {{ background: var(--page-plane); }}
.badge {{ display: inline-block; padding: 2px 8px; border-radius: 0; border: 1px solid var(--border); font-size: 11px; text-shadow: var(--glow); }}
.badge.tf {{ color: var(--tf); border-color: var(--tf); }}
.badge.dr {{ color: var(--dr); border-color: var(--dr); }}
.badge.shorter {{ color: var(--accent-ad); border-color: var(--accent-ad); }}
.badge.longer {{ color: var(--accent-control); border-color: var(--accent-control); }}
.badge.same {{ color: var(--text-primary); border-color: var(--text-primary); }}
.dot {{ display: inline-block; width: 8px; height: 8px; border-radius: 0; margin-right: 3px; }}
.table-wrap {{ max-height: 640px; overflow: auto; border: 1px solid var(--border); border-radius: 0; }}
"""

_JS = r"""
(function() {
  var DATA = window.__HITS__;
  var FIELDS = window.__FIELDS__;   // {key: label}, numeric + categorical combined
  var NUMERIC_KEYS = window.__NUMERIC_KEYS__;
  var CAT_PALETTE = window.__CAT_PALETTE__;
  var OTHER_GRAY = '#898781';
  var tip = document.getElementById('tip');

  function showTip(el, e) {
    tip.textContent = el.getAttribute('data-tip');
    tip.style.display = 'block';
    var x = e.clientX + 14, y = e.clientY + 14;
    if (x + 260 > window.innerWidth) x = e.clientX - 274;
    tip.style.left = x + 'px'; tip.style.top = y + 'px';
  }
  document.addEventListener('pointermove', function(e) {
    var t = e.target.closest('[data-tip]');
    if (!t) { tip.style.display = 'none'; return; }
    showTip(t, e);
  });
  document.addEventListener('pointerdown', function(e) {
    var t = e.target.closest('[data-file]');
    if (t) { window.open(t.getAttribute('data-file'), '_blank'); return; }
    var tt = e.target.closest('[data-tip]');
    tip.style.display = tt ? 'block' : 'none';
  });

  // ---- shared filter/search/sort/plot state ----
  var state = {
    q: '', cellType: '', group: '', bucket: '', evidence: '', sortKey: 'padj', sortDir: 1,
    xKey: 'delta_pp', yKey: 'neg_log10_padj', colorKey: 'master_group', colorScheme: 'sequential'
  };

  function applyFilters() {
    return DATA.filter(function(r) {
      if (state.q && r.gene.toLowerCase().indexOf(state.q) === -1) return false;
      if (state.cellType && r.cell_type !== state.cellType) return false;
      if (state.group && r.master_group !== state.group) return false;
      if (state.bucket && r.change_bucket !== state.bucket) return false;
      if (state.evidence === 'yes' && !r.has_evidence) return false;
      if (state.evidence === 'no' && r.has_evidence) return false;
      return true;
    });
  }

  function sortRows(rows) {
    var k = state.sortKey, d = state.sortDir;
    return rows.slice().sort(function(a, b) {
      var av = a[k], bv = b[k];
      if (av === null || av === undefined) av = d > 0 ? Infinity : -Infinity;
      if (bv === null || bv === undefined) bv = d > 0 ? Infinity : -Infinity;
      if (av < bv) return -1 * d;
      if (av > bv) return 1 * d;
      return 0;
    });
  }

  function esc(s) { var d = document.createElement('div'); d.textContent = String(s); return d.innerHTML; }
  function isNumeric(key) { return NUMERIC_KEYS.indexOf(key) !== -1; }

  function fmtTick(v) {
    if (v === null || v === undefined || isNaN(v)) return '';
    if (Math.abs(v) < 0.001 && v !== 0) return v.toExponential(1);
    return (Math.round(v * 100) / 100).toString();
  }

  function groupBadge(g) {
    return g === 'trial_failure_candidate'
      ? '<span class="badge tf">trial failure</span>'
      : '<span class="badge dr">drug repurposing</span>';
  }

  // ---- table ----
  function renderTable() {
    var rows = sortRows(applyFilters());
    document.getElementById('count').textContent = rows.length + ' / ' + DATA.length + ' hits';
    var body = rows.map(function(r) {
      var evChips = Object.keys(r.evidence).filter(function(k) { return r.evidence[k]; }).join(', ') || '&mdash;';
      var stopFlag = r.has_failed_trial
        ? ' <span class="badge" style="background:var(--status-critical); color:#fff;" title="' +
          esc(r.ot_trials_terminated + '/' + r.ot_trials_total + ' trials stopped early') + '">stopped trial</span>'
        : '';
      return '<tr data-file="' + esc(r.file) + '">' +
        '<td><strong>' + esc(r.gene) + '</strong></td>' +
        '<td>' + esc(r.cell_type) + '</td>' +
        '<td>' + groupBadge(r.master_group) + '</td>' +
        '<td><span class="badge ' + r.change_bucket + '">' + esc(r.protein_change_type) + '</span></td>' +
        '<td>' + (r.padj !== null ? r.padj.toExponential(2) : '&mdash;') + '</td>' +
        '<td style="color:var(--accent-control);">' + r.control_pct.toFixed(1) + '%</td>' +
        '<td style="color:var(--accent-ad);">' + r.ad_pct.toFixed(1) + '%</td>' +
        '<td>' + (r.delta_pp > 0 ? '+' : '') + r.delta_pp.toFixed(1) + ' pp</td>' +
        '<td>' + r.n_alts + '</td>' +
        '<td style="white-space:normal; max-width:220px;">' + evChips + stopFlag + '</td>' +
        '</tr>';
    }).join('');
    document.getElementById('tbody').innerHTML = body;
  }

  // ---- color helpers ----
  function lerp(a, b, t) { return a + (b - a) * t; }
  function rgbHex(rgb) {
    return '#' + rgb.map(function(v) {
      return Math.max(0, Math.min(255, Math.round(v))).toString(16).padStart(2, '0');
    }).join('');
  }
  function seqColor(t) {           // light -> dark blue, t in [0,1]
    return rgbHex([lerp(205, 13, t), lerp(226, 54, t), lerp(251, 107, t)]);
  }
  function divColor(t) {           // blue <-> red through neutral gray, t in [-1,1]
    var target = t < 0 ? [42, 120, 214] : [227, 73, 72];
    var a = Math.min(Math.abs(t), 1);
    return rgbHex([lerp(240, target[0], a), lerp(239, target[1], a), lerp(236, target[2], a)]);
  }

  // Scientific colormaps (published-source stops, 9-point piecewise-linear
  // approximation -- plenty for a UI color scale). Unlike the single-hue
  // seq/div ramps above, these are a deliberate exception to "sequential =
  // one hue": the user asked for them by name. viridis/plasma/magma/mako/
  // rocket are perceptually-uniform, monotonic-lightness maps (colorblind-
  // safe by construction). turbo is different in kind -- a rainbow-family
  // map (Google's fix for jet's worst artifacts) chosen for punchy visual
  // range rather than monotonic lightness; it does NOT carry the same
  // colorblind-safety guarantee as the other five.
  var COLORMAPS = {
    viridis: [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[110,206,88],[253,231,37]],
    plasma:  [[13,8,135],[71,3,159],[115,1,168],[156,23,158],[189,55,134],[216,87,107],[237,121,83],[250,159,58],[240,249,33]],
    magma:   [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],[229,80,100],[251,135,97],[254,194,135],[252,253,191]],
    mako:    [[11,5,5],[35,22,58],[35,58,89],[26,89,107],[24,116,110],[41,145,101],[93,171,89],[168,196,90],[222,222,142]],
    rocket:  [[3,5,26],[43,20,74],[95,27,109],[144,33,105],[192,47,84],[224,79,63],[241,127,67],[247,179,110],[250,234,190]],
    turbo:   [[48,18,59],[70,107,227],[39,168,224],[55,207,157],[151,222,73],[233,196,49],[247,131,44],[211,60,32],[122,4,3]]
  };
  function colormapColor(name, t) {
    var stops = COLORMAPS[name];
    t = Math.max(0, Math.min(1, t));
    var pos = t * (stops.length - 1);
    var i = Math.min(Math.floor(pos), stops.length - 2);
    var frac = pos - i;
    return rgbHex([
      lerp(stops[i][0], stops[i + 1][0], frac),
      lerp(stops[i][1], stops[i + 1][1], frac),
      lerp(stops[i][2], stops[i + 1][2], frac)
    ]);
  }

  var FIXED_COLORS = {
    master_group: { trial_failure_candidate: 'var(--tf)', drug_repurposing_candidate: 'var(--dr)' },
    change_bucket: { shorter: 'var(--accent-ad)', longer: 'var(--accent-control)', same: 'var(--bar-black)' }
  };
  var DISPLAY_NAMES = {
    trial_failure_candidate: 'trial failure', drug_repurposing_candidate: 'drug repurposing',
    same: 'same length'
  };
  function displayName(v) { return DISPLAY_NAMES[v] || v; }

  function categoricalColorMap(rows, key) {
    var counts = {};
    rows.forEach(function(r) { var v = String(r[key]); counts[v] = (counts[v] || 0) + 1; });
    var sorted = Object.keys(counts).sort(function(a, b) { return counts[b] - counts[a]; });
    // Fixed hue order, never cycled: past 8 distinct values, the rest fold
    // into one shared "Other" grey instead of generating a 9th hue.
    var realSlots = sorted.length <= CAT_PALETTE.length ? CAT_PALETTE.length : CAT_PALETTE.length - 1;
    var map = {};
    sorted.forEach(function(v, i) {
      map[v] = i < realSlots ? CAT_PALETTE[i] : OTHER_GRAY;
    });
    return map;
  }

  // Dispatches on state.colorScheme. Sequential/diverging take t already
  // normalized to their own domain ([0,1] / [-1,1] respectively); the
  // colormaps just want [0,1], same as sequential.
  function scaleColor(scheme, t01, tSigned) {
    if (scheme === 'diverging') return divColor(tSigned);
    if (COLORMAPS[scheme]) return colormapColor(scheme, t01);
    return seqColor(t01);
  }

  function colorFor(r, rows, key) {
    if (!key || key === 'none') return 'var(--bar-black)';
    if (key === 'has_evidence') return r.has_evidence ? 'var(--status-good)' : 'var(--status-critical)';
    if (key === 'has_failed_trial') return r.has_failed_trial ? 'var(--status-critical)' : 'var(--gridline)';
    if (FIXED_COLORS[key]) return FIXED_COLORS[key][r[key]] || OTHER_GRAY;
    if (isNumeric(key)) {
      var vals = rows.map(function(x) { return x[key]; }).filter(function(v) { return v !== null && v !== undefined; });
      if (!vals.length) return OTHER_GRAY;
      var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
      var v = r[key];
      if (v === null || v === undefined) return 'var(--gridline)';
      var extent = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
      return scaleColor(state.colorScheme, hi > lo ? (v - lo) / (hi - lo) : 0.5, v / extent);
    }
    return categoricalColorMap(rows, key)[String(r[key])] || OTHER_GRAY;
  }

  function renderLegend(rows) {
    var key = state.colorKey;
    if (!key || key === 'none') return '';
    if (key === 'has_evidence') {
      return '<div class="legend-row"><span><span class="swatch" style="background:var(--status-good);"></span>has evidence</span>' +
        '<span><span class="swatch" style="background:var(--status-critical);"></span>no evidence</span></div>';
    }
    if (key === 'has_failed_trial') {
      return '<div class="legend-row"><span><span class="swatch" style="background:var(--status-critical);"></span>trial stopped early</span>' +
        '<span><span class="swatch" style="background:var(--gridline);"></span>no stopped trial</span></div>';
    }
    if (FIXED_COLORS[key]) {
      return '<div class="legend-row">' + Object.keys(FIXED_COLORS[key]).map(function(k) {
        return '<span><span class="swatch" style="background:' + FIXED_COLORS[key][k] + ';"></span>' + esc(displayName(k)) + '</span>';
      }).join('') + '</div>';
    }
    if (isNumeric(key)) {
      var stops = state.colorScheme === 'diverging'
        ? [scaleColor('diverging', 0, -1), scaleColor('diverging', 0.5, 0), scaleColor('diverging', 1, 1)]
        : [0, 0.25, 0.5, 0.75, 1].map(function(t) { return scaleColor(state.colorScheme, t, t); });
      var grad = 'linear-gradient(to right,' + stops.join(',') + ')';
      return '<div style="margin-top:8px;"><div style="width:180px; height:10px; border-radius: 0; ' +
        'background:' + grad + '; border:1px solid var(--border);"></div>' +
        '<div class="subtle" style="margin-top:3px;">' + esc(FIELDS[key]) + ' (low &rarr; high)</div></div>';
    }
    var map = categoricalColorMap(rows, key);
    return '<div class="legend-row" style="flex-wrap:wrap;">' + Object.keys(map).map(function(k) {
      return '<span><span class="swatch" style="background:' + map[k] + ';"></span>' + esc(k) + '</span>';
    }).join('') + '</div>';
  }

  // ---- scatter plot ----
  function fmtAxisTick(v) { return Math.round(v).toString(); }

  function renderPlot() {
    var rows = applyFilters();
    var xKey = state.xKey, yKey = state.yKey;
    var width = 700, height = 380, padL = 56, padR = 16, padT = 14, padB = 40;
    var plotW = width - padL - padR, plotH = height - padT - padB;

    if (!rows.length) {
      document.getElementById('scatter').innerHTML = '<div class="subtle">no hits match current filters</div>';
      document.getElementById('scatter-legend').innerHTML = '';
      return;
    }

    var xs = rows.map(function(r) { return r[xKey]; }).filter(function(v) { return v !== null && v !== undefined; });
    var ys = rows.map(function(r) { return r[yKey]; }).filter(function(v) { return v !== null && v !== undefined; });
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
    // Always include the origin so x=0/y=0 are guaranteed to fall inside
    // the plotted range, not just when the data happens to straddle it.
    xMin = Math.min(xMin, 0); xMax = Math.max(xMax, 0);
    yMin = Math.min(yMin, 0); yMax = Math.max(yMax, 0);
    var xPad = (xMax - xMin) * 0.08 || 1, yPad = (yMax - yMin) * 0.08 || 1;
    xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

    function xOf(v) { return padL + (v - xMin) / (xMax - xMin) * plotW; }
    function yOf(v) { return padT + plotH * (1 - (v - yMin) / (yMax - yMin)); }

    var grid = '';
    for (var i = 0; i <= 4; i++) {
      var frac = i / 4;
      var gy = padT + plotH * (1 - frac), gx = padL + plotW * frac;
      grid += '<line x1="' + padL + '" y1="' + gy.toFixed(1) + '" x2="' + (padL + plotW) + '" y2="' + gy.toFixed(1) +
        '" stroke="var(--gridline)" stroke-width="1"/>';
      grid += '<text x="' + (padL - 8) + '" y="' + (gy + 4).toFixed(1) + '" font-size="10" text-anchor="end">' +
        fmtAxisTick(yMin + (yMax - yMin) * frac) + '</text>';
      grid += '<text x="' + gx.toFixed(1) + '" y="' + (height - 14) + '" font-size="10" text-anchor="middle">' +
        fmtAxisTick(xMin + (xMax - xMin) * frac) + '</text>';
    }
    // x=0 / y=0 reference axes -- distinct from the hairline gridlines above.
    var zx = xOf(0), zy = yOf(0);
    grid += '<line x1="' + zx.toFixed(1) + '" y1="' + padT + '" x2="' + zx.toFixed(1) + '" y2="' + (padT + plotH) +
      '" stroke="var(--baseline)" stroke-width="1.5"/>';
    grid += '<line x1="' + padL + '" y1="' + zy.toFixed(1) + '" x2="' + (padL + plotW) + '" y2="' + zy.toFixed(1) +
      '" stroke="var(--baseline)" stroke-width="1.5"/>';

    var dots = rows.map(function(r) {
      var xv = r[xKey], yv = r[yKey];
      if (xv === null || xv === undefined || yv === null || yv === undefined) return '';
      var x = xOf(xv), y = yOf(yv);
      var color = colorFor(r, rows, state.colorKey);
      var tip = r.gene + ' / ' + r.cell_type + '&#10;' + r.hit_transcript + '&#10;' +
        FIELDS[xKey] + ': ' + fmtTick(xv) + '&#10;' + FIELDS[yKey] + ': ' + fmtTick(yv) +
        '&#10;(click to open dossier)';
      return '<circle cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) + '" r="4" fill="' + color + '" fill-opacity="0.8" ' +
        'stroke="var(--bar-black)" stroke-width="1" data-tip="' + tip + '" data-file="' + esc(r.file) + '" style="cursor:pointer;"/>';
    }).join('');

    var colorLabel = (state.colorKey && state.colorKey !== 'none') ? FIELDS[state.colorKey] : null;
    var title = esc(FIELDS[xKey]) + ' vs ' + esc(FIELDS[yKey]) + (colorLabel ? ', colored by ' + esc(colorLabel) : '');

    document.getElementById('scatter').innerHTML =
      '<div style="font-size:13px; font-weight:700; margin-bottom:6px;">' + title + '</div>' +
      '<svg width="' + width + '" height="' + height + '" viewBox="0 0 ' + width + ' ' + height + '" style="width:100%; height:auto;">' +
      grid + dots +
      '<text x="' + (padL + plotW / 2) + '" y="' + (height - 2) + '" font-size="10" text-anchor="middle" fill="var(--text-secondary)">' + esc(FIELDS[xKey]) + '</text>' +
      '<text x="12" y="' + (padT + plotH / 2) + '" font-size="10" text-anchor="middle" fill="var(--text-secondary)" transform="rotate(-90 12 ' + (padT + plotH / 2) + ')">' + esc(FIELDS[yKey]) + '</text>' +
      '</svg>';
    document.getElementById('scatter-legend').innerHTML = renderLegend(rows);
  }

  // ---- cell-type bar (also filter-reactive) ----
  function renderCellTypeBar() {
    var rows = applyFilters();
    var counts = {};
    rows.forEach(function(r) { counts[r.cell_type] = (counts[r.cell_type] || 0) + 1; });
    var items = Object.keys(counts).map(function(k) { return [k, counts[k]]; }).sort(function(a, b) { return b[1] - a[1]; });
    if (!items.length) { document.getElementById('celltype-bar').innerHTML = '<div class="subtle">no hits match current filters</div>'; return; }

    var width = 420, rowH = 26, gap = 6, padL = 130, padR = 40, padT = 6;
    var height = padT + items.length * (rowH + gap);
    var maxN = Math.max.apply(null, items.map(function(i) { return i[1]; }));
    var plotW = width - padL - padR;

    var parts = items.map(function(item, i) {
      var ct = item[0], n = item[1];
      var y = padT + i * (rowH + gap);
      var w = Math.max(n / maxN * plotW, 2);
      return '<text x="' + (padL - 8) + '" y="' + (y + rowH / 2 + 4).toFixed(1) + '" font-size="11" text-anchor="end">' + esc(ct) + '</text>' +
        '<rect x="' + padL + '" y="' + y + '" width="' + w.toFixed(1) + '" height="' + rowH + '" fill="var(--bar-black)" data-tip="' + esc(ct) + '&#10;' + n + ' hits"/>' +
        '<text x="' + (padL + w + 6).toFixed(1) + '" y="' + (y + rowH / 2 + 4).toFixed(1) + '" font-size="11" fill="var(--text-primary)">' + n + '</text>';
    }).join('');
    document.getElementById('celltype-bar').innerHTML =
      '<svg width="' + width + '" height="' + height + '" viewBox="0 0 ' + width + ' ' + height + '" style="width:100%; height:auto;">' + parts + '</svg>';
  }

  function renderAll() { renderTable(); renderPlot(); renderCellTypeBar(); }

  document.getElementById('search').addEventListener('input', function(e) {
    state.q = e.target.value.trim().toLowerCase(); renderAll();
  });
  ['cellType', 'group', 'bucket', 'evidence'].forEach(function(key) {
    document.getElementById('f-' + key).addEventListener('change', function(e) {
      state[key] = e.target.value; renderAll();
    });
  });
  document.querySelectorAll('th[data-key]').forEach(function(th) {
    th.addEventListener('click', function() {
      var k = th.getAttribute('data-key');
      if (state.sortKey === k) { state.sortDir *= -1; } else { state.sortKey = k; state.sortDir = 1; }
      document.querySelectorAll('th[data-key]').forEach(function(t) { t.classList.remove('sorted'); });
      th.classList.add('sorted');
      renderTable();
    });
  });
  // ---- plot config panels (scrollable option lists, not <select>) ----
  var XY_GROUPS = window.__XY_GROUPS__;         // numeric fields only
  var COLOR_GROUPS = window.__COLOR_GROUPS__;   // categorical + numeric, grouped by subject
  var SCALE_OPTIONS = [
    { key: 'sequential', label: 'sequential (blue)' },
    { key: 'diverging', label: 'diverging (centered at 0)' },
    { key: 'viridis', label: 'viridis' },
    { key: 'plasma', label: 'plasma' },
    { key: 'magma', label: 'magma' },
    { key: 'mako', label: 'mako' },
    { key: 'rocket', label: 'rocket' },
    { key: 'turbo', label: 'turbo' }
  ];

  function wireOptionRows(el, onPick) {
    el.querySelectorAll('.option-row').forEach(function(row) {
      row.addEventListener('click', function() { onPick(row.getAttribute('data-value')); });
    });
  }

  function renderFieldPanel(listId, groups, selectedKey, includeNone, onPick) {
    var el = document.getElementById(listId);
    var html = includeNone
      ? '<div class="option-row' + (selectedKey === 'none' ? ' selected' : '') + '" data-value="none">none (monochrome)</div>'
      : '';
    groups.forEach(function(group) {
      html += '<div class="option-group-label">' + esc(group[0]) + '</div>';
      group[1].forEach(function(item) {
        html += '<div class="option-row' + (item[0] === selectedKey ? ' selected' : '') + '" data-value="' + esc(item[0]) + '">' +
          esc(item[1]) + '</div>';
      });
    });
    el.innerHTML = html;
    wireOptionRows(el, onPick);
  }

  function scaleGradientCss(key) {
    var stops = key === 'diverging'
      ? [scaleColor('diverging', 0, -1), scaleColor('diverging', 0.5, 0), scaleColor('diverging', 1, 1)]
      : [0, 0.25, 0.5, 0.75, 1].map(function(t) { return scaleColor(key, t, t); });
    return 'linear-gradient(to right,' + stops.join(',') + ')';
  }

  function renderScalePanel(listId, selectedKey, onPick) {
    var el = document.getElementById(listId);
    el.innerHTML = SCALE_OPTIONS.map(function(opt) {
      return '<div class="option-row' + (opt.key === selectedKey ? ' selected' : '') + '" data-value="' + opt.key + '">' +
        '<span class="swatch-preview" style="background:' + scaleGradientCss(opt.key) + ';"></span>' + esc(opt.label) + '</div>';
    }).join('');
    wireOptionRows(el, onPick);
  }

  function refreshPanels() {
    document.getElementById('p-colorScheme-wrap').style.display = isNumeric(state.colorKey) ? '' : 'none';
    renderFieldPanel('p-xKey-list', XY_GROUPS, state.xKey, false, function(v) { state.xKey = v; refreshPanels(); renderPlot(); });
    renderFieldPanel('p-yKey-list', XY_GROUPS, state.yKey, false, function(v) { state.yKey = v; refreshPanels(); renderPlot(); });
    renderFieldPanel('p-colorKey-list', COLOR_GROUPS, state.colorKey, true, function(v) { state.colorKey = v; refreshPanels(); renderPlot(); });
    renderScalePanel('p-colorScheme-list', state.colorScheme, function(v) { state.colorScheme = v; refreshPanels(); renderPlot(); });
  }

  refreshPanels();
  renderAll();
})();
"""


def _field_groups(fields: list[tuple[str, str, str]]) -> list[list]:
    """[key, label, group] triples -> [[group, [[key, label], ...]], ...],
    order preserved. JS builds the scrollable option panels from this
    instead of native <optgroup>/<option> HTML."""
    groups: dict[str, list[list[str]]] = {}
    for k, v, g in fields:
        groups.setdefault(g, []).append([k, v])
    return [[g, items] for g, items in groups.items()]


def render_index(records: list[dict]) -> str:
    n = len(records)
    n_tf = sum(1 for r in records if r["master_group"] == "trial_failure_candidate")
    n_dr = n - n_tf
    n_ev = sum(1 for r in records if r["has_evidence"])

    cell_types = sorted({r["cell_type"] for r in records})
    cell_type_options = "".join(f'<option value="{c}">{c}</option>' for c in cell_types)

    stat_tiles = f"""
<div class="stat-row">
  <div class="stat-tile"><div class="value">{n}</div><div class="label">total hits</div></div>
  <div class="stat-tile"><div class="value" style="color:var(--tf);">{n_tf}</div><div class="label">trial failure</div></div>
  <div class="stat-tile"><div class="value" style="color:var(--dr);">{n_dr}</div><div class="label">drug repurposing</div></div>
  <div class="stat-tile"><div class="value">{n_ev}</div><div class="label">have drug evidence</div></div>
</div>"""

    charts = f"""
<div class="card" style="margin-top:32px;">
  <h2>Explore</h2>
  <div class="control-panel">
    <div class="control-group">
      <div class="control-label">x-axis</div>
      <div class="option-list" id="p-xKey-list"></div>
    </div>
    <div class="control-group">
      <div class="control-label">y-axis</div>
      <div class="option-list" id="p-yKey-list"></div>
    </div>
    <div class="control-group">
      <div class="control-label">color by</div>
      <div class="option-list" id="p-colorKey-list"></div>
    </div>
    <div class="control-group" id="p-colorScheme-wrap" style="display:none;">
      <div class="control-label">color scale</div>
      <div class="option-list" id="p-colorScheme-list"></div>
    </div>
  </div>
  <div class="charts-row">
    <div>
      <div id="scatter"></div>
      <div id="scatter-legend"></div>
    </div>
    <div>
      <div class="subtle" style="margin-bottom:8px; font-weight:700;">hits per cell type</div>
      <div id="celltype-bar"></div>
    </div>
  </div>
</div>"""

    table = f"""
<div class="card">
  <h2>All candidates</h2>
  <div class="controls">
    <input type="text" id="search" placeholder="search gene name...">
    <select id="f-cellType"><option value="">all cell types</option>{cell_type_options}</select>
    <select id="f-group">
      <option value="">all groups</option>
      <option value="trial_failure_candidate">trial failure</option>
      <option value="drug_repurposing_candidate">drug repurposing</option>
    </select>
    <select id="f-bucket">
      <option value="">any length change</option>
      <option value="shorter">shorter</option>
      <option value="longer">longer</option>
      <option value="same">same length</option>
    </select>
    <select id="f-evidence">
      <option value="">any drug evidence</option>
      <option value="yes">has evidence</option>
      <option value="no">no evidence</option>
    </select>
    <span class="count" id="count"></span>
  </div>
  <div class="table-wrap">
    <table>
      <thead>
        <tr>
          <th data-key="gene">Gene</th>
          <th data-key="cell_type">Cell type</th>
          <th data-key="master_group">Group</th>
          <th data-key="protein_change_type">Protein change</th>
          <th data-key="padj" class="sorted">padj</th>
          <th data-key="control_pct">Control</th>
          <th data-key="ad_pct">AD</th>
          <th data-key="delta_pp">&Delta; pp</th>
          <th data-key="n_alts"># alts</th>
          <th>Drug evidence</th>
        </tr>
      </thead>
      <tbody id="tbody"></tbody>
    </table>
  </div>
</div>"""

    field_labels = {k: v for k, v, _ in NUMERIC_FIELDS + CATEGORICAL_FIELDS}
    numeric_keys = [k for k, _, _ in NUMERIC_FIELDS]
    xy_groups = _field_groups(NUMERIC_FIELDS)
    color_groups = _field_groups(CATEGORICAL_FIELDS + NUMERIC_FIELDS)

    return f"""<!doctype html>
<html><head><meta charset="utf-8"/>
<title>Master Surveyor &middot; Candidate Dossiers</title>
<style>{_font_face_css()}{_CSS}</style>
</head><body>
{THEME_TOGGLE_HTML}
<div class="viz-root">
  <div class="card">
    <h1>Master Surveyor</h1>
    <div class="subtle">{n} candidates that passed J4's gate into master_surveyor &middot; click any row or chart point to open its dossier</div>
  </div>
  {stat_tiles}
  {charts}
  {table}
</div>
<div class="tip" id="tip"></div>
<script>{THEME_TOGGLE_JS}</script>
<script>
  window.__HITS__ = {json.dumps(records)};
  window.__FIELDS__ = {json.dumps(field_labels)};
  window.__NUMERIC_KEYS__ = {json.dumps(numeric_keys)};
  window.__CAT_PALETTE__ = {json.dumps(palette.isoform_colors(8))};
  window.__XY_GROUPS__ = {json.dumps(xy_groups)};
  window.__COLOR_GROUPS__ = {json.dumps(color_groups)};
</script>
<script>{_JS}</script>
</body></html>"""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    records = build_manifest()
    out_path = OUT_DIR / "index.html"
    out_path.write_text(render_index(records))
    print(f"wrote {out_path} ({len(records)} hits)")


if __name__ == "__main__":
    main()
