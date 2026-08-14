"""Interactive index page -- live cutoff filters, per-hit cart selection
(with per-drug checkboxes), and an export-to-Discovery-Studio flow with job
progress polling. Vanilla HTML/CSS/JS, no build step, same philosophy as
dossier/generate_index.py's static index (which this supersedes as the
*interactive* homepage; the static one still exists for offline browsing).
"""

from __future__ import annotations

_CSS = """
@font-face {
  font-family: "Tamzen";
  src: url("/static/fonts/Tamzen8x16r.ttf") format("truetype");
  font-weight: 400; font-style: normal; font-display: swap;
}
@font-face {
  font-family: "Tamzen";
  src: url("/static/fonts/Tamzen8x16b.ttf") format("truetype");
  font-weight: 700; font-style: normal; font-display: swap;
}
:root { color-scheme: light dark; }
body {
  font-family: "Tamzen", "SF Mono", ui-monospace, monospace; margin: 0; padding: 24px 20px 64px;
  background: #ffffff; color: #1a1a1a; max-width: 1500px; margin-inline: auto;
}
@media (prefers-color-scheme: dark) { body { background: #0b0b0f; color: #e6e6e6; } }
.brand { display: flex; align-items: center; gap: 12px; }
.brand img { width: 40px; height: 40px; flex-shrink: 0; }
h1 { font-size: 22px; margin: 0; }
.subtle { color: #767671; font-size: 12px; }
.layout { display: grid; grid-template-columns: 1fr 320px; gap: 20px; align-items: start; margin-top: 16px; }
.layout.graph-tab-active { grid-template-columns: 1fr 460px; }
@media (max-width: 1000px) { .layout, .layout.graph-tab-active { grid-template-columns: 1fr; } }
.card { border: 1px solid #d8d8d2; border-radius: 4px; padding: 14px 16px; margin-bottom: 16px; }
@media (prefers-color-scheme: dark) { .card { border-color: #2a2a30; } }
.filters { display: flex; flex-wrap: wrap; gap: 14px; align-items: flex-end; }
.filter-group { display: flex; flex-direction: column; gap: 4px; font-size: 11px; }
.filter-group input[type=text] { font-family: inherit; font-size: 12px; padding: 5px 8px; min-width: 200px; }
.filter-group input[type=number] { font-family: inherit; font-size: 12px; padding: 5px 8px; width: 90px; }
.filter-group label { text-transform: uppercase; letter-spacing: 0.04em; color: #767671; font-size: 10px; }
.toggle-group { display: flex; flex-wrap: wrap; gap: 6px; }
.toggle-btn {
  font-family: inherit; font-size: 11px; padding: 4px 9px; border: 1px solid #d8d8d2;
  background: none; border-radius: 3px; cursor: pointer; color: #767671;
}
@media (prefers-color-scheme: dark) { .toggle-btn { border-color: #2a2a30; } }
.toggle-btn.active { background: #1a1a1a; color: #ffffff; border-color: #1a1a1a; }
@media (prefers-color-scheme: dark) { .toggle-btn.active { background: #e6e6e6; color: #0b0b0f; border-color: #e6e6e6; } }
table { border-collapse: collapse; width: 100%; font-size: 12px; }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid #eee; white-space: nowrap; vertical-align: top; line-height: 1.4; }
@media (prefers-color-scheme: dark) { th, td { border-color: #222; } }
th { position: sticky; top: 0; background: inherit; cursor: pointer; vertical-align: middle; }
.table-wrap { max-height: 640px; overflow: auto; border: 1px solid #d8d8d2; }
@media (prefers-color-scheme: dark) { .table-wrap { border-color: #2a2a30; } }
.badge { display: inline-block; padding: 1px 6px; border: 1px solid currentColor; border-radius: 3px; font-size: 10px; }
.tf { color: #4a3aa7; } .dr { color: #1baf7a; }
button { font-family: inherit; cursor: pointer; }
select { font-family: inherit; font-size: 12px; }
.drug-list { display: flex; flex-direction: column; gap: 2px; }
.drug-list label {
  font-size: 11px; white-space: nowrap; line-height: 1.4;
  display: flex; align-items: center; gap: 4px; margin: 0;
}
.drug-list input[type=checkbox] { margin: 0; }
.cart-item { border-bottom: 1px solid #eee; padding: 6px 0; font-size: 12px; }
@media (prefers-color-scheme: dark) { .cart-item { border-color: #222; } }
.cart-item button { font-size: 11px; }
#export-status { font-size: 11px; }
.job-item { margin-bottom: 6px; }
.count { font-size: 11px; color: #767671; }

#splash {
  position: fixed; inset: 0; z-index: 100; display: flex; flex-direction: column;
  align-items: center; justify-content: center; gap: 14px;
  background: #ffffff; transition: opacity 0.25s ease;
}
@media (prefers-color-scheme: dark) { #splash { background: #0b0b0f; } }
#splash.hidden { opacity: 0; pointer-events: none; }
#splash img { width: 96px; height: 96px; }
#splash .spinner {
  width: 28px; height: 28px; border-radius: 50%;
  border: 3px solid #d8d8d2; border-top-color: #1a1a1a;
  animation: splash-spin 0.8s linear infinite;
}
@media (prefers-color-scheme: dark) { #splash .spinner { border-color: #2a2a30; border-top-color: #e6e6e6; } }
#splash .label { font-size: 12px; color: #767671; letter-spacing: 0.03em; }
@keyframes splash-spin { to { transform: rotate(360deg); } }

.plot-toolbar { display: flex; gap: 14px; align-items: center; margin-bottom: 10px; flex-wrap: wrap; }
.plot-toolbar .filter-group { min-width: 160px; }
.plot-toolbar select { font-family: inherit; font-size: 12px; padding: 5px 8px; width: 100%; }
#scatter-wrap { position: relative; }
#scatter { user-select: none; }
#scatter svg { display: block; width: 100%; height: auto; touch-action: none; }
#scatter circle.dot { cursor: pointer; transition: stroke-width 0.1s ease, r 0.1s ease; }
#scatter circle.dot.selected { stroke: #1baf7a; stroke-width: 2.5px; }
#scatter circle.dot.hovered { stroke: #e8a33d; stroke-width: 3.5px; r: 6; }
#scatter rect.select-rect { fill: rgba(27,175,122,0.15); stroke: #1baf7a; stroke-width: 1; stroke-dasharray: 3,2; pointer-events: none; }
.plot-actions { display: flex; gap: 8px; align-items: center; margin-top: 10px; flex-wrap: wrap; }
.plot-actions .count { font-size: 12px; }
.plot-actions button { font-size: 11px; }

.tab-bar { display: flex; gap: 4px; margin-bottom: 12px; border-bottom: 1px solid #d8d8d2; }
@media (prefers-color-scheme: dark) { .tab-bar { border-color: #2a2a30; } }
.tab-btn {
  font-family: inherit; font-size: 12px; padding: 8px 14px; background: none; border: none;
  border-bottom: 2px solid transparent; color: #767671; cursor: pointer; margin-bottom: -1px;
}
.tab-btn.active { color: #1a1a1a; border-bottom-color: #1baf7a; font-weight: 700; }
@media (prefers-color-scheme: dark) { .tab-btn.active { color: #e6e6e6; } }
.tab-content[hidden] { display: none; }

.legend-row { display: flex; flex-wrap: wrap; gap: 8px; font-size: 11px; color: #767671; margin-top: 4px; }
.legend-row .swatch { display: inline-block; width: 9px; height: 9px; border-radius: 2px; margin-right: 4px; vertical-align: -1px; }
.legend-gradient { width: 100%; height: 10px; border: 1px solid #d8d8d2; margin-top: 6px; }
@media (prefers-color-scheme: dark) { .legend-gradient { border-color: #2a2a30; } }

.graph-tip {
  position: fixed; pointer-events: none; z-index: 90; display: none;
  background: #ffffff; color: #0b0b0b; border: 1px solid rgba(11,11,11,0.15);
  box-shadow: 0 2px 10px rgba(11,11,11,0.18); font-size: 11px; padding: 6px 9px;
  border-radius: 3px; line-height: 1.5; white-space: pre; max-width: 260px;
}
.selection-item {
  display: flex; align-items: center; justify-content: space-between; gap: 8px;
  border-bottom: 1px solid #eee; padding: 5px 0; font-size: 12px;
}
@media (prefers-color-scheme: dark) { .selection-item { border-color: #222; } }
.selection-item button { font-size: 11px; }
"""

_JS = r"""
(function() {
  var currentHits = [];
  var thresholds = {tf_min_abs_delta_usage: 0, dr_min_chembl_or_ot_phase: 0};
  var SPLASH_MIN_MS = 1500;   // floor so the splash is actually perceptible,
                               // not just a flash on a fast local fetch
  var splashShownAt = Date.now();

  // Plottable numeric fields -- same set dossier/generate_index.py's static
  // scatter plot exposes (kept in sync manually; both read from the same
  // hits_deep.csv-derived record shape via dossier/manifest.py).
  var NUMERIC_FIELDS = [
    ['delta_pp', 'Δ usage Control→AD (pp)'],
    ['neg_log10_padj', '-log10(padj)'],
    ['padj', 'padj (raw)'],
    ['control_pct', 'Control usage (%)'],
    ['ad_pct', 'AD usage (%)'],
    ['n_alts', '# ranked alts'],
    ['protein_length_diff', 'Protein length change (aa)'],
    ['pct_identity', 'Protein identity to canonical'],
    ['changed_aa_fraction', 'Fraction of protein changed'],
    ['n_domains', 'Pfam domains (canonical)'],
    ['chembl_max_phase', 'ChEMBL max phase'],
    ['n_drugs', 'ChEMBL drugs (n)'],
    ['chembl_bioactive_compounds', 'ChEMBL bioactive compounds (n)'],
    ['chembl_best_pchembl', 'ChEMBL best pChEMBL'],
    ['ot_max_phase', 'Open Targets max phase'],
    ['ot_n_drugs', 'Open Targets drugs (n)'],
    ['ot_trials_total', 'Clinical trials on record (n)'],
    ['ot_trials_terminated', 'Trials stopped early (n)'],
    ['pharos_n_ligands', 'Pharos ligands (n)'],
    ['pharos_n_drugs', 'Pharos drugs (n)'],
    ['dgidb_interactions', 'DGIdb interactions (n)'],
  ];
  // Categorical fields available for coloring -- same set generate_index.py exposes.
  var CATEGORICAL_FIELDS = [
    ['master_group', 'Group'],
    ['cell_type', 'Cell type'],
    ['protein_change_type', 'Protein change type'],
    ['change_bucket', 'Length change'],
    ['has_evidence', 'Has drug evidence'],
    ['has_failed_trial', 'Has a stopped clinical trial'],
    ['pharos_tdl', 'Pharos target development level'],
  ];
  var FIELD_LABEL = {};
  NUMERIC_FIELDS.forEach(function(f) { FIELD_LABEL[f[0]] = f[1]; });
  CATEGORICAL_FIELDS.forEach(function(f) { FIELD_LABEL[f[0]] = f[1]; });

  // ---- color helpers (ported from dossier/generate_index.py's scatter plot) ----
  var CAT_PALETTE = ['#527fb4', '#4a9f7e', '#c39737', '#3a8637', '#655b98', '#be6461', '#c18197', '#c27656'];
  var OTHER_GRAY = '#898781';
  var FIXED_COLORS = {
    master_group: { trial_failure_candidate: '#4a3aa7', drug_repurposing_candidate: '#1baf7a' },
    change_bucket: { shorter: '#1baf7a', longer: '#4a3aa7', same: '#000000' },
  };
  var COLORMAPS = {
    viridis: [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[110,206,88],[253,231,37]],
    plasma:  [[13,8,135],[71,3,159],[115,1,168],[156,23,158],[189,55,134],[216,87,107],[237,121,83],[250,159,58],[240,249,33]],
    magma:   [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],[229,80,100],[251,135,97],[254,194,135],[252,253,191]],
    mako:    [[11,5,5],[35,22,58],[35,58,89],[26,89,107],[24,116,110],[41,145,101],[93,171,89],[168,196,90],[222,222,142]],
    rocket:  [[3,5,26],[43,20,74],[95,27,109],[144,33,105],[192,47,84],[224,79,63],[241,127,67],[247,179,110],[250,234,190]],
    turbo:   [[48,18,59],[70,107,227],[39,168,224],[55,207,157],[151,222,73],[233,196,49],[247,131,44],[211,60,32],[122,4,3]],
  };
  var COLOR_SCALES = [
    ['sequential', 'sequential (blue)'], ['diverging', 'diverging (centered at 0)'],
    ['viridis', 'viridis'], ['plasma', 'plasma'], ['magma', 'magma'],
    ['mako', 'mako'], ['rocket', 'rocket'], ['turbo', 'turbo'],
  ];

  function lerp(a, b, t) { return a + (b - a) * t; }
  function rgbHex(rgb) {
    return '#' + rgb.map(function(v) {
      var h = Math.max(0, Math.min(255, Math.round(v))).toString(16);
      return h.length === 1 ? '0' + h : h;
    }).join('');
  }
  function seqColor(t) { return rgbHex([lerp(205, 13, t), lerp(226, 54, t), lerp(251, 107, t)]); }
  function divColor(t) {
    var target = t < 0 ? [42, 120, 214] : [227, 73, 72];
    var a = Math.min(Math.abs(t), 1);
    return rgbHex([lerp(240, target[0], a), lerp(239, target[1], a), lerp(236, target[2], a)]);
  }
  function colormapColor(name, t) {
    var stops = COLORMAPS[name];
    t = Math.max(0, Math.min(1, t));
    var pos = t * (stops.length - 1);
    var i = Math.min(Math.floor(pos), stops.length - 2);
    var frac = pos - i;
    return rgbHex([
      lerp(stops[i][0], stops[i + 1][0], frac),
      lerp(stops[i][1], stops[i + 1][1], frac),
      lerp(stops[i][2], stops[i + 1][2], frac),
    ]);
  }
  function scaleColor(scheme, t01, tSigned) {
    if (scheme === 'diverging') return divColor(tSigned);
    if (COLORMAPS[scheme]) return colormapColor(scheme, t01);
    return seqColor(t01);
  }
  function categoricalColorMap(rows, key) {
    var counts = {};
    rows.forEach(function(r) { var v = String(r[key]); counts[v] = (counts[v] || 0) + 1; });
    var sorted = Object.keys(counts).sort(function(a, b) { return counts[b] - counts[a]; });
    var realSlots = sorted.length <= CAT_PALETTE.length ? CAT_PALETTE.length : CAT_PALETTE.length - 1;
    var map = {};
    sorted.forEach(function(v, i) { map[v] = i < realSlots ? CAT_PALETTE[i] : OTHER_GRAY; });
    return map;
  }
  function isNumericField(key) { return !!FIELD_LABEL[key] && NUMERIC_FIELDS.some(function(f) { return f[0] === key; }); }
  function colorFor(r, rows, key) {
    if (!key || key === 'none') return '#000000';
    if (key === 'has_evidence') return r.has_evidence ? '#0ca30c' : '#d03b3b';
    if (key === 'has_failed_trial') return r.has_failed_trial ? '#d03b3b' : '#cccccc';
    if (FIXED_COLORS[key]) return FIXED_COLORS[key][r[key]] || OTHER_GRAY;
    if (isNumericField(key)) {
      var vals = rows.map(function(x) { return x[key]; }).filter(function(v) { return v !== null && v !== undefined; });
      if (!vals.length) return OTHER_GRAY;
      var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
      var v = r[key];
      if (v === null || v === undefined) return '#cccccc';
      var extent = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
      return scaleColor(plotState.colorScale, hi > lo ? (v - lo) / (hi - lo) : 0.5, v / extent);
    }
    return categoricalColorMap(rows, key)[String(r[key])] || OTHER_GRAY;
  }
  var plotState = { xKey: 'delta_pp', yKey: 'neg_log10_padj', colorKey: 'master_group', colorScale: 'sequential' };
  // Display-only filter (no refetch) -- all three toggled on by default,
  // though novel_target_candidate never actually appears in /api/hits
  // (out of master_surveyor's scope), so that toggle is a no-op in practice.
  var groupFilter = { trial_failure_candidate: true, drug_repurposing_candidate: true, novel_target_candidate: true };
  function groupVisible(h) { return !!groupFilter[h.master_group]; }
  var selectedIdx = {};       // idx (into currentHits) -> true, persists across re-renders until cleared
  var selectionHistory = [];  // stack of prior selectedIdx snapshots, one per completed drag -- "undo"
  var plotPoints = [];        // [{idx, x, y}] in SVG pixel space, from the last renderPlot() call

  function qs(obj) {
    return Object.keys(obj).map(function(k) { return encodeURIComponent(k) + '=' + encodeURIComponent(obj[k]); }).join('&');
  }

  function hideSplash() {
    var splash = document.getElementById('splash');
    if (!splash) return;
    var elapsed = Date.now() - splashShownAt;
    var remaining = SPLASH_MIN_MS - elapsed;
    if (remaining > 0) {
      setTimeout(function() { splash.classList.add('hidden'); }, remaining);
    } else {
      splash.classList.add('hidden');
    }
  }

  function fetchHits() {
    var params = Object.assign({}, thresholds);
    var q = document.getElementById('search').value.trim();
    if (q) params.q = q;
    fetch('/api/hits?' + qs(params)).then(function(r) { return r.json(); }).then(function(data) {
      currentHits = data.hits;
      document.getElementById('count').textContent = data.total + ' hits';
      selectedIdx = {};        // a fresh fetch (filter/search change) invalidates old point indices
      selectionHistory = [];
      renderTable();
      renderPlot();
      renderGraphSelectionList();
      hideSplash();
    }).catch(hideSplash);
  }

  function renderTable() {
    var tbody = document.getElementById('hits-body');
    tbody.innerHTML = '';
    currentHits.forEach(function(h, idx) {
      if (!groupVisible(h)) return;
      var tr = document.createElement('tr');
      var groupClass = h.master_group === 'trial_failure_candidate' ? 'tf' : 'dr';
      var groupLabel = h.master_group === 'trial_failure_candidate' ? 'TF' : 'DR';
      var drugs = (h.drug_names || []).concat(h.ot_drug_names || []).filter(function(v, i, a) { return a.indexOf(v) === i; });
      var drugInputs = drugs.map(function(d, i) {
        return '<label><input type="checkbox" data-drug="' + idx + '" value="' + d + '"> ' + d + '</label>';
      }).join('');
      var detailUrl = '/hit/' + encodeURIComponent(h.gene) + '/' + encodeURIComponent(h.cell_type) + '?hit_enst=' + encodeURIComponent(h.hit_enst);
      var isTf = h.master_group === 'trial_failure_candidate';
      function foldSpan(ok, text) {
        return '<span style="color:' + (ok ? '#0ca30c' : '#d03b3b') + ';">' + text + '</span>';
      }
      function boolCell(v) {
        return foldSpan(!!v, v ? 'Yes' : 'No');
      }
      function fractionCell(count) {
        if (!h.alt_total) return '<span>-</span>';
        return foldSpan(count === h.alt_total, count + '/' + h.alt_total);
      }
      function foldCell(count) {
        return isTf ? fractionCell(count) : boolCell(count);
      }
      tr.innerHTML =
        '<td><span class="badge ' + groupClass + '">' + groupLabel + '</span></td>' +
        '<td><a href="' + detailUrl + '">' + h.gene + '</a></td>' +
        '<td>' + h.cell_type + '</td>' +
        '<td>' + h.protein_change_type + '</td>' +
        '<td>' + h.delta_pp + '</td>' +
        '<td>' + boolCell(h.canonical_colabfold) + '</td>' +
        '<td>' + boolCell(h.canonical_esmfold) + '</td>' +
        '<td>' + foldCell(h.alt_colabfold_folded) + '</td>' +
        '<td>' + foldCell(h.alt_esmfold_folded) + '</td>' +
        '<td class="drug-list">' + (drugInputs || '<span class="subtle">none</span>') + '</td>' +
        '<td><button data-add="' + idx + '">+ cart</button></td>';
      tbody.appendChild(tr);
    });
    tbody.querySelectorAll('button[data-add]').forEach(function(btn) {
      btn.addEventListener('click', function() {
        var idx = parseInt(btn.getAttribute('data-add'), 10);
        addToCart(idx);
      });
    });
  }

  function addToCart(idx) {
    var h = currentHits[idx];
    var checked = document.querySelectorAll('input[data-drug="' + idx + '"]:checked');
    var selected = Array.prototype.map.call(checked, function(c) { return c.value; });
    fetch('/api/cart/items', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({gene: h.gene, cell_type: h.cell_type, hit_enst: h.hit_enst, selected_drugs: selected}),
    }).then(function(r) { return r.json(); }).then(function(data) { renderCart(data.items); });
  }

  // ---- scatter plot: axis-selectable, drag-to-box-select multiple points ----
  function initPlotSelectors() {
    ['plot-x', 'plot-y'].forEach(function(id, i) {
      var sel = document.getElementById(id);
      sel.innerHTML = NUMERIC_FIELDS.map(function(f) {
        return '<option value="' + f[0] + '">' + f[1] + '</option>';
      }).join('');
      sel.value = i === 0 ? plotState.xKey : plotState.yKey;
      sel.addEventListener('change', function() {
        plotState[i === 0 ? 'xKey' : 'yKey'] = sel.value;
        renderPlot();
      });
    });

    var colorSel = document.getElementById('plot-color-key');
    colorSel.innerHTML = '<option value="none">none (monochrome)</option>' +
      '<optgroup label="Categorical">' +
      CATEGORICAL_FIELDS.map(function(f) { return '<option value="' + f[0] + '">' + f[1] + '</option>'; }).join('') +
      '</optgroup><optgroup label="Numeric">' +
      NUMERIC_FIELDS.map(function(f) { return '<option value="' + f[0] + '">' + f[1] + '</option>'; }).join('') +
      '</optgroup>';
    colorSel.value = plotState.colorKey;
    colorSel.addEventListener('change', function() { plotState.colorKey = colorSel.value; renderPlot(); });

    var scaleSel = document.getElementById('plot-color-scale');
    scaleSel.innerHTML = COLOR_SCALES.map(function(s) { return '<option value="' + s[0] + '">' + s[1] + '</option>'; }).join('');
    scaleSel.value = plotState.colorScale;
    scaleSel.addEventListener('change', function() { plotState.colorScale = scaleSel.value; renderPlot(); });
  }

  function renderLegend(rows) {
    var key = plotState.colorKey;
    var el = document.getElementById('plot-legend');
    if (!key || key === 'none') { el.innerHTML = ''; return; }

    if (key === 'has_evidence') {
      el.innerHTML = '<div class="legend-row"><span><span class="swatch" style="background:#0ca30c;"></span>has evidence</span>' +
        '<span><span class="swatch" style="background:#d03b3b;"></span>no evidence</span></div>';
      return;
    }
    if (key === 'has_failed_trial') {
      el.innerHTML = '<div class="legend-row"><span><span class="swatch" style="background:#d03b3b;"></span>trial stopped early</span>' +
        '<span><span class="swatch" style="background:#cccccc;"></span>no stopped trial</span></div>';
      return;
    }
    if (FIXED_COLORS[key]) {
      el.innerHTML = '<div class="legend-row">' + Object.keys(FIXED_COLORS[key]).map(function(k) {
        return '<span><span class="swatch" style="background:' + FIXED_COLORS[key][k] + ';"></span>' + k + '</span>';
      }).join('') + '</div>';
      return;
    }
    if (isNumericField(key)) {
      var stops = plotState.colorScale === 'diverging'
        ? [scaleColor('diverging', 0, -1), scaleColor('diverging', 0.5, 0), scaleColor('diverging', 1, 1)]
        : [0, 0.25, 0.5, 0.75, 1].map(function(t) { return scaleColor(plotState.colorScale, t, t); });
      el.innerHTML = '<div class="legend-gradient" style="background:linear-gradient(to right,' + stops.join(',') + ');"></div>' +
        '<div class="subtle" style="margin-top:3px;">' + FIELD_LABEL[key] + ' (low &rarr; high)</div>';
      return;
    }
    var map = categoricalColorMap(rows, key);
    el.innerHTML = '<div class="legend-row">' + Object.keys(map).map(function(k) {
      return '<span><span class="swatch" style="background:' + map[k] + ';"></span>' + k + '</span>';
    }).join('') + '</div>';
  }

  function hitUrl(h) {
    return '/hit/' + encodeURIComponent(h.gene) + '/' + encodeURIComponent(h.cell_type) + '?hit_enst=' + encodeURIComponent(h.hit_enst);
  }

  function renderGraphSelectionList() {
    var indices = Object.keys(selectedIdx);
    document.getElementById('graph-sel-count').textContent = indices.length ? (indices.length + ' selected') : '';
    var el = document.getElementById('graph-selection-items');
    if (!indices.length) {
      el.innerHTML = '<div class="subtle">drag on the graph to select points</div>';
      return;
    }
    el.innerHTML = indices.map(function(idxStr) {
      var idx = parseInt(idxStr, 10);
      var h = currentHits[idx];
      if (!h) return '';
      return '<div class="selection-item" data-idx="' + idx + '">' +
        '<span><b>' + h.gene + '</b> / ' + h.cell_type + '</span>' +
        '<span><a href="' + hitUrl(h) + '">open</a> ' +
        '<button data-unselect="' + idx + '" title="remove from selection">&times;</button></span></div>';
    }).join('');
    el.querySelectorAll('button[data-unselect]').forEach(function(btn) {
      btn.addEventListener('click', function() {
        var idx = parseInt(btn.getAttribute('data-unselect'), 10);
        delete selectedIdx[idx];
        var circle = document.querySelector('#scatter circle[data-idx="' + idx + '"]');
        if (circle) circle.classList.remove('selected');
        renderGraphSelectionList();
      });
    });
    el.querySelectorAll('.selection-item').forEach(function(row) {
      var idx = row.getAttribute('data-idx');
      row.addEventListener('pointerenter', function() {
        var circle = document.querySelector('#scatter circle[data-idx="' + idx + '"]');
        if (circle) circle.classList.add('hovered');
      });
      row.addEventListener('pointerleave', function() {
        var circle = document.querySelector('#scatter circle[data-idx="' + idx + '"]');
        if (circle) circle.classList.remove('hovered');
      });
    });
  }

  function undoSelection() {
    if (!selectionHistory.length) return;
    selectedIdx = selectionHistory.pop();
    document.querySelectorAll('#scatter circle.dot').forEach(function(c) {
      var idx = parseInt(c.getAttribute('data-idx'), 10);
      c.classList.toggle('selected', !!selectedIdx[idx]);
    });
    renderGraphSelectionList();
  }

  function renderPlot() {
    var xKey = plotState.xKey, yKey = plotState.yKey;
    var width = 700, height = 380, padL = 56, padR = 16, padT = 14, padB = 34;
    var plotW = width - padL - padR, plotH = height - padT - padB;
    var container = document.getElementById('scatter');
    plotPoints = [];

    var usable = currentHits
      .map(function(h, idx) { return {h: h, idx: idx}; })
      .filter(function(p) {
        return groupVisible(p.h) &&
          p.h[xKey] !== null && p.h[xKey] !== undefined && p.h[yKey] !== null && p.h[yKey] !== undefined;
      });

    if (!usable.length) {
      container.innerHTML = '<div class="subtle">no hits match current filters</div>';
      return;
    }

    var xs = usable.map(function(p) { return p.h[xKey]; });
    var ys = usable.map(function(p) { return p.h[yKey]; });
    var xMin = Math.min.apply(null, xs.concat([0])), xMax = Math.max.apply(null, xs.concat([0]));
    var yMin = Math.min.apply(null, ys.concat([0])), yMax = Math.max.apply(null, ys.concat([0]));
    var xPad = (xMax - xMin) * 0.08 || 1, yPad = (yMax - yMin) * 0.08 || 1;
    xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

    function xOf(v) { return padL + (v - xMin) / (xMax - xMin) * plotW; }
    function yOf(v) { return padT + plotH * (1 - (v - yMin) / (yMax - yMin)); }

    var grid = '';
    for (var i = 0; i <= 4; i++) {
      var frac = i / 4;
      var gy = padT + plotH * (1 - frac), gx = padL + plotW * frac;
      grid += '<line x1="' + padL + '" y1="' + gy.toFixed(1) + '" x2="' + (padL + plotW) + '" y2="' + gy.toFixed(1) + '" stroke="#e5e5e0" stroke-width="1"/>';
      grid += '<text x="' + (padL - 8) + '" y="' + (gy + 4).toFixed(1) + '" font-size="10" text-anchor="end" fill="#767671">' + Math.round(yMin + (yMax - yMin) * frac) + '</text>';
      grid += '<text x="' + gx.toFixed(1) + '" y="' + (height - 8) + '" font-size="10" text-anchor="middle" fill="#767671">' + Math.round(xMin + (xMax - xMin) * frac) + '</text>';
    }

    // x=0 / y=0 reference axes (solid black, distinct from the hairline
    // gridlines above) plus a red marker at the origin itself.
    var zx = xOf(0), zy = yOf(0);
    grid += '<line x1="' + zx.toFixed(1) + '" y1="' + padT + '" x2="' + zx.toFixed(1) + '" y2="' + (padT + plotH) + '" stroke="#000000" stroke-width="1.3"/>';
    grid += '<line x1="' + padL + '" y1="' + zy.toFixed(1) + '" x2="' + (padL + plotW) + '" y2="' + zy.toFixed(1) + '" stroke="#000000" stroke-width="1.3"/>';
    grid += '<circle cx="' + zx.toFixed(1) + '" cy="' + zy.toFixed(1) + '" r="3.5" fill="#d03b3b" stroke="#000000" stroke-width="0.75"/>';

    var rowsForColor = usable.map(function(p) { return p.h; });
    var dots = usable.map(function(p) {
      var x = xOf(p.h[xKey]), y = yOf(p.h[yKey]);
      plotPoints.push({idx: p.idx, x: x, y: y});
      var cls = 'dot' + (selectedIdx[p.idx] ? ' selected' : '');
      var color = colorFor(p.h, rowsForColor, plotState.colorKey);
      return '<circle class="' + cls + '" data-idx="' + p.idx + '" cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) + '" r="4" ' +
        'fill="' + color + '" fill-opacity="0.8" stroke="#000" stroke-width="1"/>';
    }).join('');

    container.innerHTML =
      '<svg id="scatter-svg" viewBox="0 0 ' + width + ' ' + height + '">' + grid + dots +
      '<text x="' + (padL + plotW / 2) + '" y="' + height + '" font-size="10" text-anchor="middle" fill="#767671">' + FIELD_LABEL[xKey] + '</text>' +
      '<text x="12" y="' + (padT + plotH / 2) + '" font-size="10" text-anchor="middle" fill="#767671" transform="rotate(-90 12 ' + (padT + plotH / 2) + ')">' + FIELD_LABEL[yKey] + '</text>' +
      '</svg>';
    renderLegend(rowsForColor);

    wireDragSelect(document.getElementById('scatter-svg'));
    container.querySelectorAll('circle.dot').forEach(function(c) {
      var idx = parseInt(c.getAttribute('data-idx'), 10);
      var h = currentHits[idx];
      c.addEventListener('click', function(e) {
        if (dragMoved) return;   // a drag that ended on top of a dot must not also navigate away
        window.location.href = hitUrl(h);
      });
      c.addEventListener('pointerenter', function(e) { showGraphTip(e, h, xKey, yKey); });
      c.addEventListener('pointermove', function(e) { positionGraphTip(e); });
      c.addEventListener('pointerleave', hideGraphTip);
    });
  }

  function showGraphTip(evt, h, xKey, yKey) {
    var tip = document.getElementById('graph-tip');
    tip.textContent = h.gene + ' / ' + h.cell_type + '\n' + h.hit_transcript +
      '\n' + FIELD_LABEL[xKey] + ': ' + h[xKey] + '\n' + FIELD_LABEL[yKey] + ': ' + h[yKey] +
      '\n(click to open, drag to multi-select)';
    tip.style.display = 'block';
    positionGraphTip(evt);
  }

  function positionGraphTip(evt) {
    var tip = document.getElementById('graph-tip');
    if (tip.style.display !== 'block') return;
    tip.style.left = (evt.clientX + 14) + 'px';
    tip.style.top = (evt.clientY + 14) + 'px';
  }

  function hideGraphTip() {
    document.getElementById('graph-tip').style.display = 'none';
  }

  var dragStart = null, dragMoved = false, dragRectEl = null;

  function svgPoint(svg, evt) {
    var pt = svg.createSVGPoint();
    pt.x = evt.clientX; pt.y = evt.clientY;
    return pt.matrixTransform(svg.getScreenCTM().inverse());
  }

  var dragHistoryPushed = false;

  function wireDragSelect(svg) {
    if (!svg) return;
    svg.addEventListener('pointerdown', function(e) {
      dragStart = svgPoint(svg, e);
      dragMoved = false;
      dragHistoryPushed = false;
      hideGraphTip();   // a drag shouldn't leave a stale tooltip from the point it started on
      svg.setPointerCapture(e.pointerId);
    });
    svg.addEventListener('pointermove', function(e) {
      if (!dragStart) return;
      var cur = svgPoint(svg, e);
      var dx = Math.abs(cur.x - dragStart.x), dy = Math.abs(cur.y - dragStart.y);
      if (dx < 3 && dy < 3 && !dragMoved) return;   // ignore sub-pixel jitter so plain clicks still work
      dragMoved = true;
      if (!dragHistoryPushed) {
        // Snapshot the pre-drag selection ONCE per drag, so "undo" reverts
        // this whole drag in one step rather than one point at a time.
        selectionHistory.push(Object.assign({}, selectedIdx));
        dragHistoryPushed = true;
      }

      var x0 = Math.min(dragStart.x, cur.x), x1 = Math.max(dragStart.x, cur.x);
      var y0 = Math.min(dragStart.y, cur.y), y1 = Math.max(dragStart.y, cur.y);

      if (!dragRectEl) {
        dragRectEl = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        dragRectEl.setAttribute('class', 'select-rect');
        svg.appendChild(dragRectEl);
      }
      dragRectEl.setAttribute('x', x0); dragRectEl.setAttribute('y', y0);
      dragRectEl.setAttribute('width', x1 - x0); dragRectEl.setAttribute('height', y1 - y0);

      plotPoints.forEach(function(p) {
        var inBox = p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1;
        var circle = svg.querySelector('circle[data-idx="' + p.idx + '"]');
        if (!circle) return;
        if (inBox) { selectedIdx[p.idx] = true; circle.classList.add('selected'); }
      });
      renderGraphSelectionList();
    });
    svg.addEventListener('pointerup', function(e) {
      dragStart = null;
      if (dragRectEl) { dragRectEl.remove(); dragRectEl = null; }
      // dragMoved intentionally stays true through the immediately-following
      // click event (browsers fire click right after pointerup on the same
      // target) so a drag that released over a dot doesn't also navigate;
      // it's reset on the NEXT pointerdown instead.
      setTimeout(function() { dragMoved = false; }, 0);
    });
  }

  function addSelectedToCart() {
    var indices = Object.keys(selectedIdx);
    if (!indices.length) return;
    Promise.all(indices.map(function(idxStr) {
      var h = currentHits[parseInt(idxStr, 10)];
      return fetch('/api/cart/items', {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({gene: h.gene, cell_type: h.cell_type, hit_enst: h.hit_enst, selected_drugs: []}),
      });
    })).then(function() { return fetch('/api/cart'); })
      .then(function(r) { return r.json(); })
      .then(function(data) {
        renderCart(data.items);
        clearPlotSelection();   // the batch has been committed to the cart
      });
  }

  function clearPlotSelection() {
    selectedIdx = {};
    selectionHistory = [];
    document.querySelectorAll('#scatter circle.dot.selected').forEach(function(c) { c.classList.remove('selected'); });
    renderGraphSelectionList();
  }

  function renderCart(items) {
    var el = document.getElementById('cart-items');
    el.innerHTML = '';
    items.forEach(function(item) {
      var div = document.createElement('div');
      div.className = 'cart-item';
      div.innerHTML = '<div><b>' + item.gene + '</b> / ' + item.cell_type + '</div>' +
        '<div class="subtle">' + (item.selected_drugs.join(', ') || 'no drugs selected') + '</div>' +
        '<button data-remove-gene="' + item.gene + '" data-remove-ct="' + item.cell_type + '" data-remove-enst="' + item.hit_enst + '">remove</button>';
      el.appendChild(div);
    });
    el.querySelectorAll('button[data-remove-gene]').forEach(function(btn) {
      btn.addEventListener('click', function() {
        fetch('/api/cart/items', {
          method: 'DELETE', headers: {'Content-Type': 'application/json'},
          body: JSON.stringify({
            gene: btn.getAttribute('data-remove-gene'),
            cell_type: btn.getAttribute('data-remove-ct'),
            hit_enst: btn.getAttribute('data-remove-enst'),
          }),
        }).then(function(r) { return r.json(); }).then(function(data) { renderCart(data.items); });
      });
    });
    document.getElementById('cart-count').textContent = items.length + ' in cart';
  }

  function loadCart() {
    fetch('/api/cart').then(function(r) { return r.json(); }).then(function(data) { renderCart(data.items); });
  }

  var pollTimer = null;
  function startExport() {
    fetch('/api/export', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({thresholds_used: thresholds})})
      .then(function(r) { return r.json(); })
      .then(function(data) {
        if (data.error) { document.getElementById('export-status').textContent = data.error; return; }
        pollJob(data.job_id);
      });
  }

  function pollJob(jobId) {
    if (pollTimer) clearInterval(pollTimer);
    function poll() {
      fetch('/api/jobs/' + jobId).then(function(r) { return r.json(); }).then(function(status) {
        renderJobStatus(jobId, status);
        if (status.status === 'done' || status.status === 'done_with_errors') clearInterval(pollTimer);
      });
    }
    poll();
    pollTimer = setInterval(poll, 3000);
  }

  function renderJobStatus(jobId, status) {
    var el = document.getElementById('export-status');
    var html = '<div><b>Job ' + jobId + '</b> — ' + status.status + '</div>';
    (status.items || []).forEach(function(item, idx) {
      html += '<div class="job-item">' + item.gene + ' / ' + item.cell_type + ': ' + item.stage;
      if (item.item_status === 'done') {
        html += ' <a href="/api/export/' + jobId + '/download/' + idx + '">[download zip]</a>' +
          ' (confidence: ' + (item.result ? item.result.structure_confidence : '?') + ')';
      }
      if (item.item_status === 'error') {
        html += ' <span style="color:#d03b3b">FAILED: ' + (item.error || '').slice(0, 200) + '</span>';
      }
      html += '</div>';
    });
    el.innerHTML = html;
  }

  document.getElementById('search').addEventListener('input', fetchHits);
  document.getElementById('tf-cutoff').addEventListener('change', function(e) {
    thresholds.tf_min_abs_delta_usage = parseFloat(e.target.value) || 0;
    fetchHits();
  });
  document.getElementById('dr-cutoff').addEventListener('change', function(e) {
    thresholds.dr_min_chembl_or_ot_phase = parseInt(e.target.value, 10) || 0;
    fetchHits();
  });
  document.getElementById('export-btn').addEventListener('click', startExport);
  document.getElementById('plot-add-selected-btn').addEventListener('click', addSelectedToCart);
  document.getElementById('plot-clear-selection-btn').addEventListener('click', clearPlotSelection);
  document.getElementById('plot-undo-btn').addEventListener('click', undoSelection);

  document.querySelectorAll('.tab-btn').forEach(function(btn) {
    btn.addEventListener('click', function() {
      var target = btn.getAttribute('data-tab');
      document.querySelectorAll('.tab-btn').forEach(function(b) { b.classList.toggle('active', b === btn); });
      document.getElementById('tab-table').hidden = target !== 'table';
      document.getElementById('tab-graph').hidden = target !== 'graph';
      document.getElementById('layout').classList.toggle('graph-tab-active', target === 'graph');
      if (target === 'graph') renderPlot();   // cheap; guards against any stale layout from being hidden
    });
  });

  document.querySelectorAll('#group-toggle .toggle-btn').forEach(function(btn) {
    btn.addEventListener('click', function() {
      var g = btn.getAttribute('data-group');
      groupFilter[g] = !groupFilter[g];
      btn.classList.toggle('active', groupFilter[g]);
      renderTable();
      renderPlot();
    });
  });

  initPlotSelectors();
  fetchHits();
  loadCart();
  renderGraphSelectionList();
})();
"""

INDEX_HTML = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Isthmus — master_surveyor</title>
<style>{_CSS}</style>
</head>
<body>
<div id="splash">
  <img src="/static/images/Isthmus_logo.png" alt="Isthmus">
  <div class="spinner"></div>
  <div class="label">loading hits&hellip;</div>
</div>
<div class="brand">
  <img src="/static/images/Isthmus_logo.png" alt="Isthmus logo">
  <h1>Isthmus</h1>
</div>
<div class="subtle">Interactive master_surveyor dossier — filter, select transcripts/drugs, export to BIOVIA Discovery Studio.</div>

<div class="layout" id="layout">
  <div>
    <div class="tab-bar">
      <button class="tab-btn active" id="tab-btn-table" data-tab="table">Table</button>
      <button class="tab-btn" id="tab-btn-graph" data-tab="graph">Graph Explorer</button>
    </div>

    <div class="tab-content" id="tab-table">
      <div class="card table-wrap">
        <table>
          <thead><tr>
            <th>Group</th><th>Gene</th><th>Cell type</th><th>Change</th><th>&Delta;pp</th>
            <th>ColabFold<br>canonical</th><th>ESMFold<br>canonical</th><th>ColabFold<br>isoform</th><th>ESMFold<br>isoform</th>
            <th>Drugs</th><th></th>
          </tr></thead>
          <tbody id="hits-body"></tbody>
        </table>
      </div>
    </div>

    <div class="tab-content" id="tab-graph" hidden>
      <div class="card">
        <div class="plot-toolbar">
          <div class="filter-group">
            <label>X axis</label>
            <select id="plot-x"></select>
          </div>
          <div class="filter-group">
            <label>Y axis</label>
            <select id="plot-y"></select>
          </div>
        </div>
        <div id="scatter-wrap">
          <div id="scatter"></div>
        </div>
        <div class="subtle" style="margin-top:8px;">Drag to select multiple points &middot; click one point to open its hit page &middot; hover for details.</div>
      </div>
    </div>
  </div>

  <div>
    <div class="card">
      <h2 style="font-size:13px;margin:0 0 8px;">Filters</h2>
      <div class="filter-group" style="margin-bottom:10px;">
        <label>Search gene</label>
        <input type="text" id="search" placeholder="e.g. KLC1">
      </div>
      <div class="filter-group" style="margin-bottom:10px;">
        <label>trial_failure: min |&Delta; usage|</label>
        <input type="number" id="tf-cutoff" step="0.01" value="0">
      </div>
      <div class="filter-group" style="margin-bottom:10px;">
        <label>drug_repurposing: min ChEMBL/OT phase</label>
        <input type="number" id="dr-cutoff" step="1" min="0" max="4" value="0">
      </div>
      <div class="filter-group" style="margin-bottom:10px;">
        <label>Group type</label>
        <div class="toggle-group" id="group-toggle">
          <button type="button" class="toggle-btn active" data-group="trial_failure_candidate">Trial failure</button>
          <button type="button" class="toggle-btn active" data-group="drug_repurposing_candidate">Drug repurposing</button>
          <button type="button" class="toggle-btn active" data-group="novel_target_candidate">Novel target</button>
        </div>
        <div class="subtle" style="margin-top:4px;">Novel target is out of master_surveyor's scope (excluded upstream), so that toggle won't surface any hits here.</div>
      </div>
      <div class="count" id="count"></div>
    </div>
    <div class="card">
      <h2 style="font-size:13px;margin:0 0 8px;">Color by</h2>
      <div class="filter-group" style="margin-bottom:10px;">
        <label>Field</label>
        <select id="plot-color-key"></select>
      </div>
      <div class="filter-group" id="plot-color-scale-group" style="margin-bottom:10px;">
        <label>Color scale (numeric fields only)</label>
        <select id="plot-color-scale"></select>
      </div>
      <div id="plot-legend"></div>
    </div>
    <div class="card">
      <h2 style="font-size:13px;margin:0 0 8px;">Graph selection <span class="count" id="graph-sel-count"></span></h2>
      <div id="graph-selection-items"></div>
      <div class="plot-actions">
        <button id="plot-add-selected-btn">+ Add to cart</button>
        <button id="plot-undo-btn">undo</button>
        <button id="plot-clear-selection-btn">clear</button>
      </div>
    </div>
    <div class="card">
      <h2 style="font-size:13px;margin:0 0 8px;">Cart <span class="count" id="cart-count"></span></h2>
      <div id="cart-items"></div>
      <button id="export-btn" style="margin-top:10px;">Export selected &rarr; Discovery Studio</button>
    </div>
    <div class="card" id="export-status"></div>
  </div>
</div>

<div class="graph-tip" id="graph-tip"></div>

<script>{_JS}</script>
</body>
</html>
"""
