"""Live per-hit detail page -- reuses dossier/render.py's existing sections
(header, stacked bar, donor dots, sequence diff) as-is, and adds two things
the static dossier doesn't have: a drug PICKER (checkboxes -> add to cart,
replacing the read-only drug_table_section) and a structure-confidence
readout once the hit has actually been through an export job.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier import data
from dossier.render import (
    THEME_TOGGLE_HTML, THEME_TOGGLE_JS, _BASE_CSS, _TIP_JS, _esc, _font_face_css,
    donor_dot_section, header_section, sequence_section, stacked_bar_section,
)
from master_surveyor.m3_export import hit_folder_name
from master_surveyor.config import DOCKING_DIR


def _drug_checklist(label: str, names: list[str]) -> str:
    if not names:
        return ""
    boxes = "".join(
        f'<label style="display:block; font-size:12px; margin:3px 0;">'
        f'<input type="checkbox" class="drug-check" value="{_esc(n)}"> {_esc(n)}</label>'
        for n in names
    )
    return f'<div style="margin-top:8px;"><div class="subtle" style="margin-bottom:4px;">{_esc(label)}</div>{boxes}</div>'


def drug_picker_section(row: pd.Series) -> str:
    drug_names = [d for d in str(row.get("drug_names") or "").split("|") if d]
    ot_names = [d for d in str(row.get("ot_drug_names") or "").split("|") if d]
    all_names = list(dict.fromkeys(drug_names + ot_names))

    body = (
        _drug_checklist("ChEMBL / Open Targets candidates", all_names)
        if all_names else '<div class="subtle">no named drug candidates on record</div>'
    )
    return f"""
<div class="card">
  <h2>Select drugs &amp; add to cart</h2>
  {body}
  <button id="add-to-cart-btn" style="margin-top:10px;">+ Add to cart</button>
  <div id="cart-add-status" class="subtle" style="margin-top:6px;"></div>
</div>"""


def structure_viewer_section(gene: str, cell_type: str, hit_enst: str) -> str:
    """Raw viewer over whatever ColabFold/ESMFold structures are already
    cached for this hit's canonical + every ranked alt isoform -- reads
    directly from STRUCTURE_CACHE_DIR via /api/structures and /api/structure_pdb,
    with no dependency on the cart/export flow (unlike structure_confidence_section
    below, which needs the QC pass that only runs post-export).
    """
    return f"""
<div class="card">
  <h2>Structures</h2>
  {_structure_viewer_html()}
</div>"""


def structure_confidence_section(gene: str, cell_type: str, transcript: str) -> str:
    folder = DOCKING_DIR / hit_folder_name(gene, cell_type, transcript)
    manifest_path = folder / "manifest.json"
    if not manifest_path.exists():
        return """
<div class="card">
  <h2>Structure confidence</h2>
  <div class="subtle">Not yet exported -- add this hit to the cart and run an export to see
  numeric ColabFold/ESMFold confidence metrics over the altered region (the structures
  themselves are viewable above as soon as they're folded, no export needed).</div>
</div>"""

    manifest = json.loads(manifest_path.read_text())
    conf = manifest.get("structure_confidence", {})
    region = conf.get("region_confidence", {})
    sup = conf.get("superposition", {})
    verdict = conf.get("verdict", "unknown")
    verdict_color = {"high": "var(--status-good)", "medium": "#c98a1a", "low": "var(--status-critical)"}.get(verdict, "var(--muted)")

    rows = [
        ("Verdict", f'<span style="color:{verdict_color}; font-weight:700; text-transform:uppercase;">{_esc(verdict)}</span>'),
        ("Region pLDDT (altered span only)", region.get("region_mean_plddt")),
        ("Whole-protein pLDDT", region.get("whole_protein_mean_plddt")),
        ("Region-vs-rest PAE", region.get("region_vs_rest_mean_pae")),
        ("Ensemble spread (RMSD, &Aring;)", conf.get("ensemble_spread_rmsd")),
        ("ColabFold vs ESMFold (RMSD, &Aring;)", conf.get("colabfold_vs_esmfold_rmsd")),
        ("Alt-vs-canonical, changed region (mean local dist, &Aring;)", sup.get("changed_region_mean_local_dist")),
        ("Alt-vs-canonical, outside changed region (mean local dist, &Aring;)", sup.get("outside_span_mean_local_dist")),
    ]
    rows_html = "".join(
        f'<div style="display:flex; justify-content:space-between; padding:4px 0; border-bottom:1px solid var(--gridline); font-size:12px;">'
        f'<span class="subtle">{label}</span><span>{value if value is not None else "&mdash;"}</span></div>'
        for label, value in rows
    )

    return f"""
<div class="card">
  <h2>Structure confidence</h2>
  {rows_html}
  <div class="subtle" style="margin-top:8px;">Canonical structure source: {_esc(manifest.get('canonical_structure_source', '?'))}</div>
</div>"""


def _structure_viewer_html() -> str:
    return """
<div class="ngl-card">
  <div class="viewer-split">
    <div class="viewer-panel">
      <div class="viewer-panel-header">Canonical</div>
      <div class="viewer-toolbar">
        <label class="viewer-select-label">structure
          <select id="ngl-canonical-select"></select>
        </label>
        <label class="viewer-select-label">style
          <select id="ngl-canonical-style-select">
            <option value="cartoon">cartoon</option>
            <option value="spacefill">space-filling</option>
          </select>
        </label>
      </div>
      <div class="viewer-canvas"><div class="viewer-mount" id="ngl-viewer-canonical"></div><div class="viewer-loading">loading structure&hellip;</div></div>
      <div class="viewer-toolbar">
        <button type="button" data-ngl-action="spin" data-ngl-panel="canonical">spin</button>
        <button type="button" data-ngl-action="reset" data-ngl-panel="canonical">reset view</button>
      </div>
    </div>
    <div class="viewer-panel">
      <div class="viewer-panel-header">Alternative isoform</div>
      <div class="viewer-toolbar">
        <label class="viewer-select-label">isoform
          <select id="ngl-alt-rank-select"></select>
        </label>
        <label class="viewer-select-label">method
          <select id="ngl-alt-method-select"></select>
        </label>
        <label class="viewer-select-label">style
          <select id="ngl-alt-style-select">
            <option value="cartoon">cartoon</option>
            <option value="spacefill">space-filling</option>
          </select>
        </label>
      </div>
      <div class="viewer-canvas"><div class="viewer-mount" id="ngl-viewer-alt"></div><div class="viewer-loading">loading structure&hellip;</div></div>
      <div class="viewer-toolbar">
        <button type="button" data-ngl-action="spin" data-ngl-panel="alt">spin</button>
        <button type="button" data-ngl-action="reset" data-ngl-panel="alt">reset view</button>
      </div>
    </div>
  </div>
  <div class="plddt-legend">
    <div class="bar"></div>
    <div class="ticks"><span>low confidence</span><span>high confidence</span></div>
  </div>
</div>"""


_STRUCTURE_VIEWER_CSS = """
.ngl-card { margin-bottom: 16px; }
.viewer-split { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; }
@media (max-width: 860px) { .viewer-split { grid-template-columns: 1fr; } }
.viewer-panel-header { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; margin-bottom: 6px; }
.viewer-toolbar { display: flex; gap: 6px; flex-wrap: wrap; margin-bottom: 8px; }
.viewer-toolbar button {
  font-family: inherit; font-size: 12px; cursor: pointer;
  background: var(--page-plane); color: var(--text-primary);
  border: 1px solid var(--border); border-radius: 0; padding: 5px 10px;
}
.viewer-toolbar button:hover { border-color: var(--accent-control); color: var(--accent-control); }
.viewer-toolbar button.is-active { background: var(--accent-control); color: var(--surface-1); border-color: var(--accent-control); }
.viewer-select-label { display: flex; align-items: center; gap: 6px; font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.03em; }
.viewer-select-label select { font-family: inherit; font-size: 12px; text-transform: none; letter-spacing: normal; background: var(--page-plane); color: var(--text-primary); border: 1px solid var(--border); border-radius: 0; padding: 3px 6px; }
.viewer-canvas {
  position: relative; height: 320px; overflow: hidden;
  background: var(--page-plane); border: 1px solid var(--border);
}
/* NGL's Stage breaks irrecoverably if its mount element has any child
   present at construction time, so the loading label lives as a sibling
   overlay rather than nesting inside the element it owns. */
.viewer-mount { position: absolute; inset: 0; }
.viewer-canvas .viewer-loading {
  position: absolute; inset: 0; display: flex; align-items: center; justify-content: center;
  font-size: 11px; color: var(--muted); letter-spacing: 0.04em; text-transform: uppercase;
  pointer-events: none;
}
.plddt-legend { display: flex; flex-direction: column; gap: 4px; margin-top: 8px; }
.plddt-legend .bar {
  height: 6px; border-radius: 100px;
  background: linear-gradient(90deg, var(--accent-ad) 0%, #c98a1a 45%, #d6c22a 60%, var(--accent-control) 100%);
}
.plddt-legend .ticks { display: flex; justify-content: space-between; font-size: 10px; color: var(--muted); }
"""

_STRUCTURE_VIEWER_JS = """
(function () {
  var canonMount = document.getElementById("ngl-viewer-canonical");
  var altMount = document.getElementById("ngl-viewer-alt");
  if (!canonMount || !altMount || !window.NGL) return;

  var gene = window.__HIT__.gene, cellType = window.__HIT__.cell_type, hitEnst = window.__HIT__.hit_enst;
  var listUrl = "/api/structures/" + encodeURIComponent(gene) + "/" + encodeURIComponent(cellType)
    + "?hit_enst=" + encodeURIComponent(hitEnst);

  var canonicalSelect = document.getElementById("ngl-canonical-select");
  var altRankSelect = document.getElementById("ngl-alt-rank-select");
  var altMethodSelect = document.getElementById("ngl-alt-method-select");
  var canonicalStyleSelect = document.getElementById("ngl-canonical-style-select");
  var altStyleSelect = document.getElementById("ngl-alt-style-select");

  var altIsoforms = [];

  // Both panels render with the same pLDDT confidence coloring -- they no
  // longer share one stage (that's what the red alt overlay used to
  // distinguish), so identical coloring makes them directly comparable.
  var plddtScheme = NGL.ColormakerRegistry.addScheme(function () {
    this.atomColor = function (atom) {
      var t = Math.max(0, Math.min(1, (atom.bfactor - 50) / (95 - 50)));
      var stops = [[227, 73, 72], [214, 194, 42], [42, 120, 214]];
      var seg = t < 0.5 ? [stops[0], stops[1], t / 0.5] : [stops[1], stops[2], (t - 0.5) / 0.5];
      var r = seg[0][0] + (seg[1][0] - seg[0][0]) * seg[2];
      var g = seg[0][1] + (seg[1][1] - seg[0][1]) * seg[2];
      var b = seg[0][2] + (seg[1][2] - seg[0][2]) * seg[2];
      return (Math.round(r) << 16) | (Math.round(g) << 8) | Math.round(b);
    };
  }, "plddt");

  function pdbUrl(opt) {
    var u = "/api/structure_pdb?seq_hash=" + encodeURIComponent(opt.seq_hash) + "&method=" + encodeURIComponent(opt.method);
    if (opt.file) u += "&file=" + encodeURIComponent(opt.file);
    return u;
  }

  function fillSelect(select, options, valueFn, labelFn) {
    select.innerHTML = "";
    options.forEach(function (o, i) {
      var el = document.createElement("option");
      el.value = String(i);
      el.textContent = labelFn(o);
      select.appendChild(el);
    });
    select.disabled = options.length === 0;
    if (options.length === 0) {
      var el = document.createElement("option");
      el.textContent = "not folded yet";
      select.appendChild(el);
    }
  }

  // Each panel owns an independent NGL.Stage so canonical and alt render
  // (and spin/reset) side by side without one hiding behind the other.
  function makePanel(mount) {
    var stage = null, comp = null, spinning = false, anyLoaded = false, style = "cartoon";

    function clearLoading() {
      var l = mount.parentElement && mount.parentElement.querySelector(".viewer-loading");
      if (l) l.remove();
    }
    function markLoaded() {
      if (!anyLoaded) { anyLoaded = true; clearLoading(); }
    }
    function showError(msg) {
      if (anyLoaded) return;
      var l = mount.parentElement && mount.parentElement.querySelector(".viewer-loading");
      if (l) { l.textContent = msg; l.style.display = "flex"; }
    }
    function ensureStage() {
      if (stage) return stage;
      stage = new NGL.Stage(mount.id, { backgroundColor: "white" });
      // See the CSS comment above the mount element: only the DOM
      // background needs overriding post-construction -- NGL always clears
      // the WebGL canvas itself with alpha 0.
      stage.viewer.renderer.domElement.style.background = "transparent";
      var resizeObs = new ResizeObserver(function () { stage.handleResize(); });
      resizeObs.observe(mount);
      return stage;
    }

    // quality:"low" avoids a WebGL context-loss bug in NGL's default impostor
    // shader path when combined with per-residue coloring, seen on
    // software/remote-rendered GPUs -- plausible here since this is a
    // single-box internal tool that may be viewed over a remote session.
    function applyRepresentation() {
      if (!comp) return;
      comp.removeAllRepresentations();
      comp.addRepresentation(style, { color: plddtScheme, quality: "low" });
    }

    function load(opt) {
      if (!opt) { if (comp) { stage.removeComponent(comp); comp = null; } return; }
      ensureStage();
      fetch(pdbUrl(opt)).then(function (r) { return r.json(); }).then(function (data) {
        if (!data.pdb) throw new Error("no pdb");
        if (comp) { stage.removeComponent(comp); comp = null; }
        var blob = new Blob([data.pdb], { type: "text/plain" });
        return stage.loadFile(blob, { ext: "pdb" });
      }).then(function (c) {
        if (!c) return;
        comp = c;
        applyRepresentation();
        comp.autoView();
        markLoaded();
      }).catch(function () { showError("structure not available"); });
    }

    // Swaps representation in place on the already-loaded component --
    // switching style shouldn't re-fetch the PDB or reset the camera.
    function setStyle(newStyle) {
      style = newStyle;
      applyRepresentation();
    }

    function spinToggle(btn) {
      spinning = !spinning;
      if (stage) stage.setSpin(spinning ? [0, 1, 0] : null, spinning ? 0.01 : 0);
      if (btn) btn.classList.toggle("is-active", spinning);
    }

    function reset() {
      spinning = false;
      if (stage) { stage.setSpin(null, 0); stage.autoView(); }
    }

    return { load: load, setStyle: setStyle, spinToggle: spinToggle, reset: reset, error: showError };
  }

  var canonPanel = makePanel(canonMount);
  var altPanel = makePanel(altMount);

  canonicalStyleSelect.addEventListener("change", function () { canonPanel.setStyle(canonicalStyleSelect.value); });
  altStyleSelect.addEventListener("change", function () { altPanel.setStyle(altStyleSelect.value); });

  function currentAltRank() {
    return altIsoforms[Number(altRankSelect.value)];
  }

  function onAltRankChange() {
    var rank = currentAltRank();
    fillSelect(altMethodSelect, rank ? rank.options : [], null, function (o) { return o.label; });
    onAltMethodChange();
  }

  function onAltMethodChange() {
    var rank = currentAltRank();
    var opt = rank && rank.options[Number(altMethodSelect.value)];
    altPanel.load(opt);
  }

  fetch(listUrl).then(function (r) {
    if (!r.ok) throw new Error("hit not found");
    return r.json();
  }).then(function (data) {
    altIsoforms = data.alt_isoforms || [];
    var hasCanonical = data.canonical.options.length > 0;
    var hasAnyAlt = altIsoforms.some(function (r) { return r.options.length > 0; });

    fillSelect(canonicalSelect, data.canonical.options, null, function (o) { return o.label; });
    canonicalSelect.addEventListener("change", function () {
      canonPanel.load(data.canonical.options[Number(canonicalSelect.value)]);
    });
    if (hasCanonical) canonPanel.load(data.canonical.options[0]);
    else canonPanel.error("no canonical structure folded yet");

    fillSelect(altRankSelect, altIsoforms, null, function (r) {
      var tag = r.is_gate_driver ? " \\u2605" : "";
      var delta = (typeof r.usage_delta === "number") ? " (\\u0394" + r.usage_delta.toFixed(2) + ")" : "";
      return "rank " + r.alt_rank + tag + delta;
    });
    altRankSelect.addEventListener("change", onAltRankChange);
    altMethodSelect.addEventListener("change", onAltMethodChange);
    if (altIsoforms.length) onAltRankChange();
    if (!hasAnyAlt) altPanel.error("no alt isoform folded yet");
  }).catch(function () {
    canonPanel.error("structures not available");
    altPanel.error("structures not available");
  });

  document.addEventListener("click", function (e) {
    var btn = e.target.closest("button[data-ngl-action]");
    if (!btn) return;
    var action = btn.getAttribute("data-ngl-action");
    var panel = btn.getAttribute("data-ngl-panel") === "alt" ? altPanel : canonPanel;
    if (action === "spin") panel.spinToggle(btn);
    else if (action === "reset") panel.reset();
  });
})();
"""


_CART_JS = """
(function() {
  var gene = window.__HIT__.gene, cellType = window.__HIT__.cell_type, hitEnst = window.__HIT__.hit_enst;
  var btn = document.getElementById('add-to-cart-btn');
  if (!btn) return;
  btn.addEventListener('click', function() {
    var checked = document.querySelectorAll('.drug-check:checked');
    var selected = Array.prototype.map.call(checked, function(c) { return c.value; });
    fetch('/api/cart/items', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({gene: gene, cell_type: cellType, hit_enst: hitEnst, selected_drugs: selected}),
    }).then(function(r) { return r.json(); }).then(function() {
      document.getElementById('cart-add-status').textContent = 'Added to cart (' + selected.length + ' drug(s) selected).';
    });
  });
})();
"""


def render_detail(gene: str, cell_type: str, hit_enst: str | None = None) -> str:
    rows = data.get_hit_rows(gene, cell_type, hit_enst)
    hit_row = rows.iloc[0]
    transcript = hit_row.get("hit_transcript_name") or hit_row["hit_ENST_ID"]

    body = (
        header_section(hit_row)
        + stacked_bar_section(gene, cell_type, hit_row["hit_ENST_ID"])
        + donor_dot_section(gene, cell_type, hit_row["hit_ENST_ID"])
        + sequence_section(rows)
        + drug_picker_section(hit_row)
        + structure_viewer_section(gene, cell_type, hit_row["hit_ENST_ID"])
        + structure_confidence_section(gene, cell_type, transcript)
    )

    hit_json = json.dumps({
        "gene": gene, "cell_type": cell_type, "hit_enst": hit_row["hit_ENST_ID"], "transcript": transcript,
    })

    # No font-face override needed here -- render.py's own _BASE_CSS already
    # sets .viz-root to Tamzen (with the face embedded by _font_face_css()),
    # which is exactly the font this page should use too.
    isthmus_font_css = """
.isthmus-brand { display: flex; align-items: center; gap: 8px; margin-bottom: 12px; }
.isthmus-brand img { width: 24px; height: 24px; }
"""

    return f"""<!doctype html>
<html><head><meta charset="utf-8"/>
<title>{_esc(gene)} &middot; {_esc(cell_type)} &mdash; Isthmus</title>
<style>{_font_face_css()}{_BASE_CSS}{isthmus_font_css}{_STRUCTURE_VIEWER_CSS}</style>
</head><body>
{THEME_TOGGLE_HTML}
<div class="viz-root">
<div class="isthmus-brand">
  <img src="/static/images/Isthmus_logo.png" alt="Isthmus">
  <a class="subtle" href="/">&larr; back to Isthmus index</a>
</div>
{body}
</div>
<div class="tip"></div>
<script>window.__HIT__ = {hit_json};</script>
<script>{THEME_TOGGLE_JS}</script>
<script>{_TIP_JS}</script>
<script>{_CART_JS}</script>
<script src="/static/vendor/ngl.js"></script>
<script>{_STRUCTURE_VIEWER_JS}</script>
</body></html>"""
