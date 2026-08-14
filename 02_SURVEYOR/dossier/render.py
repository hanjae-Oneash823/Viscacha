"""DOSSIER — self-contained HTML report for one (gene, cell_type) candidate.

Monochrome base + two accents (Control=blue, AD=red), used consistently
across every chart. Inline SVG marks, one shared JS tooltip layer (data-tip
attributes + a single pointermove listener), Tamzen bitmap monospace type
embedded as base64 so the file opens standalone with no network dependency.
"""

from __future__ import annotations

import base64
import html
import math

import pandas as pd

from dossier import data, palette, sequence_diff
from dossier.config import ASSETS_DIR, CANDIDATE_TYPE_LABELS, STRUCTURE_CACHE_DIR


def _esc(s) -> str:
    return html.escape("" if pd.isna(s) else str(s))


def _tip_esc(s: str) -> str:
    """Escape a multi-line tooltip string for an HTML attribute, turning
    real newlines into &#10; -- the JS layer reads the attribute back with
    getAttribute() (not innerHTML), so &#10; decodes to an actual line
    break in the resulting textContent, which is what `white-space: pre`
    on .tip needs to render lines instead of a literal backslash-n.
    """
    return html.escape(s).replace("\n", "&#10;")


def _canonical_name(canonical_enst: str) -> str:
    tx_map = data.load_tx_id_map()
    hit = tx_map.loc[tx_map["ENST_ID"] == canonical_enst, "transcript_name"]
    return hit.iloc[0] if len(hit) else canonical_enst


# ---------------------------------------------------------------------------
# Fonts / shell
# ---------------------------------------------------------------------------

def _font_face_css() -> str:
    faces = []
    for fname, weight in [("Tamzen8x16r.ttf", 400), ("Tamzen8x16b.ttf", 700)]:
        b64 = base64.b64encode((ASSETS_DIR / "fonts" / fname).read_bytes()).decode()
        faces.append(f"""
@font-face {{
  font-family: "Tamzen";
  src: url(data:font/ttf;base64,{b64}) format("truetype");
  font-weight: {weight};
  font-style: normal;
  font-display: swap;
}}""")
    return "\n".join(faces)


_BASE_CSS = f"""
:root {{ color-scheme: light; }}
.viz-root {{
  --surface-1: {palette.SURFACE}; --page-plane: {palette.PAGE_PLANE};
  --text-primary: {palette.TEXT_PRIMARY}; --text-secondary: {palette.TEXT_SECONDARY};
  --muted: {palette.MUTED}; --gridline: {palette.GRIDLINE}; --baseline: {palette.BASELINE};
  --border: {palette.BORDER}; --match-gray: {palette.MATCH_GRAY}; --bar-black: #000000;
  --accent-control: {palette.ACCENT_CONTROL_HEX}; --accent-ad: {palette.ACCENT_AD_HEX};
  --pure-red: {palette.PURE_RED_HEX};
  --status-good: #0ca30c; --status-critical: #d03b3b; --glow: none;
}}
/* Dark theme is a deliberate, distinct world (pure black + neon), not an
   inverted tint of the light one -- see palette.py's dark-theme note. */
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{ {palette.NEON_DARK_CSS_VARS} }}
}}
:root[data-theme="dark"] .viz-root {{ {palette.NEON_DARK_CSS_VARS} }}
.viz-root {{
  font-family: "Tamzen", ui-monospace, monospace;
  background: var(--page-plane); color: var(--text-primary);
  max-width: 1400px; margin: 0 auto; padding: 24px 20px 64px;
  line-height: 1.5;
}}
.card {{
  background: var(--surface-1); border: 1px solid var(--border);
  border-radius: 0; padding: 20px 24px; margin-bottom: 20px;
}}
h1 {{ font-size: 28px; margin: 0 0 4px; letter-spacing: 0.02em; text-shadow: var(--glow); }}
h2 {{ font-size: 15px; text-transform: uppercase; letter-spacing: 0.08em;
      color: var(--text-secondary); margin: 0 0 14px; font-weight: 700; }}
.badge {{ display: inline-block; padding: 3px 10px; border-radius: 0;
          border: 1px solid var(--border); font-size: 12px; margin-right: 8px; text-shadow: var(--glow); }}
.badge.control {{ color: var(--accent-control); border-color: var(--accent-control); }}
.badge.ad {{ color: var(--accent-ad); border-color: var(--accent-ad); }}
.badge.change-decrease {{ color: var(--accent-ad); border-color: var(--accent-ad); }}
.badge.change-increase {{ color: var(--accent-control); border-color: var(--accent-control); }}
.badge.change-same {{ color: var(--text-primary); border-color: var(--text-primary); }}
.subtle {{ color: var(--text-secondary); font-size: 13px; }}
.muted  {{ color: var(--muted); }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
td, th {{ text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--gridline); }}
th {{ color: var(--text-secondary); font-weight: 700; }}
tr:last-child td {{ border-bottom: none; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 12px; margin-top: 10px; font-size: 12px; color: var(--text-secondary); }}
.legend .swatch {{ display: inline-block; width: 10px; height: 10px; border-radius: 0; margin-right: 5px; vertical-align: -1px; }}
.cell-type-block {{ margin-bottom: 22px; }}
.cell-type-block:last-child {{ margin-bottom: 0; }}
.ct-label {{ font-size: 13px; font-weight: 700; margin-bottom: 8px; }}
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
"""

# Shared across every page this package renders (dossiers + the index
# homepage) -- one toggle, one storage key, so flipping the theme on one
# page and then opening another (same origin, e.g. via `python -m
# http.server`) keeps the choice. Defined once here so generate_index.py
# imports it instead of hand-syncing a second copy.
THEME_TOGGLE_HTML = '<button class="theme-toggle" id="theme-toggle" title="Toggle theme" aria-label="Toggle theme">&#9788;</button>'
THEME_TOGGLE_JS = """
(function() {
  var root = document.documentElement;
  var btn = document.getElementById('theme-toggle');
  var stored = localStorage.getItem('dossier-theme');
  if (stored) { root.setAttribute('data-theme', stored); }
  function current() {
    return root.getAttribute('data-theme') ||
      (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  }
  function updateIcon() { btn.innerHTML = current() === 'dark' ? '&#9728;' : '&#9788;'; }
  updateIcon();
  btn.addEventListener('click', function() {
    var next = current() === 'dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    localStorage.setItem('dossier-theme', next);
    updateIcon();
  });
})();
"""

_TIP_JS = """
(function() {
  var tip = document.querySelector('.tip');
  document.addEventListener('pointermove', function(e) {
    var t = e.target.closest('[data-tip]');
    if (!t) { tip.style.display = 'none'; return; }
    tip.textContent = t.getAttribute('data-tip');
    tip.style.display = 'block';
    var x = e.clientX + 14, y = e.clientY + 14;
    if (x + 260 > window.innerWidth) x = e.clientX - 274;
    tip.style.left = x + 'px'; tip.style.top = y + 'px';
  });
  document.addEventListener('pointerdown', function(e) {
    var t = e.target.closest('[data-tip]');
    tip.style.display = t ? 'block' : 'none';
    if (t) { tip.textContent = t.getAttribute('data-tip'); }
  });
})();
"""


# ---------------------------------------------------------------------------
# Section 1 — header
# ---------------------------------------------------------------------------

def _structure_image_html(uniprot_acc: str) -> str:
    """Canonical AlphaFold cartoon, for visual purposes only -- pre-rendered
    and cached by fetch_structures.py (per uniprot_acc, one gene has one
    canonical structure regardless of which cell type or alt transcript
    triggered the hit). Silently omitted if the accession was never fetched
    (script hasn't run, no AlphaFold entry, or fetch failed) -- a dossier
    never depends on this cache existing.
    """
    if not uniprot_acc or pd.isna(uniprot_acc):
        return ""
    png_path = STRUCTURE_CACHE_DIR / f"{uniprot_acc}.png"
    if not png_path.exists():
        return ""
    b64 = base64.b64encode(png_path.read_bytes()).decode()
    # Caption overlays the image (bottom-left) rather than sitting below it,
    # so the thumbnail can run larger without growing the header's height.
    # Fixed dark-on-white regardless of site theme -- the render itself is
    # always on a white background (3Dmol.js render, not theme-aware), so a
    # dark scrim + white text reads reliably in both light and dark mode.
    return f"""
<div style="flex-shrink:0; position:relative; width:260px; line-height:0;">
  <img src="data:image/png;base64,{b64}" alt="AlphaFold canonical structure"
       style="width:260px; height:auto; display:block; border:1px solid var(--border);"/>
  <div style="position:absolute; left:0; bottom:0; padding:4px 8px; font-size:10px; line-height:1.4;
              background:rgba(0,0,0,0.6); color:#ffffff;">
    AlphaFold model<br/>canonical &middot; visual only
  </div>
</div>"""


def header_section(hit_row: pd.Series) -> str:
    gene = hit_row["gene_name"]
    hit_enst = hit_row["hit_ENST_ID"]
    hit_name = hit_row["hit_transcript_name"]
    canonical_enst = hit_row["canonical_enst"]
    is_mane = hit_enst == canonical_enst
    group = hit_row["master_group"]
    group_label = CANDIDATE_TYPE_LABELS.get(group, group)

    canonical_line = ""
    if not is_mane:
        canonical_line = (
            f'<div class="subtle">canonical (MANE Select): '
            f'<strong>{_esc(_canonical_name(canonical_enst))}</strong> '
            f'({_esc(canonical_enst)})</div>'
        )

    padj = hit_row.get("chi_padj")
    ctrl_pct = float(hit_row.get("Control") or 0) * 100
    ad_pct = float(hit_row.get("AD") or 0) * 100
    delta_pp = ad_pct - ctrl_pct
    stats_line = (
        f'<div class="subtle" style="margin-top:8px;">'
        f'DTU adjusted p-value: <strong>{padj:.2e}</strong>'
        f'&nbsp;&middot;&nbsp; usage '
        f'<strong style="color:{palette.ACCENT_CONTROL};">{ctrl_pct:.1f}%</strong> (Control) &rarr; '
        f'<strong style="color:{palette.ACCENT_AD};">{ad_pct:.1f}%</strong> (AD)'
        f'&nbsp;(&Delta; {delta_pp:+.1f} pp)</div>'
        if pd.notna(padj) else ""
    )

    structure_image = _structure_image_html(hit_row.get("uniprot_acc"))

    return f"""
<div class="card">
  <div style="display:flex; justify-content:space-between; gap:20px;">
    <div>
      <h1>{_esc(gene)}</h1>
      <div class="subtle" style="margin-bottom:10px;">
        transcript: <strong>{_esc(hit_name)}</strong> ({_esc(hit_enst)})
        &nbsp;&middot;&nbsp; cell type: <strong>{_esc(hit_row['cell_type'])}</strong>
      </div>
      {canonical_line}
      {stats_line}
      <div style="margin-top:12px;">
        <span class="badge {'control' if is_mane else 'ad'}">
          {'MANE canonical' if is_mane else 'alternative to MANE canonical'}
        </span>
        <span class="badge">{_esc(group_label)}</span>
      </div>
    </div>
    {structure_image}
  </div>
</div>"""


# ---------------------------------------------------------------------------
# Section 2 — stacked isoform-usage bars, per cell type
# ---------------------------------------------------------------------------

def _isoform_positions(segments: list[tuple[str, float, str, bool]], width: float) -> dict[str, tuple[float, float, str]]:
    """name -> (x_start, x_end, color) for every segment with usage > 0."""
    x = 0.0
    pos: dict[str, tuple[float, float, str]] = {}
    for name, frac, color, _is_hit in segments:
        if frac > 0:
            pos[name] = (x, x + frac * width, color)
        x += frac * width
    return pos


def _stacked_bar_layer(y: float, bar_h: float, segments: list[tuple[str, float, str, bool]], width: int) -> str:
    x = 0.0
    gap = 2
    parts = [f'<rect x="0" y="{y}" width="{width}" height="{bar_h}" fill="var(--gridline)" opacity="0.5"/>']
    for name, frac, color, is_hit in segments:
        w = max(frac * width - gap, 0)
        if w > 0:
            tip = _tip_esc(f"{name}{' (hit transcript)' if is_hit else ''}\n{frac * 100:.1f}%")
            stroke = 'stroke="var(--text-primary)" stroke-width="2"' if is_hit else ""
            parts.append(
                f'<rect data-tip="{tip}" '
                f'x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{color}" {stroke}/>'
            )
        x += frac * width
    return "".join(parts)


def _stacked_bars_svg(ctrl_segments: list[tuple[str, float, str, bool]],
                       ad_segments: list[tuple[str, float, str, bool]], width: int) -> str:
    """Control bar, a transitioning-fill ribbon per isoform present in BOTH
    conditions, then the AD bar -- same bezier-ribbon technique as the
    sequence-alignment diagram, colored per isoform instead of match/mismatch.
    An isoform only in one condition (usage 0 in the other) gets no ribbon --
    there's nothing to transition from/to.
    """
    bar_h, ribbon_h = 24, 44
    ctrl_y, ad_y = 0, bar_h + ribbon_h
    height = ad_y + bar_h

    ctrl_pos = _isoform_positions(ctrl_segments, width)
    ad_pos = _isoform_positions(ad_segments, width)
    y_top, y_bot = ctrl_y + bar_h, ad_y
    y_mid = (y_top + y_bot) / 2

    ribbons = []
    for name in set(ctrl_pos) & set(ad_pos):
        cx1, cx2, color = ctrl_pos[name]
        ax1, ax2, _ = ad_pos[name]
        path = (f"M {cx1:.1f},{y_top} C {cx1:.1f},{y_mid:.1f} {ax1:.1f},{y_mid:.1f} {ax1:.1f},{y_bot} "
                f"L {ax2:.1f},{y_bot} C {ax2:.1f},{y_mid:.1f} {cx2:.1f},{y_mid:.1f} {cx2:.1f},{y_top} Z")
        tip = _tip_esc(f"{name}\ncontrol {(cx2 - cx1) / width * 100:.1f}%\nAD {(ax2 - ax1) / width * 100:.1f}%")
        ribbons.append(f'<path d="{path}" fill="{color}" fill-opacity="0.3" data-tip="{tip}"/>')

    return f"""
<div style="display:flex; gap:10px;">
  <div style="display:flex; flex-direction:column; justify-content:space-between; width:64px; font-size:12px; font-weight:700; height:{height}px;">
    <div style="color:{palette.ACCENT_CONTROL};">Control</div>
    <div style="color:{palette.ACCENT_AD};">AD</div>
  </div>
  <svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" style="flex:1;">
    {_stacked_bar_layer(ctrl_y, bar_h, ctrl_segments, width)}
    {''.join(ribbons)}
    {_stacked_bar_layer(ad_y, bar_h, ad_segments, width)}
  </svg>
</div>"""


def stacked_bar_section(gene: str, cell_type: str | None = None, hit_enst: str | None = None) -> str:
    cell_types = data.get_gene_cell_types(gene)
    width = 1000
    blocks = []
    for ct in cell_types:
        # resolve_hit_enst always returns something, even if this OTHER
        # cell_type of the gene independently happens to carry two hits --
        # only the current (cell_type, hit_enst) this dossier is about needs
        # the exact pick; elsewhere any one representative hit is fine.
        resolved = data.resolve_hit_enst(gene, ct, hit_enst if ct == cell_type else None)
        rows = data.get_hit_rows(gene, ct, resolved).sort_values("alt_rank")
        this_hit_enst = rows.iloc[0]["hit_ENST_ID"]
        colors = palette.isoform_colors(len(rows))
        named_names = set(rows["alt_transcript_name"])
        legend_items = []
        # Pooled usage (sum raw / sum gene-total per condition) for every
        # isoform of the gene -- the same convention the DTU test itself
        # uses for chi_padj/delta_usage, so a named segment's height here
        # matches the header's percentages exactly. NOT hits_deep.csv's
        # alt_usage_pct_AD/control (an unweighted per-donor mean computed
        # separately for isoform ranking) -- the two can diverge sharply
        # when one donor's depth dominates.
        isoform_usage = data.get_gene_isoform_usage(gene, ct).set_index("transcript_name")

        ad_segments, ctrl_segments = [], []
        for color, (_, r) in zip(colors, rows.iterrows()):
            name = r["alt_transcript_name"]
            is_hit = r["alt_ENST_ID"] == this_hit_enst
            # The hit transcript gets a striking pure-red fill instead of its
            # ordinal isoform color -- it's the one segment the whole section
            # is about, and it should be unmistakable at a glance.
            fill = palette.PURE_RED if is_hit else color
            usage = isoform_usage.loc[name] if name in isoform_usage.index else None
            ad_segments.append((name, float(usage["usage_pct_AD"]) if usage is not None else 0.0, fill, is_hit))
            ctrl_segments.append((name, float(usage["usage_pct_control"]) if usage is not None else 0.0, fill, is_hit))
            weight = "700" if is_hit else "400"
            legend_items.append(
                f'<span><span class="swatch" style="background:{fill};"></span>'
                f'<span style="font-weight:{weight};">{_esc(name)}{" (hit)" if is_hit else ""}</span></span>'
            )

        # Every other isoform of the gene, not ranked/named by J1c -- drawn
        # as its own grey segment (individually hoverable) so the bar shows
        # the gene's real full isoform composition, not just the shortlist.
        other = isoform_usage.reset_index()
        other = other[~other["transcript_name"].isin(named_names)].sort_values(
            "usage_pct_AD", ascending=False
        )
        n_other = len(other)
        for _, o in other.iterrows():
            ad_segments.append((o["transcript_name"], float(o["usage_pct_AD"]), "var(--muted)", False))
            ctrl_segments.append((o["transcript_name"], float(o["usage_pct_control"]), "var(--muted)", False))
        if n_other:
            legend_items.append(
                f'<span><span class="swatch" style="background:var(--muted);"></span>'
                f'{n_other} other isoform{"s" if n_other != 1 else ""} (hover for detail)</span>'
            )

        blocks.append(f"""
<div class="cell-type-block">
  <div class="ct-label">{_esc(ct)}</div>
  {_stacked_bars_svg(ctrl_segments, ad_segments, width)}
  <div class="legend">{''.join(legend_items)}</div>
</div>""")
    return f"""
<div class="card">
  <h2>Isoform usage &middot; Control vs AD</h2>
  <div class="subtle" style="margin-top:-8px; margin-bottom:14px;">pooled usage (reads summed across donors per condition) &mdash; matches the DTU test statistic in the header</div>
  {''.join(blocks)}
</div>"""


# ---------------------------------------------------------------------------
# Section 3 — per-donor raw usage dot plot, per cell type
# ---------------------------------------------------------------------------

def _depth_legend_svg(x: float, y_top: float, bar_h: float, max_depth: float, chart_id: str) -> str:
    """One monochrome light→dark gradient swatch (fill = depth, shared
    across both conditions) plus a small key showing what the border color
    means (Control/AD identity lives on the stroke, not the fill).
    """
    bar_w = 12
    gid = f"depth-{chart_id}"
    defs = f"""<linearGradient id="{gid}" x1="0" y1="1" x2="0" y2="0">
      <stop offset="0%" stop-color="{palette.depth_gray(0.0)}"/>
      <stop offset="100%" stop-color="{palette.depth_gray(1.0)}"/>
    </linearGradient>"""
    bar = (f'<rect x="{x:.1f}" y="{y_top}" width="{bar_w}" height="{bar_h:.1f}" '
           f'fill="url(#{gid})" stroke="var(--border)" stroke-width="1"/>')
    border_key = "".join(
        f'<circle cx="{x + bar_w / 2:.1f}" cy="{y_top + bar_h + 22 + i * 16}" r="5" '
        f'fill="var(--surface-1)" stroke="{color}" stroke-width="1.5"/>'
        f'<text x="{x + bar_w + 8:.1f}" y="{y_top + bar_h + 26 + i * 16}" font-size="9" '
        f'fill="{color}">{cond}</text>'
        for i, (cond, color) in enumerate(
            (("Control", palette.ACCENT_CONTROL), ("AD", palette.ACCENT_AD))
        )
    )
    return f"""
    <text x="{x}" y="{y_top - 6}" font-size="9" fill="var(--text-secondary)">gene reads</text>
    <defs>{defs}</defs>
    {bar}
    <text x="{x + bar_w + 4}" y="{y_top + 8}" font-size="9" fill="var(--text-secondary)">{int(max_depth):,}</text>
    <text x="{x + bar_w + 4}" y="{y_top + bar_h}" font-size="9" fill="var(--text-secondary)">0</text>
    {border_key}"""


def _dot_plot_svg(df: pd.DataFrame, chart_id: str) -> str:
    pad_l, pad_t, pad_b = 42, 10, 78
    slot_w, group_gap, plot_h = 36, 46, 170
    legend_gap, legend_w = 26, 80
    colors = {"Control": palette.ACCENT_CONTROL, "AD": palette.ACCENT_AD}
    groups = [(cond, df[df["condition"] == cond].sort_values("donor")) for cond in ("Control", "AD")]
    n_total = sum(len(sub) for _, sub in groups)
    plot_w = max(n_total * slot_w + (len(groups) - 1) * group_gap, 1)
    width = pad_l + plot_w + legend_gap + legend_w
    height = pad_t + plot_h + pad_b
    y_max = max(0.05, df["psi"].max() * 1.15) if len(df) and df["psi"].notna().any() else 1.0
    max_depth = df["gene_total_count"].max() if len(df) and df["gene_total_count"].notna().any() else 0

    def y_of(psi: float) -> float:
        return pad_t + plot_h * (1 - psi / y_max)

    def depth_t(depth: float) -> float:
        return math.log1p(depth) / math.log1p(max_depth) if max_depth > 0 else 0.0

    grid = []
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        gy = pad_t + plot_h * (1 - frac)
        grid.append(f'<line x1="{pad_l}" y1="{gy:.1f}" x2="{pad_l + plot_w}" y2="{gy:.1f}" '
                     f'stroke="var(--gridline)" stroke-width="1"/>')
        grid.append(f'<text x="{pad_l - 8}" y="{gy + 4:.1f}" font-size="10" text-anchor="end">'
                     f'{frac * y_max * 100:.0f}%</text>')

    dots, donor_labels, group_labels, dividers = [], [], [], []
    x_cursor = pad_l
    for gi, (cond, sub) in enumerate(groups):
        group_start = x_cursor
        color = colors[cond]
        for r in sub.itertuples():
            x = x_cursor + slot_w / 2
            psi = r.psi if pd.notna(r.psi) else 0
            y = y_of(psi)
            fill = palette.depth_gray(depth_t(r.gene_total_count))
            tip = _tip_esc(
                f"donor: {r.donor}\ncondition: {cond}\n"
                f"raw count: {int(r.raw_count)} / gene total {int(r.gene_total_count)}\n"
                f"usage: {psi * 100:.1f}%"
            )
            dots.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="{fill}" stroke="{color}" '
                f'stroke-width="1.5" data-tip="{tip}"/>'
            )
            donor_labels.append(
                f'<text x="{x:.1f}" y="{pad_t + plot_h + 14}" font-size="9" text-anchor="end" '
                f'transform="rotate(-55 {x:.1f} {pad_t + plot_h + 14})">{_esc(r.donor)}</text>'
            )
            x_cursor += slot_w
        group_labels.append(
            f'<text x="{(group_start + x_cursor) / 2:.1f}" y="{height - 6}" font-size="12" '
            f'text-anchor="middle" fill="{color}" font-weight="700">{cond}</text>'
        )
        if gi < len(groups) - 1:
            dx = x_cursor + group_gap / 2
            dividers.append(f'<line x1="{dx:.1f}" y1="{pad_t}" x2="{dx:.1f}" y2="{pad_t + plot_h}" '
                             f'stroke="var(--gridline)" stroke-width="1" stroke-dasharray="2,2"/>')
        x_cursor += group_gap

    legend = _depth_legend_svg(pad_l + plot_w + legend_gap, pad_t, plot_h, max_depth, chart_id)

    return f"""<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">
    {''.join(grid)}
    {''.join(dividers)}
    <line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t + plot_h}" stroke="var(--baseline)"/>
    {''.join(dots)}
    {''.join(donor_labels)}
    {''.join(group_labels)}
    {legend}
  </svg>"""


def donor_dot_section(gene: str, cell_type: str | None = None, hit_enst: str | None = None) -> str:
    cell_types = data.get_gene_cell_types(gene)
    blocks = []
    for ct in cell_types:
        resolved = data.resolve_hit_enst(gene, ct, hit_enst if ct == cell_type else None)
        hit_row = data.get_hit_rows(gene, ct, resolved).iloc[0]
        transcript_name = hit_row["hit_transcript_name"]
        df = data.get_donor_usage(gene, transcript_name, ct)
        blocks.append(f"""
<div class="cell-type-block">
  <div class="ct-label">{_esc(ct)} &middot; <span class="muted">{_esc(transcript_name)}</span></div>
  {_dot_plot_svg(df, ct)}
</div>""")
    return f"""
<div class="card">
  <h2>Per-donor raw usage</h2>
  <div class="subtle" style="margin-bottom:14px;">dot darkness = gene-level read depth (log scale) &middot; hover for raw counts</div>
  {''.join(blocks)}
</div>"""


# ---------------------------------------------------------------------------
# Section 4 — canonical vs alt sequence comparison
# ---------------------------------------------------------------------------

def _ruler_svg(y: float, seq_len: int, px_per_aa: float) -> str:
    """Amino-acid position ticks along the canonical bar's reference frame."""
    if seq_len <= 0:
        return ""
    step = 500
    for candidate in (20, 50, 100, 200, 500):
        if seq_len / candidate <= 12:
            step = candidate
            break
    parts = []
    pos = 0
    while pos <= seq_len:
        x = pos * px_per_aa
        parts.append(f'<line x1="{x:.1f}" y1="{y}" x2="{x:.1f}" y2="{y + 5}" '
                     f'stroke="var(--muted)" stroke-width="1"/>')
        if pos > 0:
            parts.append(f'<text x="{x:.1f}" y="{y - 3}" font-size="12" text-anchor="middle" '
                         f'fill="var(--muted)">{pos}</text>')
        pos += step
    return "".join(parts)


def _domain_track_svg(y: float, domains: list[dict], px_per_aa: float,
                       color_map: dict[str, str], track_h: float) -> str:
    """A caliper bracket -- end ticks + a spanning line -- under a centered
    name label, instead of a filled box. Reads as an annotation pointing at
    a range, not another data-bearing bar competing with the protein bar
    itself for visual weight.
    """
    tick_h = 8
    line_y = y + track_h - tick_h / 2 - 1
    parts = []
    for d in domains:
        x = d["start"] * px_per_aa
        x2 = (d["end"] + 1) * px_per_aa
        color = color_map.get(d["name"], "var(--muted)")
        tip = _tip_esc(f"{d['name']}\n{d['start']}-{d['end']}")
        parts.append(
            f'<rect x="{x:.1f}" y="{y}" width="{max(x2 - x, 1):.1f}" height="{track_h}" '
            f'fill="transparent" data-tip="{tip}"/>'
            f'<line x1="{x:.1f}" y1="{line_y - tick_h / 2:.1f}" x2="{x:.1f}" y2="{line_y + tick_h / 2:.1f}" '
            f'stroke="{color}" stroke-width="1.5"/>'
            f'<line x1="{x2:.1f}" y1="{line_y - tick_h / 2:.1f}" x2="{x2:.1f}" y2="{line_y + tick_h / 2:.1f}" '
            f'stroke="{color}" stroke-width="1.5"/>'
            f'<line x1="{x:.1f}" y1="{line_y:.1f}" x2="{x2:.1f}" y2="{line_y:.1f}" '
            f'stroke="{color}" stroke-width="1.5"/>'
            f'<text x="{(x + x2) / 2:.1f}" y="{y + 11:.1f}" font-size="12" text-anchor="middle" '
            f'fill="var(--text-primary)">{_esc(d["name"])}</text>'
        )
    return "".join(parts)


def _protein_bar_svg(blocks: list[tuple[int, int, str]], px_per_aa: float,
                      y: float, bar_h: float) -> str:
    parts = []
    for start, end, op in blocks:
        x = start * px_per_aa
        w = max((end - start + 1) * px_per_aa, 0.6)
        color = "var(--bar-black)" if op == "match" else palette.ACCENT_AD
        parts.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{color}"/>')
    return "".join(parts)


def _ribbons_svg(segments: list[dict], y_top: float, y_bot: float, px_per_aa: float) -> str:
    """Curved bands connecting each aligned region between the canonical and
    alt bars -- grey for a match, red for a mismatch. An indel has no
    counterpart on the other side, so it gets a small triangular flag
    pointing into the ribbon gap instead of a connecting band.
    """
    y_mid = (y_top + y_bot) / 2
    parts = []
    for s in segments:
        op = s["op"]
        if op in ("match", "mismatch"):
            cx1, cx2 = s["canon_start"] * px_per_aa, (s["canon_end"] + 1) * px_per_aa
            ax1, ax2 = s["alt_start"] * px_per_aa, (s["alt_end"] + 1) * px_per_aa
            fill = "var(--match-gray)" if op == "match" else palette.ACCENT_AD
            opacity = 0.6 if op == "match" else 0.5
            path = (f"M {cx1:.1f},{y_top} C {cx1:.1f},{y_mid:.1f} {ax1:.1f},{y_mid:.1f} {ax1:.1f},{y_bot} "
                    f"L {ax2:.1f},{y_bot} C {ax2:.1f},{y_mid:.1f} {cx2:.1f},{y_mid:.1f} {cx2:.1f},{y_top} Z")
            tip = _tip_esc(f"{op}\ncanonical {s['canon_start']}-{s['canon_end']}\n"
                            f"alt {s['alt_start']}-{s['alt_end']}")
            parts.append(f'<path d="{path}" fill="{fill}" fill-opacity="{opacity}" data-tip="{tip}"/>')
        elif op == "deleted":
            x1, x2 = s["canon_start"] * px_per_aa, (s["canon_end"] + 1) * px_per_aa
            xm = (x1 + x2) / 2
            n = s["canon_end"] - s["canon_start"] + 1
            tip = _tip_esc(f"deleted in alt\ncanonical {s['canon_start']}-{s['canon_end']} ({n} aa)")
            parts.append(f'<polygon points="{x1:.1f},{y_top} {x2:.1f},{y_top} {xm:.1f},{y_top + 14:.1f}" '
                         f'fill="{palette.ACCENT_AD}" fill-opacity="0.65" data-tip="{tip}"/>')
        elif op == "inserted":
            x1, x2 = s["alt_start"] * px_per_aa, (s["alt_end"] + 1) * px_per_aa
            xm = (x1 + x2) / 2
            n = s["alt_end"] - s["alt_start"] + 1
            tip = _tip_esc(f"inserted in alt\nalt {s['alt_start']}-{s['alt_end']} ({n} aa)")
            parts.append(f'<polygon points="{x1:.1f},{y_bot} {x2:.1f},{y_bot} {xm:.1f},{y_bot - 14:.1f}" '
                         f'fill="{palette.ACCENT_AD}" fill-opacity="0.65" data-tip="{tip}"/>')
    return "".join(parts)


def _sequence_diff_svg(canonical_seq: str, alt_seq: str,
                        canonical_domains: list[dict], alt_domains: list[dict]) -> str:
    x0, width = 74, 1340
    plot_w = width - x0 - 12
    segments = sequence_diff.align_segments(canonical_seq, alt_seq)
    canon_len, alt_len = len(canonical_seq), len(alt_seq)
    px_per_aa = plot_w / max(canon_len, alt_len, 1)

    ruler_h, track_h, bar_h, ribbon_h, gap = 18, 28, 22, 56, 5
    ruler_y = ruler_h
    canon_track_y = ruler_y + gap
    canon_bar_y = canon_track_y + track_h + gap
    ribbon_top = canon_bar_y + bar_h
    ribbon_bot = ribbon_top + ribbon_h
    alt_bar_y = ribbon_bot
    alt_track_y = alt_bar_y + bar_h + gap
    height = alt_track_y + track_h + 6

    color_map = palette.domain_colors(
        [d["name"] for d in canonical_domains] + [d["name"] for d in alt_domains]
    )
    canon_blocks = [(s["canon_start"], s["canon_end"], s["op"]) for s in segments if s["canon_start"] is not None]
    alt_blocks = [(s["alt_start"], s["alt_end"], s["op"]) for s in segments if s["alt_start"] is not None]

    inner = _ruler_svg(ruler_y, canon_len, px_per_aa)
    inner += _domain_track_svg(canon_track_y, canonical_domains, px_per_aa, color_map, track_h)
    inner += _protein_bar_svg(canon_blocks, px_per_aa, canon_bar_y, bar_h)
    inner += _ribbons_svg(segments, ribbon_top, ribbon_bot, px_per_aa)
    inner += _protein_bar_svg(alt_blocks, px_per_aa, alt_bar_y, bar_h)
    inner += _domain_track_svg(alt_track_y, alt_domains, px_per_aa, color_map, track_h)

    labels = (
        f'<text x="{x0 - 8}" y="{canon_bar_y + bar_h / 2 + 4:.1f}" font-size="11" text-anchor="end">canonical</text>'
        f'<text x="{x0 - 8}" y="{canon_bar_y + bar_h / 2 + 16:.1f}" font-size="9" text-anchor="end" '
        f'fill="var(--muted)">{canon_len} aa</text>'
        f'<text x="{x0 - 8}" y="{alt_bar_y + bar_h / 2 + 4:.1f}" font-size="11" text-anchor="end">alt</text>'
        f'<text x="{x0 - 8}" y="{alt_bar_y + bar_h / 2 + 16:.1f}" font-size="9" text-anchor="end" '
        f'fill="var(--muted)">{alt_len} aa</text>'
    )

    return (f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
            f'{labels}<g transform="translate({x0},0)">{inner}</g></svg>')


def _usage_text(usage_control: float, usage_ad: float) -> str:
    """Control usage -> AD usage as colored numbers with an arrow between --
    replaces the earlier mini-dumbbell chart; same two colors, no SVG."""
    return (
        f'<span style="color:{palette.ACCENT_CONTROL}; font-weight:700;">{usage_control * 100:.1f}%</span>'
        f'<span class="muted"> &rarr; </span>'
        f'<span style="color:{palette.ACCENT_AD}; font-weight:700;">{usage_ad * 100:.1f}%</span>'
    )


def _domain_list(enst_id: str, kind: str) -> list[dict]:
    return data.domains_for(enst_id, kind)


def _change_badge_class(length_diff: float) -> str:
    """Red = protein got shorter, blue = protein got longer, black = same
    length (substitution, identical, or a net-zero internal swap)."""
    if length_diff < 0:
        return "change-decrease"
    if length_diff > 0:
        return "change-increase"
    return "change-same"


def _one_sequence_block(row: pd.Series, usage_ad: float, usage_control: float) -> str:
    change = row["protein_change_type"]
    lost = row.get("domains_lost")
    gained = row.get("domains_gained")
    badge_class = _change_badge_class(row.get("protein_length_diff", 0) or 0)
    detail = f'<span class="badge {badge_class}">{_esc(change)}</span>'
    if pd.notna(lost) and str(lost).strip():
        detail += f'<span class="subtle">&nbsp;lost: <strong>{_esc(lost)}</strong></span>'
    if pd.notna(gained) and str(gained).strip():
        detail += f'<span class="subtle">&nbsp;gained: <strong>{_esc(gained)}</strong></span>'

    canonical_seq = row.get("canonical_protein_seq") or ""
    alt_seq = row.get("alt_protein_seq") or ""
    canonical_domains = _domain_list(row["canonical_enst"], "canonical")
    alt_domains = _domain_list(row["alt_ENST_ID"], "alt")

    svg = (_sequence_diff_svg(canonical_seq, alt_seq, canonical_domains, alt_domains)
           if canonical_seq and alt_seq else '<div class="subtle">sequence unavailable</div>')

    return f"""
<div style="padding:14px 0; border-bottom:1px solid var(--gridline);">
  <div style="display:flex; align-items:center; gap:10px; margin-bottom:8px; flex-wrap:wrap;">
    <strong>{_esc(row['alt_transcript_name'])}</strong>
    <span class="muted">rank {int(row['alt_rank'])}</span>
    {detail}
    {_usage_text(usage_control, usage_ad)}
  </div>
  {svg}
</div>"""


def sequence_section(rows: pd.DataFrame) -> str:
    gene, cell_type = rows.iloc[0]["gene_name"], rows.iloc[0]["cell_type"]
    isoform_usage = data.get_gene_isoform_usage(gene, cell_type).set_index("transcript_name")

    def _usage_for(row: pd.Series) -> tuple[float, float]:
        name = row["alt_transcript_name"]
        if name in isoform_usage.index:
            u = isoform_usage.loc[name]
            return float(u["usage_pct_AD"]), float(u["usage_pct_control"])
        return 0.0, 0.0

    group = rows.iloc[0]["master_group"]
    if group == "trial_failure_candidate":
        alt_rows = rows[rows["is_canonical"] != True].sort_values("alt_rank")  # noqa: E712
        title = "Ranked alternate isoforms vs canonical"
        body = "".join(
            _one_sequence_block(r, *_usage_for(r)) for _, r in alt_rows.iterrows()
        )
        if not len(alt_rows):
            body = '<div class="subtle">no ranked alternates with a protein sequence.</div>'
    else:
        title = "Alternative vs canonical"
        row = rows.iloc[0]
        body = _one_sequence_block(row, *_usage_for(row))
    return f"""
<div class="card">
  <h2>{title}</h2>
  <div class="subtle" style="margin-bottom:4px;">ribbons connect aligned regions (grey = match, red = changed); flags mark residues with no counterpart &middot; domain brackets are stable-colored per Pfam name &middot; badge color: red = shorter, blue = longer, black = same length &middot; Control&rarr;AD usage in colored text</div>
  {body}
</div>"""


# ---------------------------------------------------------------------------
# Section 5 — drug evidence table
# ---------------------------------------------------------------------------

def _phase_meter(phase) -> str:
    """4 filled/unfilled squares -- clinical phase as a glanceable meter
    instead of a bare integer."""
    try:
        phase = int(phase) if pd.notna(phase) else 0
    except (TypeError, ValueError):
        phase = 0
    dots = "".join(
        f'<span style="display:inline-block; width:9px; height:9px; border-radius: 0; '
        f'margin-right:3px; background:{"var(--text-primary)" if i < phase else "var(--gridline)"};"></span>'
        for i in range(4)
    )
    return f'{dots}<span class="muted" style="font-size:11px; margin-left:4px;">phase {phase}/4</span>'


def _evidence_indicator(has_evidence: bool) -> str:
    """Status color is reserved and separate from the Control/AD accent
    pair -- green/red here means good/bad evidence, not a series."""
    if has_evidence:
        bg, symbol, label = "var(--status-good)", "&#10003;", "drug evidence"
    else:
        bg, symbol, label = "var(--status-critical)", "&#10007;", "no drug evidence"
    return (
        f'<div style="margin-top:10px; padding-top:10px; border-top:1px solid var(--border); '
        f'display:flex; align-items:center; gap:7px;">'
        f'<span style="display:inline-flex; align-items:center; justify-content:center; '
        f'width:18px; height:18px; border-radius: 0; background:{bg}; color:#ffffff; '
        f'font-size:12px; font-weight:700; flex-shrink:0;">{symbol}</span>'
        f'<span class="subtle" style="font-size:11px;">{label}</span></div>'
    )


def _drug_card(title: str, body: str, has_evidence: bool) -> str:
    return f"""
<div style="background:var(--page-plane); border:1px solid var(--border); border-radius: 0; padding:14px 16px;">
  <div style="font-size:11px; text-transform:uppercase; letter-spacing:0.06em; color:var(--text-secondary); font-weight:700; margin-bottom:8px;">{_esc(title)}</div>
  {body}
  {_evidence_indicator(has_evidence)}
</div>"""


def _has_text(val) -> bool:
    return pd.notna(val) and str(val).strip() != ""


def _drug_field(row: pd.Series, label: str, col: str) -> str:
    val = row.get(col)
    if pd.isna(val) or val == "" or val == 0:
        return ""
    return f'<div style="font-size:12px; margin-top:4px;"><span class="muted">{_esc(label)}:</span> {_esc(val)}</div>'


def _trial_stop_note(row: pd.Series) -> str:
    """Terminated/withdrawn/suspended clinical trial evidence -- a drug was
    actually tried against this target and stopped, which max-phase alone
    doesn't surface. Status-critical color, separate from the has-evidence
    check mark below it."""
    n_stopped = row.get("ot_trials_terminated", 0)
    if pd.isna(n_stopped) or n_stopped <= 0:
        return ""
    n_total = row.get("ot_trials_total", 0)
    reasons = str(row.get("ot_trial_stop_reasons") or "").strip()
    example = str(row.get("ot_trial_stop_example") or "").strip()

    reasons_html = ""
    if reasons:
        tags = "".join(
            f'<span class="badge" style="background:var(--status-critical); color:#fff; '
            f'margin-right:4px; margin-top:4px;">{_esc(r.replace("_", " "))}</span>'
            for r in reasons.split("|") if r
        )
        reasons_html = f'<div style="margin-top:6px;">{tags}</div>'

    example_html = (
        f'<div class="subtle" style="font-size:11px; margin-top:6px; font-style:italic;">'
        f'&ldquo;{_esc(example)}&rdquo;</div>'
        if example else ""
    )

    return (
        f'<div style="margin-top:8px; padding:8px 10px; border-left:3px solid var(--status-critical); '
        f'background:var(--surface);">'
        f'<div style="font-size:12px; font-weight:700; color:var(--status-critical);">'
        f'{int(n_stopped)}/{int(n_total)} trial{"s" if n_total != 1 else ""} stopped early</div>'
        f'{reasons_html}{example_html}</div>'
    )


def drug_table_section(row: pd.Series) -> str:
    cards = []

    chembl = ""
    if row.get("chembl_max_phase", 0) > 0:
        chembl += f'<div style="margin-bottom:6px;">{_phase_meter(row.get("chembl_max_phase"))}</div>'
    chembl += _drug_field(row, "Target", "chembl_target_id")
    chembl += _drug_field(row, "Drugs", "drug_names")
    chembl += _drug_field(row, "Drug count", "n_drugs")
    if row.get("is_druggable"):
        chembl += '<div style="margin-top:6px;"><span class="badge control">druggable</span></div>'
    chembl += _drug_field(row, "Bioactive compounds", "chembl_bioactive_compounds")
    chembl += _drug_field(row, "Best pChEMBL", "chembl_best_pchembl")
    if chembl:
        chembl_has_drug = _has_text(row.get("drug_names")) or row.get("chembl_max_phase", 0) >= 1
        cards.append(_drug_card("ChEMBL", chembl, chembl_has_drug))

    ot = ""
    if row.get("ot_max_phase", 0) > 0:
        ot += f'<div style="margin-bottom:6px;">{_phase_meter(row.get("ot_max_phase"))}</div>'
    ot += _drug_field(row, "Drugs", "ot_drug_names")
    ot += _drug_field(row, "Drug count", "ot_n_drugs")
    ot += _trial_stop_note(row)
    if ot:
        ot_has_drug = _has_text(row.get("ot_drug_names")) or row.get("ot_max_phase", 0) >= 1
        cards.append(_drug_card("Open Targets", ot, ot_has_drug))

    pharos = ""
    tdl = row.get("pharos_tdl")
    if pd.notna(tdl) and str(tdl).strip():
        pharos += f'<div style="margin-bottom:6px;"><span class="badge change-same">{_esc(tdl)}</span></div>'
    pharos += _drug_field(row, "Ligands", "pharos_n_ligands")
    pharos += _drug_field(row, "Drugs", "pharos_n_drugs")
    if pharos:
        pharos_has_drug = row.get("pharos_n_drugs", 0) > 0
        cards.append(_drug_card("Pharos", pharos, pharos_has_drug))

    dgidb_n = row.get("dgidb_interactions")
    if pd.notna(dgidb_n) and dgidb_n > 0:
        dgidb_n = int(dgidb_n)
        dgidb = (f'<div style="font-size:24px; font-weight:700;">{dgidb_n}</div>'
                 f'<div class="subtle">interaction{"s" if dgidb_n != 1 else ""}</div>')
        cards.append(_drug_card("DGIdb", dgidb, True))

    body = (
        f'<div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:14px;">'
        f'{"".join(cards)}</div>'
        if cards else '<div class="subtle">no drug-database evidence recorded</div>'
    )
    return f"""
<div class="card">
  <h2>Drug target database findings</h2>
  {body}
</div>"""


# ---------------------------------------------------------------------------
# Top-level assembly
# ---------------------------------------------------------------------------

def render(gene: str, cell_type: str, hit_enst: str | None = None) -> str:
    rows = data.get_hit_rows(gene, cell_type, hit_enst)
    hit_row = rows.iloc[0]

    body = (
        header_section(hit_row)
        + stacked_bar_section(gene, cell_type, hit_row["hit_ENST_ID"])
        + donor_dot_section(gene, cell_type, hit_row["hit_ENST_ID"])
        + sequence_section(rows)
        + drug_table_section(hit_row)
    )

    return f"""<!doctype html>
<html><head><meta charset="utf-8"/>
<title>{_esc(gene)} &middot; {_esc(cell_type)} dossier</title>
<style>{_font_face_css()}{_BASE_CSS}</style>
</head><body>
{THEME_TOGGLE_HTML}
<div class="viz-root">{body}</div>
<div class="tip"></div>
<script>{THEME_TOGGLE_JS}</script>
<script>{_TIP_JS}</script>
</body></html>"""
