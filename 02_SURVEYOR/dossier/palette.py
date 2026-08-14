"""DOSSIER — shared color roles. Monochrome base + two accents (Control/AD),
per the dataviz skill's palette method. Values match the skill's validated
reference palette (references/palette.md): blue/red is its diverging pair,
reused here as the Control/AD identity channel used consistently across
every chart in the report.
"""

# ---- monochrome ink / chrome (light mode; dark overrides live in the CSS) --
SURFACE       = "#fcfcfb"
PAGE_PLANE    = "#f9f9f7"
TEXT_PRIMARY  = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
MUTED         = "#898781"
GRIDLINE      = "#e1e0d9"
BASELINE      = "#c3c2b7"
BORDER        = "rgba(11,11,11,0.10)"

# ---- the two accents -------------------------------------------------------
# Raw hex (light-mode values) -- used only to seed the --accent-control /
# --accent-ad CSS custom properties themselves. Everywhere else in the
# codebase, use ACCENT_CONTROL / ACCENT_AD below (the var() reference), so a
# color baked into an SVG fill at generation time still repaints correctly
# when the viewer toggles to the dark/neon theme at runtime.
ACCENT_CONTROL_HEX = "#2a78d6"   # blue
ACCENT_AD_HEX      = "#e34948"   # red
PURE_RED_HEX       = "#e60000"

ACCENT_CONTROL = "var(--accent-control)"
ACCENT_AD      = "var(--accent-ad)"      # doubles as the "changed" color in
                                           # the sequence-diff bars (mismatch/indel/insertion)
# A fully saturated red, distinct from (and stronger than) ACCENT_AD --
# reserved for one-off "this specific mark is the hit" call-outs where the
# softer accent red would blend in with the rest of the AD-colored chrome.
PURE_RED = "var(--pure-red)"

# ---- "matched" gray for the sequence-diff bars/ribbons ---------------------
# Deliberately darker/heavier than the chrome grays (gridline/baseline) --
# in this one diagram gray IS a data color (the matched-region fill), not
# furniture, so it needs to read as present, not recede like a gridline.
MATCH_GRAY = "#6a6964"

# ---- categorical identity for ranked-isoform stacked-bar segments --------
# Fixed order (never cycled), matching the dataviz skill's 8-hue theme.
# Blue and red are held back -- they're the Control/AD accent, and reusing
# them at full saturation for isoform identity would collide with that
# meaning. Toned down (blended toward neutral gray) so the segments read as
# a quiet categorical set, not a second pair of accents.
_CATEGORICAL = [
    "#2a78d6",  # blue
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
    "#e87ba4",  # magenta
    "#eb6834",  # orange
]
_TONE_GRAY = (0x8a, 0x89, 0x84)
_TONE_AMOUNT = 0.42  # fraction blended toward _TONE_GRAY


def _blend(hex_color: str, target: tuple[int, int, int], amount: float) -> str:
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    rgb = (
        round(r + amount * (target[0] - r)),
        round(g + amount * (target[1] - g)),
        round(b + amount * (target[2] - b)),
    )
    return "#%02x%02x%02x" % rgb


def isoform_colors(n: int) -> list[str]:
    """Toned-down categorical colors, fixed order, one per ranked isoform.
    Cycles only past 8 (the fixed theme's size), which should not happen in
    practice given J1c's per-gene isoform cap.
    """
    return [_blend(_CATEGORICAL[i % len(_CATEGORICAL)], _TONE_GRAY, _TONE_AMOUNT) for i in range(n)]


def domain_colors(names: list[str]) -> dict[str, str]:
    """Stable domain-name -> toned-down categorical color, first-seen order,
    so the same Pfam domain reads as the same color on the canonical track,
    the alt track, and across every ranked alt of the same gene.
    """
    unique = list(dict.fromkeys(names))
    return dict(zip(unique, isoform_colors(len(unique))))


# Monochrome depth ramp -- pure grayscale (no hue), so fill encodes only
# magnitude (gene-level read depth) and needs one shared legend, not one
# per condition. Condition identity lives on the dot's border instead
# (ACCENT_CONTROL/ACCENT_AD), keeping the two channels visually separate.
DEPTH_GRAY_LIGHT = (0xe6, 0xe5, 0xe0)
DEPTH_GRAY_DARK  = (0x2a, 0x2a, 0x28)


def depth_gray(t: float) -> str:
    """t in [0, 1]: 0 = lightest (shallow), 1 = darkest (deep)."""
    t = max(0.0, min(1.0, t))
    rgb = tuple(
        round(DEPTH_GRAY_LIGHT[c] + t * (DEPTH_GRAY_DARK[c] - DEPTH_GRAY_LIGHT[c]))
        for c in range(3)
    )
    return "#%02x%02x%02x" % rgb


# ---- dark theme: a deliberate, distinct world, not an inverted tint --------
# Defined once here so render.py (dossiers) and generate_index.py (homepage)
# can't drift into two different "dark modes." Pure black page, near-black
# card surfaces (just enough lift for edges to read without a border doing
# all the work), pure white text, saturated neon accents. --glow holds a
# text-shadow value (off in light mode via the base .viz-root's `--glow: none`)
# so `text-shadow: var(--glow);` glows only in dark mode, no duplicate rules.
NEON_DARK_CSS_VARS = (
    "--surface-1: #0a0a0f; --page-plane: #000000; "
    "--text-primary: #ffffff; --text-secondary: #b0b0c0; "
    "--muted: #6e6e80; --gridline: #1c1c26; --baseline: #2e2e3d; "
    "--border: rgba(255,255,255,0.14); --match-gray: #8f8fa3; --bar-black: #ffffff; "
    "--accent-control: #00e5ff; --accent-ad: #ff2d6b; --pure-red: #ff0044; "
    "--status-good: #39ff88; --status-critical: #ff2d6b; "
    "--glow: 0 0 10px currentColor;"
)

# Index-page-only identity colors (trial_failure / drug_repurposing), layered
# on top of NEON_DARK_CSS_VARS -- neon violet/green instead of the light
# theme's toned violet/aqua.
NEON_DARK_TF = "#b537f2"
NEON_DARK_DR = "#39ff88"
