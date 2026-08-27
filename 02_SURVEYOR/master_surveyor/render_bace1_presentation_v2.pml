# Improved presentation renders for BACE1 canonical vs BACE1-202 (401 aa)
# Two figures: (1) domain-loss overview, (2) active-site pocket close-up

reinitialize
set bg_rgb, [1.0, 1.0, 1.0]
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 1
set ray_trace_color, [0.25, 0.28, 0.32]
set ray_trace_gain, 0.15
set ray_shadows, 0
set ambient, 0.55
set direct, 0.55
set specular, 0.15
set shininess, 15
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_side_chain_helper, 1
set depth_cue, 1
set fog_start, 0.45
set ray_opaque_background, 1
set label_size, 24
set label_color, black
set label_font_id, 7
set label_outline_color, white
set dash_gap, 0.3
set dash_width, 3

load outputs/docking_campaign/systems/BACE1_verubecestat/inputs/5HU1.pdb, canonical_complex
load outputs/docking_campaign/analysis/bc_candidates/BACE1_202_401/sp_P56817-5_BACE1_HUMAN_Isoform_5_of_Beta-secretase_1_OS_Homo_sapiens_OX_9606_GN_BACE1_aligned_rank_001_alphafold2_ptm_model_1_seed_003.pdb, alternate
load outputs/docking_campaign/systems/BACE1_verubecestat_401/runs/alternate_401_seed1103_ex32/docked_poses.pdbqt, alt_docked
create alt_pose1, alt_docked, 1, 1
remove solvent
remove chain B

select canonical, canonical_complex and chain A and polymer.protein
select crystal_drug, canonical_complex and chain A and resn 66F
select alternate_core, alternate and resi 21-334
select deleted_region, canonical and resi 57-120
select missing_contacts, canonical and resi 70-75+91-96
select catalytic_lost, canonical and resi 93
select catalytic_retained, canonical and resi 289
select catalytic_retained_alt, alternate and resi 189

# ============ FIGURE 1: domain-loss overview ============
hide everything
bg_color white

show cartoon, canonical or alternate_core
color 0x8D96A1, canonical
set cartoon_transparency, 0.22, canonical
color 0x0F9B8E, alternate_core
set cartoon_transparency, 0.0, alternate_core

show cartoon, deleted_region
color 0xE0483E, deleted_region
set cartoon_transparency, 0.0, deleted_region

show sticks, crystal_drug
util.cbay crystal_drug
color 0xF4B942, crystal_drug and elem C
set stick_radius, 0.24, crystal_drug

show spheres, catalytic_lost and name CA
color 0xB0202A, catalytic_lost
set sphere_scale, 0.85, catalytic_lost

show spheres, catalytic_retained and name CA
color 0x1B4F91, catalytic_retained
set sphere_scale, 0.85, catalytic_retained

orient canonical or alternate_core
turn y, 25
turn x, -10
zoom canonical or alternate_core, 3

ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/bace1_overview_raw.png, dpi=300

# ============ FIGURE 2: active-site close-up ============
hide everything
bg_color white

show cartoon, byres ((canonical or alternate) within 9 of crystal_drug)
color 0xC7CDD4, canonical
color 0x7FD4C6, alternate
set cartoon_transparency, 0.7, canonical
set cartoon_transparency, 0.7, alternate

show sticks, crystal_drug
color 0xF4B942, crystal_drug and elem C
set stick_radius, 0.28, crystal_drug

show sticks, alt_pose1
color 0x8E44AD, alt_pose1 and elem C
set stick_radius, 0.28, alt_pose1

show sticks, catalytic_lost and not name C+N+O
color 0xB0202A, catalytic_lost and elem C
set stick_radius, 0.3, catalytic_lost
show nb_spheres, catalytic_lost and not name C+N+O
set nb_spheres_size, 0.22, catalytic_lost

pseudoatom crystal_com, selection=crystal_drug
pseudoatom alt_com, selection=alt_pose1
distance disp_dist, crystal_com, alt_com
hide labels, disp_dist
color 0x444444, disp_dist
set dash_radius, 0.05, disp_dist
hide everything, crystal_com or alt_com

orient crystal_drug or alt_pose1 or (catalytic_lost and not name C+N+O)
zoom crystal_drug or alt_pose1 or (catalytic_lost and not name C+N+O), 5
turn y, 15

ray 2400, 1800
png /tmp/claude-1813/-home-welcome3-Viscacha-pipeline/7f2faea9-593c-4a2e-ae4d-ecc12e47472b/scratchpad/bace1_closeup_raw.png, dpi=300

quit
