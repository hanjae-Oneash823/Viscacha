reinitialize
set bg_rgb, [1.0, 1.0, 1.0]
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 0
set ray_shadows, 0
set ambient, 0.48
set direct, 0.62
set specular, 0.20
set shininess, 22
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set depth_cue, 0

load outputs/docking_campaign/systems/BACE1_verubecestat/inputs/5HU1.pdb, canonical_complex
load outputs/docking_campaign/analysis/bc_candidates/BACE1_202_401/sp_P56817-5_BACE1_HUMAN_Isoform_5_of_Beta-secretase_1_OS_Homo_sapiens_OX_9606_GN_BACE1_aligned_rank_001_alphafold2_ptm_model_1_seed_003.pdb, alternate
load outputs/docking_campaign/systems/BACE1_verubecestat_401/runs/alternate_401_seed1103_ex32/docked_poses.pdbqt, alt_docked
create alt_pose1, alt_docked, 1, 1
remove solvent

select canonical, canonical_complex and chain A and polymer.protein
select crystal_drug, canonical_complex and chain A and resn 66F
select alternate_core, alternate and resi 21-334
select missing_contacts, canonical and resi 70-75+91-96
select retained_can, canonical and resi 132+133+134+135+137+168+169+171+176+179+289+290+291+292+293+396
select retained_alt, alternate and resi 32+33+34+35+37+68+69+71+76+79+189+190+191+192+193+296
select catalytic_can, canonical and resi 93+289
select catalytic_alt, alternate and resi 189

hide everything
show cartoon, canonical or alternate_core
color 0x9AA5B1, canonical
color 0x167D8D, alternate
color 0xE15554, missing_contacts
set cartoon_transparency, 0.52, canonical
set cartoon_transparency, 0.08, alternate_core
show sticks, crystal_drug
color 0xF2B134, crystal_drug
set stick_radius, 0.28, crystal_drug
show spheres, missing_contacts and name CA
set sphere_scale, 0.46, missing_contacts and name CA
color 0xE15554, missing_contacts
orient canonical or alternate_core
zoom canonical or alternate_core, 5
turn y, 10
turn x, -8
ray 2200, 1600
png outputs/docking_campaign/figures/bc_candidates/B_BACE1_isoform_overlay.png, dpi=300

hide everything
show cartoon, (canonical or alternate) within 12 of crystal_drug
set cartoon_transparency, 0.72, canonical
set cartoon_transparency, 0.48, alternate
color 0x6D7A88, canonical
color 0x167D8D, alternate
show sticks, crystal_drug or alt_pose1
color 0xF2B134, crystal_drug
color 0x8E44AD, alt_pose1
set stick_radius, 0.25
show sticks, retained_can or retained_alt or missing_contacts or catalytic_can or catalytic_alt
color 0xDCE2E8, retained_can
color 0x5FC3B5, retained_alt
color 0xE15554, missing_contacts
color 0xC62828, canonical and resi 93
color 0x3949AB, canonical and resi 289
color 0x3949AB, alternate and resi 189
show spheres, (canonical and resi 93+289 and name CA) or (alternate and resi 189 and name CA)
set sphere_scale, 0.38
zoom crystal_drug or alt_pose1, 7
ray 2200, 1600
png outputs/docking_campaign/figures/bc_candidates/B_BACE1_pose_displacement.png, dpi=300
quit
