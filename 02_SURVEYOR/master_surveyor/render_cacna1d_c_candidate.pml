reinitialize
set bg_rgb, [1.0, 1.0, 1.0]
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set ambient, 0.42
set direct, 0.68
set specular, 0.22
set shininess, 24
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set depth_cue, 0

load outputs/master_surveyor/cache/structures/d96a034227d14cdd/canonical_afdb.pdb, canonical
load outputs/master_surveyor/cache/structures/a90d14f0a0086f3f/alt_colabfold/alt_a90d14f0a0086f3f_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_003.pdb, alternate
remove solvent
align alternate and name CA, canonical and name CA, cycles=5, cutoff=2.0

select can_pocket, canonical and resi 1078+1081+1082+1085+1154+1156+1194+1198+1205+1209+1212+1489+1492+1493
select alt_pocket, alternate and resi 1098+1101+1102+1105+1174+1176+1214+1218+1225+1229+1232+1509+1512+1513
select lost_tail, canonical and resi 1606-2161

hide everything
show cartoon, canonical or alternate
color gray80, canonical
color 0x177E89, alternate
color 0xE45756, lost_tail
set cartoon_transparency, 0.62, canonical and not lost_tail
set cartoon_transparency, 0.10, lost_tail
set cartoon_transparency, 0.05, alternate
show spheres, (can_pocket or alt_pocket) and name CA
set sphere_scale, 0.62, (can_pocket or alt_pocket) and name CA
color 0x6C2BD9, can_pocket and name CA
color 0xF2B134, alt_pocket and name CA
orient canonical
zoom canonical, 3
turn x, -7
turn y, 11
ray 2600, 1900
png outputs/docking_campaign/figures/bc_candidates/C_CACNA1D_full_structure_overlay.png, dpi=300

hide everything
show cartoon, (canonical or alternate) within 13 of (can_pocket or alt_pocket)
set cartoon_transparency, 0.52, canonical
set cartoon_transparency, 0.28, alternate
show sticks, can_pocket or alt_pocket
set stick_radius, 0.20
color 0x355C7D, canonical and not can_pocket
color 0x26A69A, alternate and not alt_pocket
color 0x6C2BD9, can_pocket
color 0xF2B134, alt_pocket
show spheres, (can_pocket or alt_pocket) and name CA
set sphere_scale, 0.35, (can_pocket or alt_pocket) and name CA
zoom can_pocket or alt_pocket, 8
ray 2600, 1900
png outputs/docking_campaign/figures/bc_candidates/C_CACNA1D_retained_pocket_overlay.png, dpi=300
quit
