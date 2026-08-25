# BIOVIA Discovery Studio candidate pairs

This is the complete set of candidate comparisons selected for initial structural triage. “Pocket absent” means that a numerical docking run against the alternate is not scientifically meaningful; use a canonical re-docking and domain/pocket-loss figure instead.

| Priority | Gene | Canonical protein | Alternate protein / transcript | Drug | Comparison to run | Notes |
|---:|---|---|---|---|---|---|
| 1 | FYN | Full-length FYN (MANE ENST00000354650) | `transcript160449.chr6.nic` (novel): C-terminal truncation after residue 115; loses SH3, SH2, and kinase domains | Saracatinib (AZD0530) | Canonical FYN–saracatinib re-docking; alternate = pocket absent | Strongest AD-specific story. The alternate is AD-enriched in excitatory neurons (58.0% vs 2.15% control). Validate the novel transcript/protein before making a causal claim. |
| 2 | KIT | Full-length c-KIT | KIT-223: Ser715-minus splice isoform (one-residue deletion in kinase-insert region) | Masitinib | Canonical KIT vs KIT-223 comparative docking | Best conventional, same-pocket docking pair. Alternate is AD-enriched in inhibitory neurons (37.7% vs 3.5%). Expect any score difference to be modest. |
| 3 | GABRA2 | Full-length GABA-A receptor alpha-2 subunit | GABRA2-206: extensive C-terminal truncation, losing most ligand-binding and transmembrane regions | AZD7325 (BAER-101) | Canonical receptor docking; alternate = pocket/receptor absent | AD-only detection in astrocytes (9.3% vs 0%). Model the receptor as a pentamer; do not report an alternate docking score. |
| 4 | CACNA1D | Full-length CaV1.3 alpha-1D channel | CACNA1D-214: 536-aa C-terminal truncation | Isradipine | Canonical vs CACNA1D-214 structural/pocket comparison; optional docking | Strong AD switch in inhibitory neurons (28.1% vs 1.8%). The dihydropyridine pocket is likely retained, so use as a negative-control comparison. |
| 5 | BACE1 | BACE1-501 full-length active beta-secretase | BACE1-476 | Verubecestat (MK-8931) | Canonical vs alternate docking after alternate homology modeling | High-quality AD drug-failure case and experimental canonical co-crystal. Only retain if your data show the alternate is AD-enriched. |
| 6 | BACE1 | BACE1-501 full-length active beta-secretase | BACE1-457 | Verubecestat (MK-8931) | Canonical vs alternate docking after alternate homology modeling | More extensive deletion than BACE1-476, near the active site. Model quality must pass before interpreting docking. |
| 7 | CHRNA7 | Alpha-7 nicotinic acetylcholine receptor homopentamer | CHRNA7/CHRFAM7A mixed receptor complex | Encenicline | Canonical versus mixed-pentamer ligand-site comparison | Human-specific altered subunit rather than a standard splice isoform; use only if genotype/expression supports CHRFAM7A in the samples. |
| 8 | PDE9A | Full-length PDE9A catalytic isoform | The specific AD-enriched coding-altered PDE9A isoform from your data | BI 409306 (osoresnontrine) | Comparative docking only if the alternate changes/truncates the catalytic domain | Conditional candidate. Do not run an isoform that differs only in an N-terminal region outside the inhibitor pocket. |

## Recommended order for the presentation

1. FYN–saracatinib: canonical re-docking plus kinase-domain-loss figure.
2. KIT–masitinib: direct two-isoform docking comparison.
3. BACE1-501–BACE1-476/457–verubecestat: established AD reference comparison, conditional on your expression data.
4. GABRA2–AZD7325: second pocket-loss/domain-loss visual.
5. CACNA1D–isradipine and PDE9A–BI 409306: useful controls or follow-up candidates.

