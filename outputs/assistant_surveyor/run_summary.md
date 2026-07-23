# ASSISTANT_SURVEYOR — Run Summary
Generated: 2026-07-10 15:44:05
Runtime: 36.2 s

## Input
- Hits: 1,201 rows (845 genes, 1,140 transcripts)

## Junior-Layer Gate (by biotype_class)
| Outcome | Hits | Genes |
|---------|------|-------|
| Pass | 926 | 690 |
| Drop | 275 | 243 |

## biotype_class Distribution
| biotype_class | Count | Junior Pass? |
|---------------|-------|--------------|
| PC_UTR | 481 | yes |
| PC_CDS | 297 | yes |
| RI | 192 | no |
| novel | 133 | yes |
| NMD | 82 | no |
| PC_CDS_ND | 15 | yes |
| other | 1 | no |

## proxy_type Distribution
| proxy_type | Count |
|------------|-------|
| N | 674 |
| D | 431 |
| NMD | 82 |
| C | 14 |

## OT Label Distribution
| OT Label | Count |
|----------|-------|
| novel | 989 |
| emerging | 113 |
| supported | 99 |

## AD Prior Hits
| Category | Count |
|----------|-------|
| none | 1197 |
| gwas | 2 |
| pathway | 2 |

## Top 10 Junior-Pass Candidates (by |Δ usage|)
| gene | transcript | cell_type | dPSI | biotype | proxy | OT score | OT label | AD prior |
|------|-----------|-----------|------|---------|-------|----------|----------|----------|
| ZNF558 | ZNF558-201 | OPC | +0.921 | PC_UTR | N | 0.000 | novel | none |
| ICA1L | ICA1L-201 | Inhibitory_neuron | -0.906 | PC_UTR | N | 0.179 | emerging | none |
| ZFAND3 | transcript52114.chr6.nic | OPC | +0.890 | novel | D | 0.000 | novel | none |
| CLUH | CLUH-215 | Excitatory_neuron | -0.882 | PC_CDS | D | 0.000 | novel | none |
| DKC1 | DKC1-201 | Astrocyte | -0.864 | PC_UTR | N | 0.010 | novel | none |
| EDEM3 | EDEM3-201 | OPC | -0.858 | PC_UTR | N | 0.000 | novel | none |
| SIRT7 | SIRT7-201 | Inhibitory_neuron | -0.853 | PC_UTR | N | 0.000 | novel | none |
| ADAMTS10 | transcript17744.chr19.nnic | Excitatory_neuron | +0.848 | novel | C | 0.005 | novel | none |
| NCOA6 | NCOA6-201 | Excitatory_neuron | -0.844 | PC_UTR | N | 0.000 | novel | none |
| CRAT | CRAT-201 | Astrocyte | -0.844 | PC_UTR | N | 0.254 | supported | none |

## Cell-Type Breakdown (Junior Pass)
| Cell type | Pass hits |
|-----------|-----------|
| Inhibitory_neuron | 237 |
| Excitatory_neuron | 210 |
| Oligodendrocyte | 184 |
| OPC | 160 |
| Astrocyte | 100 |
| Microglia | 18 |
| Vascular_cell | 17 |
