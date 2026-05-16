# Phase 4 Pathway Scoring QC

Generated: 2026-05-16T12:31:52-04:00
MSigDB version: v2024.1.Hs

## Gene Set Inventory

|   hallmark |   reactome_curated |   kegg_medicus_curated |   mandatory_aggregation_sets |
|-----------:|-------------------:|-----------------------:|-----------------------------:|
|         50 |                244 |                     50 |                           29 |

## Cohort Output Shapes

| cohort    |   n_samples_pathway_scores |   n_pathway_sets |   n_samples_features |   n_features |
|:----------|---------------------------:|-----------------:|---------------------:|-------------:|
| TCGA-BRCA |                        213 |              345 |                  213 |            7 |
| GSE96058  |                       3409 |              347 |                 3409 |            7 |
| METABRIC  |                       1980 |              347 |                 1980 |            7 |
| GSE20685  |                        327 |              344 |                  327 |            7 |

## Seven-Pathway Aggregation Coverage

| cohort    | feature                  |   n_components_present | components_present                                                                                                                                                                                                        |
|:----------|:-------------------------|-----------------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| TCGA-BRCA | Pathway_Immune           |                      5 | HALLMARK_INTERFERON_GAMMA_RESPONSE;HALLMARK_INTERFERON_ALPHA_RESPONSE;HALLMARK_INFLAMMATORY_RESPONSE;HALLMARK_ALLOGRAFT_REJECTION;HALLMARK_COMPLEMENT                                                                     |
| TCGA-BRCA | Pathway_Proliferation    |                      5 | HALLMARK_E2F_TARGETS;HALLMARK_G2M_CHECKPOINT;HALLMARK_MYC_TARGETS_V1;HALLMARK_MYC_TARGETS_V2;HALLMARK_MITOTIC_SPINDLE                                                                                                     |
| TCGA-BRCA | Pathway_DNA_Repair       |                      6 | HALLMARK_DNA_REPAIR;REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR;REACTOME_NUCLEOTIDE_EXCISION_REPAIR;REACTOME_MISMATCH_REPAIR;KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION;KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR |
| TCGA-BRCA | Pathway_Metabolism       |                      4 | HALLMARK_OXIDATIVE_PHOSPHORYLATION;HALLMARK_GLYCOLYSIS;HALLMARK_FATTY_ACID_METABOLISM;KEGG_MEDICUS_REFERENCE_GLYCOLYSIS                                                                                                   |
| TCGA-BRCA | Pathway_Stromal_EMT      |                      3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;HALLMARK_ANGIOGENESIS;HALLMARK_HEDGEHOG_SIGNALING                                                                                                                              |
| TCGA-BRCA | Pathway_Apoptosis_Stress |                      3 | HALLMARK_APOPTOSIS;HALLMARK_P53_PATHWAY;HALLMARK_HYPOXIA                                                                                                                                                                  |
| TCGA-BRCA | Pathway_Hormone          |                      3 | HALLMARK_ESTROGEN_RESPONSE_EARLY;HALLMARK_ESTROGEN_RESPONSE_LATE;HALLMARK_ANDROGEN_RESPONSE                                                                                                                               |
| GSE96058  | Pathway_Immune           |                      5 | HALLMARK_INTERFERON_GAMMA_RESPONSE;HALLMARK_INTERFERON_ALPHA_RESPONSE;HALLMARK_INFLAMMATORY_RESPONSE;HALLMARK_ALLOGRAFT_REJECTION;HALLMARK_COMPLEMENT                                                                     |
| GSE96058  | Pathway_Proliferation    |                      5 | HALLMARK_E2F_TARGETS;HALLMARK_G2M_CHECKPOINT;HALLMARK_MYC_TARGETS_V1;HALLMARK_MYC_TARGETS_V2;HALLMARK_MITOTIC_SPINDLE                                                                                                     |
| GSE96058  | Pathway_DNA_Repair       |                      6 | HALLMARK_DNA_REPAIR;REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR;REACTOME_NUCLEOTIDE_EXCISION_REPAIR;REACTOME_MISMATCH_REPAIR;KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION;KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR |
| GSE96058  | Pathway_Metabolism       |                      4 | HALLMARK_OXIDATIVE_PHOSPHORYLATION;HALLMARK_GLYCOLYSIS;HALLMARK_FATTY_ACID_METABOLISM;KEGG_MEDICUS_REFERENCE_GLYCOLYSIS                                                                                                   |
| GSE96058  | Pathway_Stromal_EMT      |                      3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;HALLMARK_ANGIOGENESIS;HALLMARK_HEDGEHOG_SIGNALING                                                                                                                              |
| GSE96058  | Pathway_Apoptosis_Stress |                      3 | HALLMARK_APOPTOSIS;HALLMARK_P53_PATHWAY;HALLMARK_HYPOXIA                                                                                                                                                                  |
| GSE96058  | Pathway_Hormone          |                      3 | HALLMARK_ESTROGEN_RESPONSE_EARLY;HALLMARK_ESTROGEN_RESPONSE_LATE;HALLMARK_ANDROGEN_RESPONSE                                                                                                                               |
| METABRIC  | Pathway_Immune           |                      5 | HALLMARK_INTERFERON_GAMMA_RESPONSE;HALLMARK_INTERFERON_ALPHA_RESPONSE;HALLMARK_INFLAMMATORY_RESPONSE;HALLMARK_ALLOGRAFT_REJECTION;HALLMARK_COMPLEMENT                                                                     |
| METABRIC  | Pathway_Proliferation    |                      5 | HALLMARK_E2F_TARGETS;HALLMARK_G2M_CHECKPOINT;HALLMARK_MYC_TARGETS_V1;HALLMARK_MYC_TARGETS_V2;HALLMARK_MITOTIC_SPINDLE                                                                                                     |
| METABRIC  | Pathway_DNA_Repair       |                      6 | HALLMARK_DNA_REPAIR;REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR;REACTOME_NUCLEOTIDE_EXCISION_REPAIR;REACTOME_MISMATCH_REPAIR;KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION;KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR |
| METABRIC  | Pathway_Metabolism       |                      4 | HALLMARK_OXIDATIVE_PHOSPHORYLATION;HALLMARK_GLYCOLYSIS;HALLMARK_FATTY_ACID_METABOLISM;KEGG_MEDICUS_REFERENCE_GLYCOLYSIS                                                                                                   |
| METABRIC  | Pathway_Stromal_EMT      |                      3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;HALLMARK_ANGIOGENESIS;HALLMARK_HEDGEHOG_SIGNALING                                                                                                                              |
| METABRIC  | Pathway_Apoptosis_Stress |                      3 | HALLMARK_APOPTOSIS;HALLMARK_P53_PATHWAY;HALLMARK_HYPOXIA                                                                                                                                                                  |
| METABRIC  | Pathway_Hormone          |                      3 | HALLMARK_ESTROGEN_RESPONSE_EARLY;HALLMARK_ESTROGEN_RESPONSE_LATE;HALLMARK_ANDROGEN_RESPONSE                                                                                                                               |
| GSE20685  | Pathway_Immune           |                      5 | HALLMARK_INTERFERON_GAMMA_RESPONSE;HALLMARK_INTERFERON_ALPHA_RESPONSE;HALLMARK_INFLAMMATORY_RESPONSE;HALLMARK_ALLOGRAFT_REJECTION;HALLMARK_COMPLEMENT                                                                     |
| GSE20685  | Pathway_Proliferation    |                      5 | HALLMARK_E2F_TARGETS;HALLMARK_G2M_CHECKPOINT;HALLMARK_MYC_TARGETS_V1;HALLMARK_MYC_TARGETS_V2;HALLMARK_MITOTIC_SPINDLE                                                                                                     |
| GSE20685  | Pathway_DNA_Repair       |                      6 | HALLMARK_DNA_REPAIR;REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR;REACTOME_NUCLEOTIDE_EXCISION_REPAIR;REACTOME_MISMATCH_REPAIR;KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION;KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR |
| GSE20685  | Pathway_Metabolism       |                      4 | HALLMARK_OXIDATIVE_PHOSPHORYLATION;HALLMARK_GLYCOLYSIS;HALLMARK_FATTY_ACID_METABOLISM;KEGG_MEDICUS_REFERENCE_GLYCOLYSIS                                                                                                   |
| GSE20685  | Pathway_Stromal_EMT      |                      3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;HALLMARK_ANGIOGENESIS;HALLMARK_HEDGEHOG_SIGNALING                                                                                                                              |
| GSE20685  | Pathway_Apoptosis_Stress |                      3 | HALLMARK_APOPTOSIS;HALLMARK_P53_PATHWAY;HALLMARK_HYPOXIA                                                                                                                                                                  |
| GSE20685  | Pathway_Hormone          |                      3 | HALLMARK_ESTROGEN_RESPONSE_EARLY;HALLMARK_ESTROGEN_RESPONSE_LATE;HALLMARK_ANDROGEN_RESPONSE                                                                                                                               |

## Gene Coverage Summary

| cohort    |   pathways_scored |   median_coverage |   min_coverage |
|:----------|------------------:|------------------:|---------------:|
| GSE20685  |               344 |             0.973 |          0.106 |
| GSE96058  |               347 |             1.000 |          0.106 |
| METABRIC  |               347 |             1.000 |          0.349 |
| TCGA-BRCA |               345 |             0.980 |          0.106 |

Scores use a deterministic rank-percentile ssGSEA-style implementation. The R/GSVA entrypoint remains available for package-based cross-checking.
