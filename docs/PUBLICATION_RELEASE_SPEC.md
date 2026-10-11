# Publication Release Specification

The publication-facing repository must expose the revised 6,644-patient analysis on its default branch and in a versioned GitHub Release. It must not present the superseded 4,532-patient analysis as current.

The release version is `1.2.1`, dated 2026-10-11. The existing `v1.2.0` tag remains immutable as revision history.

The default-branch README must:

- state that TCGA-BRCA and GSE96058 had been analysed before registration;
- state that METABRIC and GSE20685 were unseen at registration;
- identify the analyses first performed after registration;
- use the revised cohort sizes, DeepSurv headline model, and primary result;
- link to the code and reproducibility release.

The repository must provide machine-checkable evidence for:

- cohort and event counts;
- all registered endpoints, including both Oncotype DX surrogate claims;
- every registered-versus-reported deviation;
- corrected stage parsing and the absence of stage in GSE96058;
- consistency among `VERSION`, `CITATION.cff`, README, and the Git tag.

The GitHub Release must expose the tagged source code and machine-readable results without attaching manuscript, Supplementary Information, or journal-correspondence files. Those documents remain in the private journal-submission package. GitHub's automatically generated source archives are the public release assets.
