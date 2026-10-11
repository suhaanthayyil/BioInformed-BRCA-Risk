# Publication Release Specification

The publication-facing repository must expose the revised 6,644-patient analysis on its default branch and in a versioned GitHub Release. It must not present the superseded 4,532-patient analysis as current.

The release version is `1.2.1`, dated 2026-10-11. The existing `v1.2.0` tag remains immutable as revision history.

The default-branch README must:

- state that TCGA-BRCA and GSE96058 had been analysed before registration;
- state that METABRIC and GSE20685 were unseen at registration;
- identify the analyses first performed after registration;
- use the revised cohort sizes, DeepSurv headline model, and primary result;
- link to the publication release and its manuscript and Supplementary Information assets.

The repository must provide machine-checkable evidence for:

- cohort and event counts;
- all registered endpoints, including both Oncotype DX surrogate claims;
- every registered-versus-reported deviation;
- corrected stage parsing and the absence of stage in GSE96058;
- consistency among `VERSION`, `CITATION.cff`, README, and the Git tag.

The GitHub Release must attach the revised manuscript and Supplementary Information in both DOCX and PDF form, plus SHA256 checksums and release notes. The response letter remains a journal correspondence file and is not a public release asset.
