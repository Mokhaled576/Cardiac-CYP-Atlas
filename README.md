# Human cardiac cytochrome P450 atlas

Analysis code and aggregate results for a human cardiac single-nucleus RNA-seq study.

## Latest release · v1.2.0

- [Download the versioned archive and cite its DOI](https://doi.org/10.5281/zenodo.24268773)
- [Browse the v1.2.0 scripts and figure reproduction notes](revision/v1.2.0/README.md)

On Zenodo, start with `Cardiac_CYP_Atlas_v1.2.0_Mentor_Revision.zip`. It contains the revision scripts and aggregate tables for CYP effect-size intervals, direct disease comparisons, sex-adjusted models, lineage subclustering, quality control and doublet sensitivity, plus aggregate inputs for revised Figures 3 and 4. The same Zenodo record includes the v1.1.0 ZIP for the baseline workflow.

## Baseline workflow and archived files

| Location | Contents |
| --- | --- |
| [analysis/](analysis/) | Baseline clustering, differential expression, pathway, donor-association and external validation code. |
| [rendering/](rendering/) | Baseline main and Extended Data figure export code. |
| [data/](data/) | Baseline aggregate contrast and cross-cohort concordance tables. |

Earlier releases are available on Zenodo as [v1.1.0](https://doi.org/10.5281/zenodo.24094024) and [v1.0.0](https://doi.org/10.5281/zenodo.24067872).

## Reproduction and access

The source accessions are in the analysis manifest. External validation used [GSE298023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298023) and [GSE116250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250). Some scripts require the archived Fir result tree or plot inputs. End-to-end reproduction from raw source files has not been verified. The Zenodo MIT license applies to the original code. Source datasets retain their original access terms.
