# Human cardiac cytochrome P450 atlas

Analysis scripts and aggregate supplementary results for the eight-figure manuscript.

| Path | Contents |
| --- | --- |
| `analysis/` | Archived clustering, differential expression, pathway, donor-association and validation code. |
| `rendering/` | Code used for the main and Extended Data figure exports. |
| `data/Supplementary_Data_1.xlsx` | 1,234 eligible pseudobulk contrasts and LV DCM and RV ARVC concordance tables. |
| `data/*.tsv` | Machine-readable aggregate estimates and figure concordance tables. |

The source accessions are listed in the analysis manifest. The external validation datasets
are [GSE298023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298023)
and [GSE116250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250).
The scripts depend on the archived Fir result tree and earlier display inputs at the paths
defined in the renderers. End-to-end reproduction from raw source files has not been verified.

Release contents include analysis code and aggregate outputs. Source-study access terms apply to the original data.

## Nature Metabolism figure exports

The [revised renderer](rendering/Nature_Metabolism_Renderer/README.md) exports the eight main figures, six Extended Data figures and 75 panels. It uses archived plot inputs from the submission package. The earlier analysis archive is [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.24067872).

## Mentor revision v1.2.0

See the [revision code and reproduction notes](revision/v1.2.0/README.md) for the discovery confidence intervals, sex-adjusted and direct disease comparisons, subtype models, doublet sensitivity, and revised figures.
