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

The repository is intended for private manuscript preparation. Raw expression matrices,
donor-level profiles, and the unpublished manuscript are outside this deposit.
