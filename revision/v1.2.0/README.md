# Cardiac CYP Atlas v1.2.0

Analysis scripts and aggregate outputs for the mentor revision. This archive extends the [v1.1.0 release](https://doi.org/10.5281/zenodo.24094024).

- `code/` contains R and Python scripts for approximate moderated Wald intervals, sex-adjusted and direct disease comparisons, subtype analysis, QC, doublet sensitivity and rendering.
- `aggregate/` contains summary tables and the aggregate inputs for revised Figures 3 and 4. The tables retain the original discovery testing families and the separate revision testing families.
- Figure rendering requires pandas, NumPy, ReportLab and an Arial font directory. Run `python code/render_with_intervals.py OUTPUT --source-data aggregate/figure_inputs.zip --ci aggregate/CYP_Discovery_With_CI.tsv --figures fig3 fig4`.

The analysis scripts expect existing source count and metadata archives. Obtain those from the original studies under their respective access terms. This release contains no expression matrices, donor profiles, manuscript or unpublished submission files. The original code is distributed under the repository's MIT license.
