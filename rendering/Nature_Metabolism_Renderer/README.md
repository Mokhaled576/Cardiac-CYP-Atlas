# Nature Metabolism figure renderer

The scripts generate eight headline-free main figures, six Extended Data figures and 75 individually exported panels. They draw from archived Fresh CYP results.

## Inputs and setup

Use the original Figure_Source_Data.zip from the submission package. Its SHA-256 is 4c85363ae711b28255a3237cc76bc9e0460290bec92425220f71fbb693e72298. The original repository and Zenodo analysis archive contain an aggregate subset. The complete renderer also needs plotting inputs from the submission archive, including donor-level expression summaries.

Use Python 3.11 or later, the pinned dependencies, Poppler and locally installed Arial fonts. Choose a new output directory for the figure files.

    python -m pip install -r requirements.txt
    python render_figures.py OUTPUT_FOLDER --source-data /path/to/Figure_Source_Data.zip
    python export_png.py OUTPUT_FOLDER

On non-Windows systems, pass --arial-dir /path/to/arial/fonts to the renderer. For a quick check, append --figures fig7 to generate Figure 7 and its eight panels. Use fig1 through fig8 or ed1 through ed6 to select other figures. SVG panels have identifiable groups. The PNG step verifies RGB color and 600-dpi metadata.

The code release includes scripts and documentation. The full plotting inputs are handled separately under the source-study access terms.
