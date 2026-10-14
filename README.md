# Fish-cracker correlation analysis

Reproducible correlation analysis and clustered heatmaps for ten quality attributes measured on seven formulations with three independent production batches each. Carbohydrate is excluded. This repository contains code and synthetic tests; raw research workbooks and generated results are deliberately ignored by Git.

## Setup and execution

Python 3.12 was used. Create an environment and install dependencies:

```sh
python -m venv .venv
# macOS/Linux:
source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m unittest -v test_correlation.py
```

Place these unchanged workbooks inside a local `source/` directory:

- `Fish_Crackers_Correlation_PCA_HCA_Preprocessed.xlsx`
- `Fish_Crackers_Proximate_Preprocessed.xlsx`
- `Fish_Crackers_Physical_Color_Preprocessed.xlsx`

Run:

```sh
python correlation_figures.py
python clustered_heatmap.py
python validate_results.py
```

Settings appear at the top of the scripts. Paths are relative to each script, not the working directory. Outputs go to `results/`.

## Input contract

`Batch_data` has one row per independent formulation–batch combination. Required identifiers are `Sample_ID`, `Batch_ID`, `Batch_number`, `Fish_species`, `Fish_inclusion_pct`, and `Is_control`. Batch numbers are 1–3 within A–G; Batch_ID uses `A_B1` etc.

Analyzed columns: `Protein_pct`, `Ash_pct`, `Fat_pct`, `Moisture_pct`, `Linear_expansion_pct`, `Reported_oil_absorption_pct`, `Specific_volume_cm3_g`, `Color_L`, `Color_a`, `Color_b`.

`Carbohydrate_pct` is checked for closure but not correlated. A–C are Pangas (50, 40, 30%); D–F are Tilapia (50, 40, 30%); G is fish-free. The script also verifies the supplied `PCA_scaled` and `HCA_means` sheets, matching by identifiers and named columns. It does not fit PCA or cluster formulations.

Missing/nonfinite data, duplicate keys, missing columns, mismatched cross-workbook batch keys and constant variables raise errors. No incomplete observations are silently discarded or imputed. Pairwise sample sizes are consequently 21 for A–G and 18 for A–F. Batch independence is an experimental-design statement supplied by the researcher, not something numerical validation can prove.

## Methods and outputs

- Pearson correlations and two-sided beta-null probabilities across 45 unique pairs.
- Benjamini–Hochberg adjustment over those 45 pairs; separate families for each sensitivity analysis.
- Spearman sensitivity using average ranks and approximate t probabilities (not exact permutation p-values).
- Pearson sensitivity excluding G.
- Full-precision CSV/JSON results, annotated heatmaps, and all 45 scatterplots colored/shaped by formulation.
- 350 dpi PNG and vector PDF exports.
- Variable clustering by average linkage on distance `1 - Pearson r`; identical row/column ordering. Both star-free and BH-adjusted p < 0.05 star versions are generated. Dendrograms are exploratory.

The optional `build_workbook.mjs` produces the formatted Excel workbook. It additionally requires Node.js and the Codex-provided `@oai/artifact-tool` module, which is not a Python dependency. Run `node build_workbook.mjs` after the analysis in an environment with that module. CSV/JSON and figures do not require Node. The validator checks Excel cells if that workbook exists and reports zero Excel checks otherwise.

## Interpretation limits

Pooled correlations reflect both between-formulation differences and within-formulation batch variability. Fixed formulation differences and compositional constraints limit the interpretation of nominal correlation probabilities. Removing carbohydrate does not eliminate compositional dependence. Excluding G changes the experimental range and can substantially change correlations. The oil index is calculated from frying weight change and is not verified direct oil uptake. These are exploratory associations, not causal mechanisms. With three batches per formulation, within-formulation correlations are poorly estimated.

## Review and validation

See [REVIEW.md](REVIEW.md). The pure numerical functions in `correlation_core.py` are import-safe and tested with known-answer cases. CI runs these tests without private data, including optimized Python mode. Full-data validation also independently checks correlations using `statistics.correlation`, probabilities using the F distribution, and BH-adjusted probabilities using a separate implementation. Results are recorded in `results/validation.json`.
