# Code review — 2026-10-14

Scope: correlation computation, source reconciliation, figures, variable clustering, Excel export and numerical validator. Reviewed locally; the connected Code Review tool exposes pull-request CI diagnostics only and was not represented as a separate local reviewer.

## Findings corrected

1. **P1 — Production validation disabled by `python -O`.** Source and output checks used assertions. Replaced with explicit `require` checks that raise ValueError in optimized execution too.
2. **P2 — Duplicate source batch identifiers could be overwritten.** Dictionary construction silently retained the last matching row. Added explicit duplicate detection and equality checks for complete source key sets.
3. **P2 — Preprocessing reconciliation depended on row/column order.** Supplied scaled data and means were compared positionally. These comparisons now align named columns and formulation/batch keys; partial records and duplicate/missing headers fail explicitly.
4. **P2 — Cluster annotations hard-coded the observed significant-pair count.** The script expected exactly 11 significant pairs. It now calculates the count from current adjusted probabilities while checking the expected comparison-family size.
5. **P2 — Core statistics could not be tested without running the full plotting script.** Extracted import-safe numerical functions and added known-answer tests, including tied ranks, BH ties, zero/perfect correlation, a known p-value, invalid data and a linkage example.
6. **P2 — Validator required a proprietary-tool-produced Excel snapshot.** Numerical and figure validation now works without the optional workbook; workbook checks are explicitly counted only when it exists. Perfect correlations are handled without division by zero.

## Validation performed

- Eight regression tests passed in normal and optimized Python modes.
- Full-data workflow ran under optimized Python.
- All 135 correlation-result records exactly matched the prior validated results.
- Independent audit checked 135 coefficient/probability/BH results and 540 Excel numeric cells.
- Eight original PNG files, three heatmap PDF pages and five scatterplot PDF pages passed integrity checks.
- Both clustered heatmap variants regenerated; the known-answer linkage test passed.

The changes fix reuse and validation problems; they did not change current correlation estimates or conclusions. The implementation remains tailored to this study's 21 complete batches and ten variables. Spearman probabilities remain explicitly approximate. Pooled formulation effects, compositional dependence and the weight-change oil index remain scientific interpretation limits, not bugs that code can remove.

Raw measurements and generated results are kept local, excluded from the public repository.
