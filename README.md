# dn-gene-panel-moea

Multi-objective evolutionary search for sparse, biologically coherent gene panels for diabetic nephropathy (DN) classification, evaluated with leakage-controlled nested cross-validation and external cohort transfer.

Six MOEAs, spanning three mechanism families, are implemented as `pymoo` extensions and benchmarked against each other and against random/top-DEG panel baselines:

| Family              | Classic baseline | Variant                       |
|---------------------|-------------------|-------------------------------|
| Dominance + crowding| NSGA-II           | NSGA-II with balanced tie-break |
| Strength + archive  | SPEA2             | SPEA2+                        |
| Decomposition       | MOEA/D            | MOEA/D-NRD                    |

Objectives: panel leanness, cross-validated balanced classification accuracy, and Gene Ontology (Biological Process) coherence.

## Install

```bash
pip install -r requirements.txt
```

## Data

All gene-expression data is public and is downloaded automatically on first use via [`GEOparse`](https://github.com/guma44/GEOparse) into `data/external/GEO/` (not checked into this repo):

- **GSE96804** — primary training/search cohort (41 DN vs. 20 control glomeruli)
- **GSE30528** — external validation cohort (9 DN vs. 13 control glomeruli)

A Gene Ontology Biological Process 2023 gene-set library (`data/external/genesets/GO_Biological_Process_2023.gmt`, from [Enrichr](https://maayanlab.cloud/Enrichr/#libraries)) is included for the pathway-coherence objective.

A third cohort, **GSE142025** (RNA-seq), is supported by `dn_moea.geo_data.load_gse142025` but requires its expression matrix to be extracted separately from that series' supplementary RNA-seq archive (GEO does not embed per-sample RNA-seq values in the SOFT file); this step is not included here.

## Usage

Run a single algorithm/seed panel search (writes the Pareto set and external-validation metrics to `--output-dir`):

```bash
python run_dk_panel_task.py --algorithm nsga2_tiebreak --seed 0 \
  --pop-size 100 --max-evals 10000 --max-pool-size 2000 --max-panel-size 30 \
  --output-dir results/panel_search
```

`--algorithm` is one of `nsga2`, `nsga2_tiebreak`, `spea2`, `spea2_plus`, `moead`, `moead_nrd`.

Run leakage-controlled nested cross-validation for one algorithm (DEG screening and panel search happen strictly inside each outer training fold):

```bash
python run_nested_cv.py --algorithm nsga2_tiebreak \
  --n-splits 5 --pop-size 100 --n-gen 100 --output-dir results/nested_cv
```

## Tests

```bash
pytest
```

The test suite is fully offline (no GEO download required).

## License

MIT — see `LICENSE`.

## Final nested experiment

`scripts/experiments/run_nested_ablation_study.py` implements the final leakage-controlled design reported in the article: four repeats of stratified five-fold outer cross-validation on GSE96804, split-specific differential-expression screening, harmonic sparse initialization without a panel-size cap, the three-objective search with real, removed, or randomized GO, the accuracy-first representative rule, and same-size random and top-DEG controls. `scripts/experiments/run_sparse_go_hpc_task.py --task-index K` runs one of the 160 tasks of the main matrix, and `scripts/experiments/run_initialization_ablation_task.py` runs the initialization ablation.

## Revision analyses

The folder `revision1/` contains the analyses added during peer review.

| Script | Analysis |
|---|---|
| `revision1/scripts/rev_task.py` | One nested run per manifest row: additional search seeds, outer repeats 5-10, permuted-label null, inner-fold DEG screening, unrestricted gene universe, bounded uniform and SparseEA initializers, and the complete SparseEA algorithm |
| `revision1/scripts/rev_baselines.py` | LASSO, elastic net, RFE, RFECV, and stability selection on the same outer splits |
| `revision1/scripts/build_panel_table.py`, `revision1/scripts/rev_uncertainty.py` | Sensitivity, specificity, confidence intervals, refitted and locked external transfer |
| `revision1/scripts/rev_go_rule_annotation.py`, `revision1/scripts/rev_enrichment_sensitivity.py` | GO-aware tolerance rule; Reactome, KEGG and Hallmark coherence and over-representation |
| `revision1/scripts/rev_aggregate.py`, `revision1/scripts/rev_external_tests.py` | Aggregation and corrected paired tests |
| `revision1/scripts/figures/` | Figures 1 and 7 and Supplementary Figure S1 (set `TNR_FONT_DIR` to a folder with Times New Roman TTF files) |

Typical order:

```bash
bash revision1/run_manifest.sh revision1/revision_main.csv revision1/results/main 16
python revision1/scripts/rev_baselines.py --per-split results/sparse_go_pairwise_final_20260920/per_split.csv \
  --output revision1/results/baselines/baselines_per_split.csv
python revision1/scripts/build_panel_table.py revision1/results/uncertainty/panels_main.csv
python revision1/scripts/rev_go_rule_annotation.py
python revision1/scripts/rev_aggregate.py
```

The main-experiment outputs (`per_split.csv`, `pareto_archives.csv`) are expected in `results/sparse_go_pairwise_final_20260920/`; they are distributed with the article as Supplementary Data. Seeds: outer resampling 20260920; search 20261000 + split (+1000 x run for additional runs); random panels 20262000; randomized GO 20263000 + split; label permutations 20264000 + permutation; random coherence panels 20265000.

The gene-set libraries in `revision1/genesets/` are the Enrichr files GO_Biological_Process_2023, Reactome_2022, KEGG_2021_Human and MSigDB_Hallmark_2020.

Software used for the article: Python 3.14.6, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6, scikit-learn 1.9.1, pymoo 0.6.2, statsmodels 0.15.0, GEOparse 2.0.4, Matplotlib 3.11.2.
