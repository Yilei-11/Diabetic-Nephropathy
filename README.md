# dn-gene-panel-moea

Multi-objective evolutionary search for sparse, biologically coherent gene panels for diabetic nephropathy (DN) classification, evaluated with leakage-controlled nested cross-validation and external cohort transfer.

Six MOEAs, spanning three mechanism families, are implemented as `pymoo` extensions and benchmarked against each other, against random/top-DEG panel controls, and against standard feature-selection baselines:

| Family               | Classic baseline | Variant                          |
|----------------------|------------------|----------------------------------|
| Dominance + crowding | NSGA-II          | NSGA-II with balanced tie-break  |
| Strength + archive   | SPEA2            | SPEA2+                           |
| Decomposition        | MOEA/D           | MOEA/D-NRD                       |

Objectives: panel leanness `1/(n+1)`, inner-CV balanced accuracy of an L2 logistic regression, and size-normalized Gene Ontology (Biological Process) coherence (mean pairwise Jaccard similarity of GO memberships).

## Repository layout

```
src/dn_moea/
  geo_data.py           GEO download, probe-to-gene mapping, labels
  preprocessing.py      Welch t-test + BH-FDR candidate-pool screening
  objectives.py         leanness, inner-CV balanced accuracy, GO coherence
  sampling.py           harmonic sparse initializer and ablation controls
  algorithms/           NSGA-II tie-break, SPEA2+, MOEA/D-NRD
scripts/
  download_genesets.py  fetch Enrichr GMT libraries
  experiments/          main nested experiment and initialization ablation
  analysis/             robustness experiments, baselines, uncertainty,
                        annotation analyses, aggregation and tests
  figures/              study flowchart, gene recurrence, robustness figures
tests/                  offline unit and integration tests
```

## Install

```bash
pip install -r requirements.txt
```

## Data

All data are public and are not stored in this repository.

- **GSE96804**: discovery cohort (41 DN vs. 20 control glomeruli). Downloaded automatically with [`GEOparse`](https://github.com/guma44/GEOparse) into `data/external/GEO/`.
- **GSE30528**: external cohort (9 DN vs. 13 control glomeruli). Downloaded automatically.
- **GSE142025**: RNA-seq external cohort. GEO does not embed per-sample RNA-seq values in the SOFT file. Extract the gene-symbol matrix from `GSE142025_RAW.tar` and save it to `data/external/GEO/GSE142025_expr.pkl` (see `dn_moea.geo_data.load_gse142025`). Only `scripts/analysis/uncertainty.py` needs it.
- **Gene-set libraries** (Enrichr GO_Biological_Process_2023, Reactome_2022, KEGG_2021_Human, MSigDB_Hallmark_2020):

```bash
python scripts/download_genesets.py   # -> data/external/genesets/*.gmt
```

## Main experiment

Design: four repeats of stratified five-fold outer CV on GSE96804 (20 outer splits). Differential-expression screening is repeated inside each outer training split. The search uses harmonic sparse initialization with no panel-size cap. It is run with real GO, without GO, or with randomized GO (annotation profiles permuted among genes). A representative panel is chosen accuracy-first and GO-blind, and same-size random and top-DEG panels serve as controls.

The 160-task matrix covers six algorithms × 20 splits with real GO, plus NSGA-II without GO and NSGA-II with randomized GO. Run one task per index and merge the outputs:

```bash
for k in $(seq 0 159); do
  python scripts/experiments/run_sparse_go_hpc_task.py --task-index $k \
    --output-dir results/sparse_go_pairwise_final_20260920/tasks
done
python scripts/experiments/collect_tasks.py \
  --task-dir results/sparse_go_pairwise_final_20260920/tasks \
  --output-dir results/sparse_go_pairwise_final_20260920
```

This writes `per_split.csv`, `pareto_archives.csv`, `summary.csv`, `go_ablation.csv`, `feature_stability.csv` and `archive_summary.csv`.

Initialization ablation for NSGA-II has 80 tasks: harmonic, mean-matched Bernoulli, i.i.d. harmonic and dense Bernoulli initialization, each on 20 splits.

```bash
python scripts/experiments/run_initialization_ablation_task.py --task-index K \
  --output-dir results/initialization_ablation_20260920
```

`scripts/experiments/run_nested_ablation_study.py` runs the same protocol in one process with configurable settings (`--help`).

## Robustness and comparison analyses

| Script | Analysis |
|---|---|
| `scripts/analysis/nested_task.py` | One nested run per manifest row: additional search seeds, outer repeats 5–10, permuted-label null, inner-fold DEG screening, unrestricted gene universe, bounded uniform and SparseEA initializers, and the complete SparseEA algorithm |
| `scripts/analysis/make_manifest.py` | Writes the 960-row task manifest for `nested_task.py` |
| `scripts/analysis/baselines.py` | LASSO, elastic net, RFE, RFECV and stability selection on the same outer splits |
| `scripts/analysis/build_panel_table.py`, `scripts/analysis/uncertainty.py` | Sensitivity, specificity, confidence intervals, refitted and locked external transfer |
| `scripts/analysis/go_rule_annotation.py`, `scripts/analysis/enrichment_sensitivity.py` | GO-aware tolerance rule; Reactome, KEGG and Hallmark coherence and over-representation |
| `scripts/analysis/aggregate.py`, `scripts/analysis/external_tests.py` | Aggregation and corrected paired tests |
| `scripts/figures/` | Study flowchart, gene recurrence heatmap and robustness figure (set `TNR_FONT_DIR` to a folder with Times New Roman TTF files) |

Typical order, after the main experiment and the initialization ablation:

```bash
python scripts/analysis/make_manifest.py --output results/manifest.csv
bash scripts/analysis/run_manifest.sh results/manifest.csv results/nested_tasks 16
python scripts/analysis/baselines.py \
  --per-split results/sparse_go_pairwise_final_20260920/per_split.csv \
  --output results/baselines/baselines_per_split.csv
python scripts/analysis/build_panel_table.py results/uncertainty/panels_main.csv
python scripts/analysis/uncertainty.py --panels results/uncertainty/panels_main.csv \
  --output-dir results/uncertainty/all
python scripts/analysis/go_rule_annotation.py
python scripts/analysis/aggregate.py
python scripts/analysis/enrichment_sensitivity.py
python scripts/analysis/external_tests.py
python scripts/figures/study_flowchart.py study_flowchart.pdf
python scripts/figures/gene_recurrence.py gene_recurrence.pdf
python scripts/figures/robustness.py robustness.pdf
```

Parsed GEO matrices are cached in `data/cache/`. Results are written below `results/`.

## Seeds

| Purpose | Seed |
|---|---|
| Outer resampling | 20260920 |
| Search | 20261000 + split (+1000 × run for additional runs) |
| Random panels | 20262000 |
| Randomized GO | 20263000 + split |
| Label permutations | 20264000 + permutation |
| Random coherence panels | 20265000 |

## Tests

```bash
pytest
```

The test suite runs offline on synthetic data and does not download anything from GEO.

## Software

Python 3.14.6, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6, scikit-learn 1.9.1, pymoo 0.6.2, statsmodels 0.15.0, GEOparse 2.0.4, Matplotlib 3.11.2.

## License

MIT — see `LICENSE`.
