# Analysis of Truncated Gaussian Copula Methods for Dimension Reduction

Code for the simulations and data analysis in the research comprehensive


## Layout

```
R/
  tlgc-pca.R              TLGC-PCA: latent correlation (bridge functions), latent vectors, TLGC_PCA()
  latent-mixture.R        Latent Gaussian mixture extension: fit_latent_mixture()
simulations/
  common.R                Shared helpers (chordal distance, transformations, mixture data, ...)
  replication.R           Section 3: TLGC model (3.2) and Poisson model (3.3) studies
  replication-results.R   Tables 1-3 and 6-8
  replication-figure.py   Figure 2
  extension.R             Section 4.3 mixture study and Appendix D.2 transformation study
  extension-results.R     Tables 4 and 9, Figures 4 and 6
  extension-example.R     Figure 3
panc8/
  panc8.R                 Section 4.4: Figure 5, Figure 7 panels, Table 5
data/celseq/              CEL-seq batch of the panc8 data (Grün et al., 2016), 10x format with cell-type annotations
```

Results are written to `output/`.

## Reproducing the report

| Report item | Command | Output |
|---|---|---|
| Section 3.2 simulation | `Rscript simulations/replication.R --study=tlgc --cores=8` | `output/replication.csv` |
| Section 3.3 simulation | `Rscript simulations/replication.R --study=poisson --cores=8` | `output/poisson.csv` |
| Tables 1, 2, 3, 6-8 | `Rscript simulations/replication-results.R` | `output/replication-tables.tex` |
| Figure 2 | `python3 simulations/replication-figure.py` | `output/figure2.pdf` |
| Section 4.3 simulation (0-40% zeros) | `Rscript simulations/extension.R --study=mixture --cores=7` | `output/mixture/` |
| Section 4.3 simulation (50-60% zeros) | `Rscript simulations/extension.R --study=mixture --cores=7 --zero_prop=0.5,0.6 --seed=2 --out_dir=output/mixture-high-zeros` | `output/mixture-high-zeros/` |
| Appendix D.2 simulation | `Rscript simulations/extension.R --study=transformations --cores=7` | `output/transformations.csv` |
| Table 4, Figure 4, Table 9, Figure 6 | `Rscript simulations/extension-results.R` | `output/table4.tex`, `figure4.pdf`, `table9.tex`, `figure6.pdf` |
| Figure 3 | `Rscript simulations/extension-example.R` | `output/figure3.pdf` |
| Section 4.4, Figure 5, Table 5, Figure 7 | `Rscript panc8/panc8.R` | `output/panc8/` |

Figure 7 in the report is the four `figure7-*.pdf` panels arranged in a 2 x 2 grid.

Every simulation script runs with fixed seeds, and the results do not depend on
the number of cores (`--cores` uses forked processes, so on Windows use `--cores=1`).
Any argument of the `run_*` functions can be passed on the command line as
`--name=value`, with comma-separated values for vectors. For example, to run a
quick test of the Section 3.2 study:

```
Rscript simulations/replication.R --study=tlgc --num_reps=2 --transf=mixed --trunc_prop=0,0.5 --file=output/test.csv
```

The full simulation studies are slow. Most of the time goes to the Gibbs sampling of the latent zeros and to
ZIPFA at high zero proportions.

## Requirements

R (>= 4.3) with `mnormt`, `igraph`, `MASS`, `PoissonPCA`, `ZIPFA`, `ggplot2`, `scales`,
`Seurat` (>= 5), `uwot` and `mclust`; Python 3 with `pandas` and `matplotlib` for Figure 2.
Tested with R 4.6.1, mnormt 2.1.2, igraph 2.3.3,
PoissonPCA 1.0.3, ZIPFA 0.8.1, Seurat 5.5.1 and uwot 0.2.5.
