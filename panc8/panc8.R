# Application to the CEL-seq batch of the panc8 pancreatic islet data (Section 4.4).
#   Rscript panc8/panc8.R   -> output/panc8/
#     fits.rds            the data, TLGC-PCA fit and mixture fit
#     figure5-scree.pdf   Figure 5
#     figure7-*.pdf       the four panels of Figure 7
#     table5-overlap.csv  Table 5
# The random number calls (TLGC-PCA, UMAP, mixture fit) must stay in this order
# to reproduce the report.

suppressPackageStartupMessages(library(Seurat))
library(ggplot2)
source("R/latent-mixture.R")

n_types <- 4
n_genes <- 100
out_dir <- "output/panc8"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)


# ---- Data: four most common cell types, CLR within each cell, top 100 HVGs --

counts <- Read10X("data/celseq")
annotations <- read.csv("data/celseq/annotations.csv")
cells <- CreateSeuratObject(counts)
cells$cell_type <- annotations$cell_type[match(colnames(cells), annotations$barcode)]
top_types <- names(sort(table(cells$cell_type), decreasing=TRUE))[1:n_types]
cells <- cells[, cells$cell_type %in% top_types]

# HVGs selected
cells <- NormalizeData(cells, normalization.method="CLR", margin=2)
cells <- FindVariableFeatures(cells, selection.method="dispersion", nfeatures=n_genes)
X <- t(as.matrix(GetAssayData(cells, layer="data")[VariableFeatures(cells), ]))
labels <- cells$cell_type
cat(sprintf("n = %d cells, p = %d genes, %.0f%% zeros\n", nrow(X), ncol(X), 100 * mean(X == 0)))


# ---- Fits --------------------------------------------------------------------

set.seed(1)
tlgc <- TLGC_PCA(X)
set.seed(1)
umap_tlgc <- uwot::umap(tlgc$Z)

# Starting groups
neighbours <- FindNeighbors(prcomp(tlgc$Z)$x[, 1:20], verbose=FALSE)$snn
init <- as.integer(FindClusters(neighbours, resolution=0.2, verbose=FALSE)[, 1])
if (max(init) != n_types) stop("Louvain gave ", max(init), " clusters; adjust the resolution")

mixture <- fit_latent_mixture(X, K=n_types, init=init, verbose=TRUE)
saveRDS(list(X=X, labels=labels, init=init, tlgc=tlgc, mixture=mixture), file.path(out_dir, "fits.rds"))

# Each fitted group is named by the cell type it mostly contains
groups <- table(fitted_group=mixture$labels, cell_type=labels)
majority <- colnames(groups)[apply(groups, 1, which.max)]
print(groups)
cat("ARI:", mclust::adjustedRandIndex(mixture$labels, labels), "\n")
cat("Ranks: pooled TLGC-PCA", tlgc$r, "|", paste(majority, mixture$r, collapse=", "), "\n")


# ---- Figure 5: scree plots -------------------------------------------------

n_show <- 20
panels <- c(paste0("Pooled TLGC-PCA (rank ", tlgc$r, ")"),
            paste0(majority, " group (rank ", mixture$r, ")"))
scree <- c(list(tlgc$scree), mixture$scree)
ranks <- c(tlgc$r, mixture$r)
eigenvalues <- do.call(rbind, lapply(seq_along(panels), function(i)
  data.frame(panel=panels[i], component=1:n_show, eigenvalue=scree[[i]][1:n_show],
             selected=1:n_show <= ranks[i])))
eigenvalues$panel <- factor(eigenvalues$panel, levels=panels)

ggplot(eigenvalues, aes(component, eigenvalue)) +
  geom_line(colour="grey60") +
  geom_point(aes(colour=selected), size=1.8) +
  geom_vline(data=data.frame(panel=factor(panels, levels=panels), rank=ranks),
             aes(xintercept=rank + 0.5), linetype="dashed", colour="grey40") +
  facet_wrap(~ panel, scales="free_y", nrow=1) +
  scale_colour_manual(values=c(`TRUE`="#0072B2", `FALSE`="grey60"), breaks=c(TRUE, FALSE),
                      labels=c(`TRUE`="Within selected rank", `FALSE`="Beyond selected rank")) +
  labs(x="Component", y="Eigenvalue", colour=NULL) +
  theme_bw() + theme(legend.position="top")
ggsave(file.path(out_dir, "figure5-scree.pdf"), width=10, height=2.9)


# ---- Figure 7: UMAPs ---------------------------------------------------------

set.seed(1)
umap_mix <- uwot::umap(mixture$Z)
set.seed(1)
umap_obs <- uwot::umap(X)

type_colours  <- setNames(scales::hue_pal()(ncol(groups)), colnames(groups))
group_names   <- paste0(seq_len(nrow(groups)), " (", majority, ")")
group_colours <- setNames(type_colours[majority], group_names)

umap_panel <- function(embedding, colour, colours, title, legend, file){
  fig <- ggplot(data.frame(UMAP1=embedding[, 1], UMAP2=embedding[, 2], colour=colour),
                aes(UMAP1, UMAP2, colour=colour)) +
    geom_point(size=0.8) +
    scale_colour_manual(values=colours) +
    labs(title=title, colour=legend)
  ggsave(file.path(out_dir, file), fig, width=5, height=4)
}
umap_panel(umap_obs, labels, type_colours, "Observation space: UMAP, true cell types",
           "Cell type", "figure7-observed.pdf")
umap_panel(umap_tlgc, labels, type_colours, "TLGC-PCA latent vectors: UMAP",
           "Cell type", "figure7-tlgc.pdf")
umap_panel(umap_mix, labels, type_colours, "Mixture latent vectors: UMAP, true cell types",
           "Cell type", "figure7-mixture-true.pdf")
umap_panel(umap_mix, group_names[mixture$labels], group_colours,
           "Mixture latent vectors: UMAP, fitted groups", "Fitted group", "figure7-mixture-fitted.pdf")


# ---- Table 5: subspace overlap in the mixture latent space -----------------
# Overlap ||V_a' V_b||_F^2 / min(r_a, r_b): the share of the smaller subspace in the larger

global <- eigen(cov(mixture$Z), symmetric=TRUE)
global_rank <- igraph::dim_select(global$values)
mean_differences <- sweep(mixture$mu, 2, colSums(mixture$phi * mixture$mu))
subspaces <- c(list(global$vectors[, 1:global_rank, drop=FALSE]),
               list(svd(t(mean_differences))$u[, 1:(n_types - 1)]),
               mixture$V)
names(subspaces) <- c("global", "between-group", majority)
cat(sprintf("Global rank %d; its direction(s) carry %.0f%% of the variance\n",
            global_rank, 100 * sum(global$values[1:global_rank]) / sum(global$values)))

overlap <- do.call(rbind, lapply(combn(names(subspaces), 2, simplify=FALSE), function(pair){
  A <- subspaces[[pair[1]]]; B <- subspaces[[pair[2]]]
  data.frame(a=pair[1], b=pair[2], rank_a=ncol(A), rank_b=ncol(B),
             overlap=sum(crossprod(A, B)^2) / min(ncol(A), ncol(B)))
}))
write.csv(overlap, file.path(out_dir, "table5-overlap.csv"), row.names=FALSE)
print(overlap, digits=2, row.names=FALSE)
