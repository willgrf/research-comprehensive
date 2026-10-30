# Figure 3: standard PCA of the true latent vectors and of those recovered by
# TLGC-PCA and by the latent mixture, on one data set with 40% zeros.
#   Rscript simulations/extension-example.R   -> output/figure3.{pdf,png}

suppressPackageStartupMessages(source("R/latent-mixture.R"))
source("simulations/common.R")
library(ggplot2)

n <- 600; p <- 10; K <- 3; r <- 2; sig <- 0.5; phi <- c(0.25, 0.35, 0.4)
delta <- 3; zero_prop <- 0.4

set.seed(1)
d <- simulate_mixture(n, p, K, r, sig, delta, phi, zero_prop)
fit <- fit_latent_mixture(d$X, K, r)
latent <- list("True latent Z"  = d$Z,
               "TLGC-PCA"       = latent_estimation(d$X, latent_corr(d$X))$Z,
               "Latent mixture" = fit$Z)

# Leading two PCs of the standardised columns, signed to agree with the PCs of the true Z
reference <- prcomp(d$Z, scale.=TRUE)$x[, 1:2]
pca_data <- do.call(rbind, lapply(names(latent), function(m){
  pcs <- prcomp(latent[[m]], scale.=TRUE)$x[, 1:2]
  pcs <- sweep(pcs, 2, sign(diag(cor(pcs, reference))), "*")
  data.frame(zeros=sprintf("%d%% zeros", round(100 * zero_prop)), source=m,
             PC1=pcs[, 1], PC2=pcs[, 2], group=paste("Group", d$c))
}))
pca_data$source <- factor(pca_data$source, levels=names(latent))

ink <- "#52514e"
fig <- ggplot(pca_data, aes(PC1, PC2, colour=group, shape=group)) +
  geom_point(size=0.7, alpha=0.5) +
  stat_ellipse(aes(group=group), level=0.95, linewidth=0.45, show.legend=FALSE) +
  facet_grid(zeros ~ source) +
  scale_colour_manual(values=c(`Group 1`="#2a78d6", `Group 2`="#eb6834", `Group 3`="#1baf7a")) +
  scale_shape_manual(values=c(`Group 1`=16, `Group 2`=17, `Group 3`=15)) +
  guides(colour=guide_legend(override.aes=list(size=2.2, alpha=1))) +
  labs(x="PC1", y="PC2") +
  theme_minimal(base_family="serif", base_size=9) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="#efeeea", linewidth=0.3),
        axis.text=element_text(colour=ink, size=7),
        axis.title=element_text(colour=ink, size=8.5),
        strip.text=element_text(face="bold", size=8.5),
        legend.position="top", legend.title=element_blank(),
        legend.text=element_text(size=8.5),
        panel.spacing=unit(10, "pt"))

dir.create("output", showWarnings=FALSE)
ggsave("output/figure3.pdf", fig, width=6.8, height=3.1, device=cairo_pdf)
ggsave("output/figure3.png", fig, width=6.8, height=3.1, dpi=200, bg="white")
