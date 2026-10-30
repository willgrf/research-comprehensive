# Tables and figures for the extension simulations:
#   Table 4, Figure 4 from output/mixture and output/mixture-high-zeros
#   Table 9, Figure 6 from output/transformations.csv
#   Rscript simulations/extension-results.R
#   -> output/{table4.tex, figure4.pdf, figure4.png, table9.tex, figure6.pdf}

library(ggplot2)

studies <- c("output/mixture", "output/mixture-high-zeros")
studies <- studies[file.exists(file.path(studies, "replicates.csv"))]
read_study <- function(file) do.call(rbind, lapply(studies, function(s)
  read.csv(file.path(s, file), stringsAsFactors=FALSE)))
reps <- read_study("replicates.csv")
comps <- read_study("components.csv")

reps$n <- factor(reps$n, levels=c(300, 600, 1200))
reps$delta_label <- factor(sprintf("κ = %g", reps$delta), levels=sprintf("κ = %g", sort(unique(reps$delta))))
reps$zeros <- 100 * reps$zero_prop


# ---- Table 4 ---------------------------------------------------------------

deltas <- sort(unique(reps$delta))
rows <- character(0)
for (z in sort(unique(reps$zero_prop))) for (nn in c(300, 600, 1200)){
  cells <- unlist(lapply(deltas, function(dl){
    r <- reps[reps$delta == dl & reps$zero_prop == z & reps$n == nn, ]
    cm <- comps[comps$delta == dl & comps$zero_prop == z & comps$n == nn, ]
    c(sprintf("%.3f", mean(r$accuracy, na.rm=TRUE)), sprintf("%.3f", median(cm$chordal, na.rm=TRUE)),
      sprintf("%.0f\\%%", 100 * mean(cm$rank_hat == 2, na.rm=TRUE)))
  }))
  if (nn == 300 && z > 0) rows <- c(rows, "\\addlinespace")
  rows <- c(rows, sprintf("%s & %d & %s \\\\", if (nn == 300) paste0(100 * z, "\\%") else "",
                          nn, paste(cells, collapse=" & ")))
}
writeLines(c("\\begin{table}[ht]", "\\centering", "\\small",
             sprintf("\\begin{tabular}{@{}rr%s@{}}", strrep("ccc", length(deltas))), "\\toprule",
             paste0(" & & ", paste(sprintf("\\multicolumn{3}{c}{$\\kappa = %g$}", deltas), collapse=" & "), " \\\\"),
             paste(sprintf("\\cmidrule(lr){%d-%d}", 3 * seq_along(deltas), 3 * seq_along(deltas) + 2), collapse=""),
             paste0("Zeros & $n$ & ", paste(rep("Accuracy & Distance & Rank 2", length(deltas)), collapse=" & "), " \\\\"),
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
             paste("\\caption{Recovery of the latent mixture ($K = 3$, $p = 10$, $r_k = 2$, 50 replications per setting).",
                   "Accuracy is the mean over replications of the proportion of samples assigned to the correct group.",
                   "Distance is the median chordal distance between each group's fitted and true $\\mathbf{V}_k$,",
                   "over all groups and replications. Rank 2 is the percentage of groups, across all replications,",
                   "for which the true rank was selected.}"),
             "\\label{tab:mixture-recovery}", "\\end{table}"),
           "output/table4.tex")


# ---- Figure 4: separation retained at n = 1200 ------------------------------

ink <- "#52514e"; muted <- "#8a8984"
method_colours <- c(`True latent Z`=ink, `TLGC-PCA`="#2a78d6", `Latent mixture`="#eb6834")

sep <- rbind(data.frame(reps[c("delta_label", "zeros", "n")], method="TLGC-PCA", retained=reps$bh_tlgc / reps$bh_true),
             data.frame(reps[c("delta_label", "zeros", "n")], method="Latent mixture", retained=reps$bh_mix / reps$bh_true))
sep <- sep[sep$n == 1200, ]
s <- aggregate(retained ~ delta_label + zeros + n + method, sep, mean)
s$method <- factor(s$method, levels=names(method_colours))

fig <- ggplot(s, aes(zeros, retained, colour=method, shape=method)) +
  geom_hline(yintercept=1, colour=muted, linewidth=0.4, linetype="22") +
  geom_line(linewidth=0.6) + geom_point(size=1.8) +
  facet_grid(n ~ delta_label, labeller=labeller(n=function(x) paste("n =", x))) +
  scale_x_continuous(breaks=sort(unique(reps$zeros)), labels=function(x) paste0(x, "%")) +
  scale_colour_manual(values=method_colours, name=NULL) +
  scale_shape_manual(values=c(`TLGC-PCA`=16, `Latent mixture`=17), name=NULL) +
  scale_y_continuous(labels=function(x) paste0(round(100 * x), "%")) +
  labs(x="Proportion of zeros", y="Separation retained") +
  theme_minimal(base_family="serif", base_size=9) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="#efeeea", linewidth=0.3),
        axis.text=element_text(colour=ink, size=7.5),
        axis.title=element_text(colour=ink, size=8.5),
        strip.text=element_text(face="bold", size=8.5),
        legend.position="top", legend.text=element_text(size=8),
        legend.title=element_text(size=8, colour=ink),
        panel.spacing=unit(12, "pt"))
ggsave("output/figure4.pdf", fig, width=6.5, height=2.6, device=cairo_pdf)
ggsave("output/figure4.png", fig, width=6.5, height=2.6, dpi=200, bg="white")


# ---- Table 9 and Figure 6: marginal transformations ------------------------

d <- read.csv("output/transformations.csv")
d$transf <- factor(d$transf, levels=c("scaling", "cubic", "exp", "mixed"))
d$zeros  <- paste0(100 * d$zero_prop, "% zeros")
K <- 3

# One row per replicate and true group
groups <- do.call(rbind, lapply(1:K, function(k)
  data.frame(transf=d$transf, zeros=d$zeros, chordal=d[[paste0("chordal_", k)]],
             rank=d[[paste0("rank_", k)]])))

summary_table <- do.call(rbind, lapply(split(d, list(d$zero_prop, d$transf)), function(s){
  g <- groups[groups$transf == s$transf[1] & groups$zeros == s$zeros[1], ]
  data.frame(zero_prop=s$zero_prop[1], transf=s$transf[1],
             accuracy=mean(s$accuracy), accuracy_sd=sd(s$accuracy),
             rank_correct=mean(g$rank == 2),
             chordal_median_rank2=median(g$chordal[g$rank == 2]), chordal_rank2_sd=sd(g$chordal[g$rank == 2]),
             chordal_median=median(g$chordal), chordal_sd=sd(g$chordal))
}))
s <- summary_table[order(summary_table$zero_prop, summary_table$transf), ]
rows <- sprintf("%.0f\\%% & %s & %.3f (%.3f) & %.0f\\%% & %.3f (%.3f) & %.3f (%.3f) \\\\",
                100 * s$zero_prop, s$transf, s$accuracy, s$accuracy_sd, 100 * s$rank_correct,
                s$chordal_median_rank2, s$chordal_rank2_sd, s$chordal_median, s$chordal_sd)
rows <- append(rows, "\\midrule", after=sum(s$zero_prop == min(s$zero_prop)))
writeLines(c("\\begin{tabular}{@{}llcccc@{}}", "\\toprule",
             "Zeros & Transformation & Accuracy & Correct Rank & Distance (Correct) & Distance (All) \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"),
           "output/table9.tex")

separation <- rbind(data.frame(transf=d$transf, zeros=d$zeros, method="TLGC-PCA", retained=d$bh_tlgc / d$bh_true),
                    data.frame(transf=d$transf, zeros=d$zeros, method="Latent mixture", retained=d$bh_mix / d$bh_true))
separation$method <- factor(separation$method, levels=c("TLGC-PCA", "Latent mixture"))

fig <- ggplot(separation, aes(transf, retained, fill=method)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40") +
  geom_boxplot(outlier.size=0.6) +
  facet_wrap(~ zeros) +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=c("TLGC-PCA"="#E69F00", "Latent mixture"="#0072B2")) +
  labs(x="Transformation", y="Separation retained", fill=NULL) +
  theme_bw() + theme(legend.position="top")
ggsave("output/figure6.pdf", fig, width=6.5, height=3.4)
