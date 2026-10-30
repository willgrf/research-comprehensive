# LaTeX tables for Section 3 and Appendix D.1 (booktabs):
#   Tables 1, 6-8 (chordal distances), Table 2 (rank selection) from output/replication.csv
#   Table 3 (Poisson model) from output/poisson.csv
#   Rscript simulations/replication-results.R   -> output/replication-tables.tex
# Bold marks the best mean and any method within one SD (of the best method) of it.

methods <- c("Truncated copula PCA", "Standard PCA (unscaled)",
             "Standard PCA (scaled)", "PoissonPCA", "ZIPFA")
method_header <- "TLGC-PCA & PCA (unscaled) & PCA (scaled) & PoissonPCA & ZIPFA"
transf_names <- c(scaling="Scaling", cubic="Cubic", exp="Exponential", mixed="Mixed")

fmt <- function(x) sprintf("%.3f", x)

table_env <- function(position, header, rows, caption, label){
  c(sprintf("\\begin{table}[%s]", position), "\\centering", "\\small",
    "\\begin{tabular}{@{}lccccc@{}}", "\\toprule", header, "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}",
    paste0("\\caption{", caption, "}"), sprintf("\\label{%s}", label), "\\end{table}", "")
}

distance_cells <- function(s, measure){
  means <- sapply(methods, function(m) mean(s[s$method == m, measure], na.rm=TRUE))
  sds   <- sapply(methods, function(m) sd(s[s$method == m, measure], na.rm=TRUE))
  best <- which.min(means)
  cells <- sprintf("%s (%s)", fmt(means), fmt(sds))
  list(cells=cells, bold=means <= means[best] + sds[best])
}

loadings_and_scores <- function(row_fn){
  c("\\multicolumn{6}{@{}l}{\\textit{Loadings}} \\\\", row_fn("dist_V"), "\\midrule",
    "\\multicolumn{6}{@{}l}{\\textit{Scores}} \\\\", row_fn("dist_U"))
}


# ---- Tables 1, 6-8 and 2 ---------------------------------------------------

d <- read.csv("output/replication.csv", stringsAsFactors=FALSE)
zero_levels <- sort(unique(d$trunc_prop))
n_reps <- max(d$rep)

# ZIPFA fails on some replicates at high zero proportions; its mean is over the
# successful fits and a superscript gives their number
zipfa_failures <- function(transf, level)
  sum(d$transf == transf & d$trunc_prop == level & d$method == "ZIPFA" & is.na(d$dist_V))

distance_table <- function(transf){
  rows <- loadings_and_scores(function(measure) sapply(zero_levels, function(level){
    x <- distance_cells(d[d$transf == transf & d$trunc_prop == level, ], measure)
    n_failed <- zipfa_failures(transf, level)
    if (n_failed > 0) x$cells[5] <- sprintf("%s$^{%d}$", x$cells[5], n_reps - n_failed)
    x$cells[x$bold] <- sprintf("\\textbf{%s}", x$cells[x$bold])
    sprintf("%g\\%% & %s \\\\", 100*level, paste(x$cells, collapse=" & "))
  }))
  any_failed <- any(sapply(zero_levels, function(level) zipfa_failures(transf, level)) > 0)
  caption <- paste0(
    "Chordal distance between the true and estimated loading and score\n",
    "subspaces under the ", tolower(transf_names[[transf]]),
    " transformation ($n=200$, $p=20$, $r=4$,\n",
    "$\\sigma=0.1$). Mean (SD) over $", n_reps, "$ replications. Bold marks the best method\n",
    "and any method within one standard deviation of it.",
    if (any_failed) paste0(" A superscript on a ZIPFA entry gives the number of\n",
                           "replications in which ZIPFA returned a fit; the mean and SD are\n",
                           "over those replications only.") else "")
  table_env("ht", paste("Zeros &", method_header, "\\\\"), rows, caption, paste0("tab:sim-", transf))
}

rank_table <- function(){
  # A failed fit counts as not selecting the true rank
  rows <- sapply(names(transf_names), function(transf){
    pct <- sapply(methods, function(m){
      s <- d[d$transf == transf & d$method == m, ]
      round(100*mean(!is.na(s$rank_hat) & s$rank_hat == s$r))
    })
    sprintf("%-11s & %s \\\\", transf_names[[transf]], paste(pct, collapse=" & "))
  })
  level_labels <- paste0(100*zero_levels, "\\%")
  caption <- paste0(
    "Percentage of replications in which the profile likelihood criterion\n",
    "selects the true rank $r=4$, pooled over the ", length(zero_levels),
    " zero proportions from $", level_labels[1], "$\nto $", level_labels[length(level_labels)],
    "$ ($", length(zero_levels)*n_reps, "$ replications per entry). A replication in which\n",
    "ZIPFA fails to return a fit counts as not selecting the true rank.")
  table_env("t", paste("Transformation &", method_header, "\\\\"), rows, caption, "tab:sim-rank")
}


# ---- Table 3 ---------------------------------------------------------------

poisson_table <- function(){
  dp <- read.csv("output/poisson.csv", stringsAsFactors=FALSE)
  rows <- loadings_and_scores(function(measure) sapply(sort(unique(dp$c)), function(c){
    x <- distance_cells(dp[dp$c == c, ], measure)
    x$cells[x$bold] <- sprintf("\\textbf{%s}", x$cells[x$bold])
    sprintf("%g (%.0f\\%%) & %s \\\\", c, 100*mean(dp$zero_prop[dp$c == c]), paste(x$cells, collapse=" & "))
  }))
  caption <- paste0(
    "Chordal distance between the true and estimated loading and score subspaces\n",
    "under the Poisson measurement model ($n=200$, $p=20$, $r=4$). Loadings are compared\n",
    "with the column space of $W$ and scores with the centred $\\eta$. Rows give $c$ and,\n",
    "in parentheses, the average proportion of zeros. Mean (SD) over $", max(dp$rep), "$\n",
    "replications. Bold marks the best method and any method within one standard deviation of it.")
  table_env("ht", paste("$c$ (zeros) &", method_header, "\\\\"), rows, caption, "tab:sim-poisson")
}


out <- c("% Generated by simulations/replication-results.R; requires booktabs", "",
         "% Table 1", distance_table("mixed"),
         "% Table 2", rank_table(),
         "% Table 3", poisson_table(),
         "% Tables 6-8", distance_table("scaling"), distance_table("cubic"), distance_table("exp"))
writeLines(out, "output/replication-tables.tex")
