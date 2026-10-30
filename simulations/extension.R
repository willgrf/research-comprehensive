# Simulation studies of the latent mixture extension.
#   --study=mixture          Section 4.3: recovery and group separation
#                            -> out_dir/{replicates.csv, components.csv}
#   --study=transformations  Appendix D.2: the same under each marginal transformation
#                            -> output/transformations.csv

suppressPackageStartupMessages(source("R/latent-mixture.R"))
source("simulations/common.R")


# ---- Section 4.3 -----------------------------------------------------------

mixture_replicate <- function(setting, seed){
  set.seed(seed)
  K <- setting$K
  phi <- as.numeric(strsplit(setting$phi, "/")[[1]])
  d <- simulate_mixture(setting$n, setting$p, K, setting$r, setting$sig, setting$delta,
                        phi, setting$zero_prop)
  Z_tlgc <- latent_estimation(d$X, latent_corr(d$X))$Z

  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(fit_latent_mixture(d$X, K), error=function(e) e)
  time_sec <- proc.time()[["elapsed"]] - start
  failed <- inherits(fit, "error")

  replicate_row <- data.frame(setting, seed=seed, zero_prop_realised=mean(d$X == 0),
                              time_sec=time_sec, error=if (failed) conditionMessage(fit) else NA,
                              bh_true=bhattacharyya(d$Z, d$c), bh_tlgc=bhattacharyya(Z_tlgc, d$c))
  if (failed) return(list(replicates=replicate_row, components=NULL))

  m <- match_groups(fit$labels, d$c, K)
  replicate_row <- cbind(replicate_row, data.frame(
    accuracy=mean(fit$labels == m[d$c]), cycles=nrow(fit$trace), converged=fit$converged,
    bh_mix=bhattacharyya(fit$Z, d$c)))
  components <- data.frame(setting, seed=seed, group=1:K,
                           chordal=mapply(chordal_dist, fit$V[m], d$V), rank_hat=fit$r[m])
  list(replicates=replicate_row, components=components)
}

run_mixture_sims <- function(n=c(300, 600, 1200), zero_prop=c(0, 0.2, 0.4), delta=c(0.5, 1),
                             p=10, K=3, r=2, sig=0.5, phi="0.25/0.35/0.4",
                             num_reps=50, cores=1, seed=1, out_dir="output/mixture"){
  settings <- expand.grid(n=n, zero_prop=zero_prop, delta=delta, p=p, K=K, r=r, sig=sig,
                          phi=phi, stringsAsFactors=FALSE)
  tasks <- expand.grid(rep=seq_len(num_reps), setting=seq_len(nrow(settings)))
  tasks$seed <- seed * 1e6 + tasks$setting * 1e3 + tasks$rep

  # Slowest tasks (largest n) first
  order_run <- order(-settings$n[tasks$setting], -settings$zero_prop[tasks$setting])
  res <- parallel::mclapply(order_run, function(i){
    setting <- settings[tasks$setting[i], ]
    out <- mixture_replicate(cbind(setting_id=tasks$setting[i], rep=tasks$rep[i], setting),
                             tasks$seed[i])
    cat(sprintf("n=%-5d zeros=%.1f delta=%.1f rep %3d  %5.0fs %s\n", setting$n, setting$zero_prop,
                setting$delta, tasks$rep[i], out$replicates$time_sec,
                if (is.na(out$replicates$error)) "" else paste("ERROR:", out$replicates$error)))
    out
  }, mc.cores=cores, mc.preschedule=FALSE)

  crashed <- vapply(res, inherits, logical(1), "try-error")
  if (any(crashed)) warning(sum(crashed), " tasks crashed: ", paste(unique(unlist(res[crashed])), collapse="; "))
  res <- res[!crashed]

  # Failed fits have no fit columns; fill them with NA before binding
  replicates <- lapply(res, `[[`, "replicates")
  all_columns <- unique(unlist(lapply(replicates, names)))
  replicates <- do.call(rbind, lapply(replicates, function(x){ x[setdiff(all_columns, names(x))] <- NA; x[all_columns] }))
  replicates <- replicates[order(replicates$setting_id, replicates$rep), ]
  components <- do.call(rbind, lapply(res, `[[`, "components"))
  components <- components[order(components$setting_id, components$rep, components$group), ]

  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(replicates, file.path(out_dir, "replicates.csv"), row.names=FALSE)
  write.csv(components, file.path(out_dir, "components.csv"), row.names=FALSE)
}


# ---- Appendix D.2 ----------------------------------------------------------
# The same latent data (seed depends on the zero proportion and replicate only)
# observed through each transformation and rounded to counts.

observe <- function(Z, del, transf){
  # round(g(Z_ij) - g(b_j)) above the threshold, 0 below, with b_j the threshold
  # (or the column minimum when there is none)
  n <- nrow(Z)
  base <- pmin(del, apply(Z, 2, min))
  X <- round(transform_latent(Z, transf) - rep(transform_latent(base, transf), each=n))
  X * (Z > rep(del, each=n))
}

transformation_replicate <- function(seed, transf, zero_prop, n, p, K, r, sig, delta, phi){
  set.seed(seed)
  d <- simulate_mixture(n, p, K, r, sig, delta, phi, zero_prop)
  X <- observe(d$Z, d$del, transf)
  Z_tlgc <- latent_estimation(X, latent_corr(X))$Z

  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(fit_latent_mixture(X, K), error=function(e) e)
  time_sec <- proc.time()[["elapsed"]] - start
  failed <- inherits(fit, "error")

  # rank_k / chordal_k: selected rank and chordal distance for true group k
  metrics <- list(accuracy=NA, converged=NA, cycles=NA, bh_mix=NA, rank_correct=NA)
  metrics[paste0("chordal_", 1:K)] <- NA
  metrics[paste0("rank_", 1:K)] <- NA
  if (!failed){
    m <- match_groups(fit$labels, d$c, K)
    metrics$accuracy     <- mean(fit$labels == m[d$c])
    metrics$converged    <- fit$converged
    metrics$cycles       <- nrow(fit$trace)
    metrics$bh_mix       <- bhattacharyya(fit$Z, d$c)
    metrics$rank_correct <- mean(fit$r[m] == r)
    metrics[paste0("chordal_", 1:K)] <- mapply(chordal_dist, fit$V[m], d$V)
    metrics[paste0("rank_", 1:K)]    <- fit$r[m]
  }

  data.frame(transf=transf, zero_prop=zero_prop, n=n, rank_true=r, seed=seed,
             zero_prop_realised=mean(X == 0), time_sec=time_sec,
             error=if (failed) conditionMessage(fit) else NA,
             bh_true=bhattacharyya(d$Z, d$c), bh_tlgc=bhattacharyya(Z_tlgc, d$c),
             metrics)
}

run_transformation_sims <- function(transf=c("scaling", "cubic", "exp", "mixed"),
                                    zero_prop=c(0.2, 0.4), n=600, p=10, K=3, r=2,
                                    sig=0.5, delta=1, phi=c(0.25, 0.35, 0.4),
                                    num_reps=20, cores=1, seed=1,
                                    file="output/transformations.csv"){
  # Each setting is appended to file as soon as it finishes
  if (file.exists(file)) stop(file, " already exists; move it or choose another file")
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  settings <- expand.grid(transf=transf, zero_prop=zero_prop, stringsAsFactors=FALSE)

  for (i in seq_len(nrow(settings))){
    setting <- settings[i, ]
    seeds <- seed * 1e6 + match(setting$zero_prop, zero_prop) * 1e3 + seq_len(num_reps)
    rows <- parallel::mclapply(seeds, transformation_replicate, transf=setting$transf,
                               zero_prop=setting$zero_prop, n=n, p=p, K=K, r=r,
                               sig=sig, delta=delta, phi=phi,
                               mc.cores=cores, mc.preschedule=FALSE)
    crashed <- vapply(rows, inherits, logical(1), "try-error")
    if (any(crashed)) warning(sum(crashed), " replicates crashed: ", unique(unlist(rows[crashed])))
    rows <- do.call(rbind, rows[!crashed])

    write.table(rows, file, sep=",", row.names=FALSE,
                col.names=!file.exists(file), append=file.exists(file))
    cat(sprintf("%-8s zeros = %2.0f%%  accuracy = %.3f\n",
                setting$transf, 100 * setting$zero_prop, mean(rows$accuracy, na.rm=TRUE)))
  }
}


if (sys.nframe() == 0L){
  args <- parse_command_line(commandArgs(trailingOnly=TRUE))
  study <- if (is.null(args$study)) "mixture" else args$study
  args$study <- NULL
  do.call(switch(study, mixture=run_mixture_sims, transformations=run_transformation_sims), args)
}
