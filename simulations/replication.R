# Replication study (Section 3): TLGC-PCA against standard PCA, PoissonPCA and ZIPFA.
#   --study=tlgc     data from the TLGC model (Section 3.2)  -> output/replication.csv
#   --study=poisson  data from the Poisson model (Section 3.3) -> output/poisson.csv

source("R/tlgc-pca.R")
source("simulations/common.R")


# ---- Methods ---------------------------------------------------------------
# Each returns loadings V (p x r), scores U (n x r) and the scree values used for
# rank selection (on the square-root scale).

fit_tlgc_pca <- function(X, r){
  fit <- TLGC_PCA(X, r)
  list(V=fit$V, U=fit$U, scree=sqrt(fit$scree))
}

fit_pca <- function(X, r, scale){
  pc <- prcomp(X, center=TRUE, scale.=scale)
  list(V=pc$rotation[, 1:r], U=pc$x[, 1:r], scree=pc$sdev)
}

fit_poisson_pca <- function(X, r){
  # Scores from the centred X, as for the other methods
  fit <- PoissonPCA::Poisson_Corrected_PCA(X, k=r)
  list(V=fit$loadings, U=scale(X, center=TRUE, scale=FALSE) %*% fit$loadings, scree=fit$sdev)
}

fit_zipfa <- function(X, r){
  # Keeps every column (cut=1) and skips the row adjustment (Madj=FALSE), whose
  # defaults fail on heavily zero-inflated data; failed fits are retried from
  # other starting values
  default_args <- list(X=X, k=r, tau=0.1, cut=1, Madj=FALSE, display=FALSE)
  attempts <- list(list(), list(initialtau="stable"), list(tau=1))
  for (a in seq_along(attempts)){
    fit <- tryCatch(do.call(ZIPFA::ZIPFA, modifyList(default_args, attempts[[a]])),
                    error=function(e) e)
    if (!inherits(fit, "error")) break
    if (a == length(attempts)) stop(fit)
    warning("ZIPFA retry ", a, " after: ", trimws(conditionMessage(fit)))
  }
  U <- fit$Ufit[[fit$itr]]
  V <- fit$Vfit[[fit$itr]]
  list(V=V, U=U, scree=svd(U %*% t(V))$d)
}

sim_methods <- list(
  "Truncated copula PCA"    = fit_tlgc_pca,
  "Standard PCA (unscaled)" = function(X, r) fit_pca(X, r, scale=FALSE),
  "Standard PCA (scaled)"   = function(X, r) fit_pca(X, r, scale=TRUE),
  "PoissonPCA"              = fit_poisson_pca,
  "ZIPFA"                   = fit_zipfa
)

evaluate_method <- function(method, data, r){
  # Chordal distances to the true loadings and scores and the profile likelihood
  # rank. Errors give NA results; warnings are collected.
  start <- proc.time()[["elapsed"]]
  warnings <- character(0)
  result <- withCallingHandlers(
    tryCatch({
      fit <- sim_methods[[method]](data$X, r)
      list(dist_V=chordal_dist(fit$V, data$V), dist_U=chordal_dist(fit$U, data$U),
           rank_hat=igraph::dim_select(sort(fit$scree, decreasing=TRUE)), error=NA_character_)
    }, error=function(e){
      list(dist_V=NA_real_, dist_U=NA_real_, rank_hat=NA_real_, error=conditionMessage(e))
    }),
    warning=function(w){
      warnings <<- c(warnings, trimws(conditionMessage(w)))
      invokeRestart("muffleWarning")
    })

  data.frame(method=method, dist_V=result$dist_V, dist_U=result$dist_U,
             rank_hat=result$rank_hat, time_sec=proc.time()[["elapsed"]] - start,
             error=result$error,
             warning=if (length(warnings)) paste(unique(warnings), collapse="; ") else NA)
}

evaluate_all <- function(data, settings, r){
  do.call(rbind, lapply(names(sim_methods), function(m) cbind(settings, evaluate_method(m, data, r))))
}

stop_if_failed <- function(res){
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) stop(sum(failed), " tasks failed:\n", paste(unique(unlist(res[failed])), collapse="\n"))
}


# ---- TLGC model (Section 3.2) ----------------------------------------------

simulate_tlgc_data <- function(n, p, r, transf, trunc_prop, sig){
  # Z = U W' + E with unit-diagonal W W' + sig^2 I; transform, shift the minimum
  # to zero, set the smallest trunc_prop of entries to zero and round
  W <- matrix(rnorm(p*r), nrow=p, ncol=r)
  W <- sqrt(1 - sig^2)*W/sqrt(diag(W %*% t(W)))
  U <- matrix(rnorm(n*r), nrow=n, ncol=r)
  E <- matrix(rnorm(n*p, sd=sig), nrow=n, ncol=p)
  X_star <- transform_latent(U %*% t(W) + E, transf)
  X_star <- X_star - min(X_star)
  X <- round(X_star*(X_star >= as.numeric(quantile(X_star, trunc_prop))))
  V <- eigen(W %*% t(W) + sig^2*diag(p), symmetric=TRUE)$vectors[, 1:r]
  list(X=X, U=U, V=V)
}

run_tlgc_sims <- function(n=200, p=20, r=4, sig=0.1,
                          transf=c("scaling", "cubic", "exp", "mixed"),
                          trunc_prop=c(0, .05, .1, .15, .2, .3, .4, .5, .7, .9),
                          num_reps=100, seed=1, cores=1, run_levels=trunc_prop,
                          file="output/replication.csv"){
  # One stream per arm (transformation x zero proportion). run_levels runs a
  # subset of the zero proportions with the same draws as the full grid.
  arms <- expand.grid(trunc_prop=trunc_prop, transf=transf, stringsAsFactors=FALSE)
  streams <- l_ecuyer_streams(nrow(arms), seed)

  res <- parallel::mclapply(which(arms$trunc_prop %in% run_levels), function(i){
    assign(".Random.seed", streams[[i]], envir=.GlobalEnv)
    cat(sprintf("%-8s pi0 = %3.0f%%\n", arms$transf[i], 100*arms$trunc_prop[i]))
    do.call(rbind, lapply(seq_len(num_reps), function(rep){
      data <- simulate_tlgc_data(n, p, r, arms$transf[i], arms$trunc_prop[i], sig)
      settings <- data.frame(n=n, p=p, r=r, transf=arms$transf[i],
                             trunc_prop=arms$trunc_prop[i], sig=sig)
      cbind(evaluate_all(data, settings, r), rep=rep)
    }))
  }, mc.cores=cores, mc.preschedule=FALSE)
  stop_if_failed(res)

  res <- do.call(rbind, res)
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  write.csv(res, file, row.names=FALSE)
  invisible(res)
}


# ---- Poisson model (Section 3.3) -------------------------------------------

simulate_poisson_data <- function(n, p, r, c){
  W <- matrix(rgamma(p*r, shape=1, rate=1), nrow=p, ncol=r)
  eta <- matrix(rgamma(n*r, shape=2, rate=1), nrow=n, ncol=r)
  eps <- matrix(rgamma(n*p, shape=0.5, rate=1), nrow=n, ncol=p)
  X <- matrix(rpois(n*p, c * (eta %*% t(W) + eps)), nrow=n, ncol=p)
  list(X=X, V=W, U=scale(eta, center=TRUE, scale=FALSE))
}

run_poisson_sims <- function(n=200, p=20, r=4, c=base::c(0.1, 0.25, 0.5, 1), num_reps=100,
                             seed=1, cores=1, file="output/poisson.csv"){
  # One stream per replicate; the smallest c (slowest) runs first
  tasks <- expand.grid(rep=seq_len(num_reps), c=c)
  streams <- l_ecuyer_streams(nrow(tasks), seed)

  res <- parallel::mclapply(order(tasks$c), function(i){
    assign(".Random.seed", streams[[i]], envir=.GlobalEnv)
    data <- simulate_poisson_data(n, p, r, tasks$c[i])
    settings <- data.frame(n=n, p=p, r=r, c=tasks$c[i], zero_prop=mean(data$X == 0))
    out <- cbind(evaluate_all(data, settings, r), rep=tasks$rep[i])
    cat(sprintf("c = %-5g rep %3d  %5.0fs\n", tasks$c[i], tasks$rep[i], sum(out$time_sec)))
    out
  }, mc.cores=cores, mc.preschedule=FALSE)
  stop_if_failed(res)

  res <- do.call(rbind, res)
  res <- res[order(res$c, res$rep, match(res$method, names(sim_methods))), ]
  dir.create(dirname(file), recursive=TRUE, showWarnings=FALSE)
  write.csv(res, file, row.names=FALSE)
  invisible(res)
}


if (sys.nframe() == 0L){
  args <- parse_command_line(commandArgs(trailingOnly=TRUE))
  study <- if (is.null(args$study)) "tlgc" else args$study
  args$study <- NULL
  do.call(switch(study, tlgc=run_tlgc_sims, poisson=run_poisson_sims), args)
}
