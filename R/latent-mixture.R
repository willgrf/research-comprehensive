# Latent Gaussian mixture extension of TLGC-PCA (Section 4 of the report).
# Z_i | c_i = k ~ N(mu_k, Sigma_k), with each Sigma_k spiked, and transformations
# and thresholds shared by all groups. mu is a K x p matrix, Sigma a list of K
# p x p matrices.

source("R/tlgc-pca.R")


# ---- Mixture margins (Section 4.1) -----------------------------------------

mixture_cdf <- function(z, j, phi, mu, Sigma){
  s <- sqrt(sapply(Sigma, function(S) S[j, j]))
  sum(phi * pnorm((z - mu[, j]) / s))
}

mixture_quantile <- function(u, j, phi, mu, Sigma){
  # The u-quantile of G_j lies between the smallest and largest group u-quantiles
  if (u == 0) return(-Inf)
  s <- sqrt(sapply(Sigma, function(S) S[j, j]))
  group_quantiles <- mu[, j] + s * qnorm(u)
  if (diff(range(group_quantiles)) == 0) return(group_quantiles[1])
  uniroot(function(z) mixture_cdf(z, j, phi, mu, Sigma) - u,
          range(group_quantiles), tol=1e-10)$root
}

standardise_margins <- function(phi, mu, Sigma){
  # The same affine map in every group, so that each G_j has mean 0 and variance 1
  group_var <- t(sapply(Sigma, diag))
  m <- colSums(phi * mu)
  v <- colSums(phi * (group_var + mu^2)) - m^2
  list(mu    = sweep(sweep(mu, 2, m), 2, sqrt(v), "/"),
       Sigma = lapply(Sigma, function(S) S / outer(sqrt(v), sqrt(v))))
}

winsorised_ecdf <- function(X, winsor){
  U <- matrix(NA_real_, nrow=nrow(X), ncol=ncol(X))
  for (j in seq_len(ncol(X))){
    nonzero <- X[, j] != 0
    U[nonzero, j] <- pmin(pmax(ecdf(X[, j])(X[nonzero, j]), winsor), 1 - winsor)
  }
  U
}


# ---- Imputation under each sample's group (step 4) -------------------------

impute_zeros <- function(Z, is_trunc, del, labels, mu, Sigma, reps, burn){
  samples <- vector("list", nrow(Z))
  for (i in seq_len(nrow(Z))){
    censored <- which(is_trunc[i, ])
    observed <- which(!is_trunc[i, ])
    if (length(censored) == 0) next
    k <- labels[i]
    samples[[i]] <- gibbs_trunc_block(Sigma[[k]], observed, censored,
                                      Z[i, observed] - mu[k, observed],
                                      del[censored] - mu[k, censored], reps, burn)
    Z[i, censored] <- mu[k, censored] + colMeans(samples[[i]])
  }
  list(Z=Z, samples=samples)
}

conditional_covariances <- function(samples, is_trunc){
  # Row i holds the p x p matrix C_i of Equation 10: the covariance of the Gibbs
  # draws on the zero entries, 0 elsewhere
  p <- ncol(is_trunc)
  C <- matrix(0, nrow=nrow(is_trunc), ncol=p * p)
  for (i in seq_len(nrow(is_trunc))){
    censored <- which(is_trunc[i, ])
    if (length(censored) == 0) next
    C_i <- matrix(0, p, p)
    C_i[censored, censored] <- cov(samples[[i]])
    C[i, ] <- C_i
  }
  C
}


# ---- Spiked Gaussian mixture by EM (step 1) --------------------------------

e_step <- function(Z, phi, mu, Sigma){
  log_dens <- sapply(seq_along(phi), function(k)
    log(phi[k]) + dmnorm(Z, mu[k, ], Sigma[[k]], log=TRUE))
  log_max <- apply(log_dens, 1, max)
  log_total <- log_max + log(rowSums(exp(log_dens - log_max)))
  list(w=exp(log_dens - log_total), loglik=sum(log_total))
}

weighted_scatter <- function(Z, C, w){
  # N_k, mu_k and S_k of Equation 10
  N <- colSums(w)
  mu <- crossprod(w, Z) / N
  S <- lapply(seq_along(N), function(k){
    centred <- sweep(Z, 2, mu[k, ])
    (crossprod(centred * sqrt(w[, k])) + matrix(colSums(w[, k] * C), ncol(Z))) / N[k]
  })
  list(N=N, mu=mu, S=S)
}

spiked_covariance <- function(S, r){
  # V Lambda V' + sigma^2 (I - V V'), the maximum likelihood spiked covariance given S
  e <- eigen(S, symmetric=TRUE)
  V <- e$vectors[, 1:r, drop=FALSE]
  sigma2 <- mean(e$values[-(1:r)])
  tcrossprod(V * rep(sqrt(e$values[1:r] - sigma2), each=nrow(S))) + diag(sigma2, nrow(S))
}

select_ranks <- function(Z, C, w){
  sapply(weighted_scatter(Z, C, w)$S, function(S)
    igraph::dim_select(eigen(S, symmetric=TRUE, only.values=TRUE)$values))
}

m_step <- function(Z, C, w, r){
  scatter <- weighted_scatter(Z, C, w)
  list(phi=scatter$N / nrow(Z), mu=scatter$mu,
       Sigma=mapply(spiked_covariance, scatter$S, r, SIMPLIFY=FALSE))
}

mixture_em <- function(Z, C, w, r, tol=1e-5){
  # Iterates until the log-likelihood stops increasing (relative tolerance tol)
  loglik <- -Inf
  iterations <- 0
  repeat {
    par <- m_step(Z, C, w, r)
    e <- e_step(Z, par$phi, par$mu, par$Sigma)
    w <- e$w
    iterations <- iterations + 1
    if (e$loglik - loglik < tol * (1 + abs(e$loglik))) break
    loglik <- e$loglik
  }
  c(par, list(w=w, loglik=e$loglik, iterations=iterations))
}


# ---- Fitting procedure (Section 4.2) ---------------------------------------

fit_latent_mixture <- function(X, K, r=NULL, max_cycles=20, winsor=1/(2*nrow(X)),
                               reps=200, burn=100, init=NULL, verbose=FALSE){
  # r: one rank for every group, K group ranks, or NULL to select each group's
  # rank by profile likelihood once per cycle. init: starting labels (default:
  # k-means clusters of the TLGC-PCA latent vectors).
  is_trunc <- (X == 0)
  pi_hat <- colMeans(is_trunc)
  F_hat <- winsorised_ecdf(X, winsor)

  start <- latent_estimation(X, latent_corr(X), winsor, reps, burn)
  Z <- start$Z
  C <- conditional_covariances(start$samples, is_trunc)

  clusters <- if (is.null(init)) kmeans(Z, centers=K, nstart=25)$cluster else init
  w <- 1 * outer(clusters, seq_len(K), "==")

  labels <- NULL
  trace <- NULL
  for (t in seq_len(max_cycles) - 1){
    # Step 1: spiked mixture, with the ranks fixed during EM
    ranks <- if (is.null(r)) select_ranks(Z, C, w) else rep_len(r, K)
    gmm <- mixture_em(Z, C, w, ranks)
    w <- gmm$w

    # Step 2: standardise the margins and map thresholds and nonzero entries through G_j^{-1}
    std <- standardise_margins(gmm$phi, gmm$mu, gmm$Sigma)
    phi <- gmm$phi
    mu <- std$mu
    Sigma <- std$Sigma
    del <- sapply(seq_len(ncol(X)), function(j) mixture_quantile(pi_hat[j], j, phi, mu, Sigma))
    for (j in seq_len(ncol(X))){
      nonzero <- !is_trunc[, j]
      Z[nonzero, j] <- sapply(F_hat[nonzero, j], mixture_quantile,
                              j=j, phi=phi, mu=mu, Sigma=Sigma)
    }

    # Step 3: most likely group
    previous_labels <- labels
    labels <- max.col(w, ties.method="first")

    # Step 4: re-impute the zeros under each sample's group
    imputed <- impute_zeros(Z, is_trunc, del, labels, mu, Sigma, reps, burn)
    Z <- imputed$Z
    C <- conditional_covariances(imputed$samples, is_trunc)

    changed <- if (t == 0) NA else sum(labels != previous_labels)
    trace <- rbind(trace, data.frame(cycle=t, loglik=gmm$loglik, em_iterations=gmm$iterations,
                                     ranks=paste(ranks, collapse="/"),
                                     group_sizes=paste(tabulate(labels, K), collapse="/"),
                                     changed=changed))
    if (verbose) message(sprintf("Cycle %d: ranks %s, log-likelihood %.1f, group sizes %s, %s labels changed",
                                 t, trace$ranks[t + 1], gmm$loglik, trace$group_sizes[t + 1], changed))
    if (!is.na(changed) && changed == 0) break
  }
  converged <- !is.na(changed) && changed == 0

  # Step 2 leaves Sigma_k only approximately spiked, so the final parameters
  # come from one more spiked M-step
  final <- m_step(Z, C, w, ranks)
  phi <- final$phi
  mu <- final$mu
  Sigma <- final$Sigma

  # Scores U_i = V_k' (Z_i - mu_k) for the sample's group k; row i is NA beyond r_k.
  # scree: eigenvalues of the unconstrained scatter matrices the ranks are selected from.
  V <- lapply(seq_len(K), function(k)
    eigen(Sigma[[k]], symmetric=TRUE)$vectors[, 1:ranks[k], drop=FALSE])
  scree <- lapply(weighted_scatter(Z, C, w)$S, function(S)
    eigen(S, symmetric=TRUE, only.values=TRUE)$values)
  U <- matrix(NA_real_, nrow=nrow(X), ncol=max(ranks))
  for (i in seq_len(nrow(X))){
    k <- labels[i]
    U[i, 1:ranks[k]] <- crossprod(V[[k]], Z[i, ] - mu[k, ])
  }

  list(phi=phi, mu=mu, Sigma=Sigma, R=lapply(Sigma, cov2cor), V=V, r=ranks, scree=scree,
       del=del, labels=labels, responsibilities=w, Z=Z, U=U, trace=trace, converged=converged)
}
