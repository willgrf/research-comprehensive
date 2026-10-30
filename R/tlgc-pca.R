# Truncated latent Gaussian copula PCA (Section 2 of the report).
#   latent_corr()       rank-based estimate of the latent correlation matrix (Section 2.3)
#   latent_estimation() latent vectors, with zeros imputed by Gibbs sampling (Section 2.5)
#   TLGC_PCA()          loadings and scores

library(mnormt)

MVN_EPS <- 1e-10


# ---- Bridge functions (Proposition 2.2) ------------------------------------

bridge_CC <- function(r) 2 * asin(r) / pi

bridge_CT <- function(r, del1){
  # Variable 1 truncated at del1, variable 2 continuous
  Sig_A <- matrix(c(         1, -1/sqrt(2), -r/sqrt(2),
                    -1/sqrt(2),          1,          r,
                    -r/sqrt(2),          r,          1), nrow=3)
  A <- 2*sadmvn(lower=c(del1, 0, 0), upper=c(Inf, Inf, Inf),
                mean=c(0, 0, 0), varcov=Sig_A, abseps=MVN_EPS) -
       2*sadmvn(lower=c(del1, 0, -Inf), upper=c(Inf, Inf, 0),
                mean=c(0, 0, 0), varcov=Sig_A, abseps=MVN_EPS)

  Sig_B <- matrix(c(        1,          0,  r/sqrt(2),
                            0,          1, -r/sqrt(2),
                    r/sqrt(2), -r/sqrt(2),          1), nrow=3)
  B <- 2*pnorm(del1)*(1 - pnorm(del1)) -
       4*sadmvn(lower=c(del1, -Inf, -Inf), upper=c(Inf, del1, 0),
                mean=c(0, 0, 0), varcov=Sig_B, abseps=MVN_EPS)

  as.numeric(A + B)
}

bridge_TT <- function(r, del1, del2){
  Sig_A <- matrix(c(        1,          0,  1/sqrt(2), -r/sqrt(2),
                            0,          1, -r/sqrt(2),  1/sqrt(2),
                    1/sqrt(2), -r/sqrt(2),          1,         -r,
                   -r/sqrt(2),  1/sqrt(2),         -r,          1), nrow=4)
  Sig_B <- matrix(c(        1,          r,  1/sqrt(2),  r/sqrt(2),
                            r,          1,  r/sqrt(2),  1/sqrt(2),
                    1/sqrt(2),  r/sqrt(2),          1,          r,
                    r/sqrt(2),  1/sqrt(2),          r,          1), nrow=4)
  lower <- rep(-Inf, 4)
  upper <- c(-del1, -del2, 0, 0)
  as.numeric(
    -2*sadmvn(lower=lower, upper=upper, mean=rep(0, 4), varcov=Sig_A, abseps=MVN_EPS) +
     2*sadmvn(lower=lower, upper=upper, mean=rep(0, 4), varcov=Sig_B, abseps=MVN_EPS))
}


# ---- Latent correlation ----------------------------------------------------

solve_bridge <- function(bridge, tau, cuts, eps=1e-7, tol=1e-10){
  # Inverts the bridge by root finding on [-1 + eps, 1 - eps], where the bridge
  # covariance matrices are non-singular; tau is clipped to the attainable range
  rmax <- 1 - eps
  f <- function(r) do.call(bridge, c(list(r), as.list(unname(cuts))))
  tau <- min(max(tau, f(-rmax)), f(rmax))
  uniroot(function(r) f(r) - tau, c(-rmax, rmax), tol=tol)$root
}

latent_pair_corr <- function(tau, truncated, del, eps, tol){
  if (!any(truncated)) sin(pi * tau / 2)
  else if (all(truncated)) solve_bridge(bridge_TT, tau, del, eps, tol)
  else solve_bridge(bridge_CT, tau, del[truncated], eps, tol)
}

kendall_tau_a <- function(X){
  # cor() returns tau-b; the bridges are derived for tau-a, which also counts tied pairs
  n_pairs <- nrow(X) * (nrow(X) - 1) / 2
  tied_pairs <- apply(X, 2, function(x){
    run_lengths <- rle(sort(x))$lengths
    sum(run_lengths * (run_lengths - 1) / 2)
  })
  untied_pairs <- sqrt(n_pairs - tied_pairs)
  cor(X, method="kendall") * outer(untied_pairs, untied_pairs) / n_pairs
}

project_psd <- function(R, ridge){
  # Projection onto the PSD cone, then shrinkage towards the identity (Equation 5)
  e <- eigen(R, symmetric=TRUE)
  R <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  (1 - ridge)*R + diag(x=ridge, nrow=nrow(R))
}

latent_corr <- function(X, eps=1e-7, tol=1e-10, ridge=1e-5){
  # A column is truncated if it contains a zero, with threshold qnorm(P(X = 0))
  Tau <- kendall_tau_a(X)
  truncated <- apply(X == 0, 2, any)
  del <- qnorm(colMeans(X == 0))

  p <- ncol(X)
  R <- diag(p)
  for (j in seq_len(p - 1)){
    for (k in (j + 1):p){
      pair <- c(j, k)
      R[j, k] <- R[k, j] <- latent_pair_corr(Tau[j, k], truncated[pair], del[pair], eps, tol)
    }
  }
  project_psd(R, ridge)
}


# ---- Latent vectors --------------------------------------------------------

rtnorm_upper <- function(mu, sd, upper){
  # One draw from N(mu, sd^2) truncated above at upper, by inverse CDF on the log scale
  a <- (upper - mu) / sd
  mu + sd * qnorm(log(runif(1)) + pnorm(a, log.p=TRUE), log.p=TRUE)
}

gibbs_trunc_block <- function(R, observed, censored, z_observed, del_censored,
                              reps=100, burn=50, thin=1){
  # Gibbs draws (reps x length(censored)) of a row's censored latent entries,
  # given its observed entries, under N(0, R)
  k <- length(censored)
  if (k == 0) return(matrix(numeric(0), nrow=reps, ncol=0))

  if (length(observed) > 0){
    A  <- R[censored, observed, drop=FALSE] %*% solve(R[observed, observed, drop=FALSE])
    mu <- as.vector(A %*% z_observed)
    S  <- R[censored, censored, drop=FALSE] - A %*% R[observed, censored, drop=FALSE]
  } else {
    mu <- rep(0, k)
    S  <- R[censored, censored, drop=FALSE]
  }
  S <- (S + t(S)) / 2
  Q <- solve(S)
  sd_j <- 1 / sqrt(diag(Q))

  z <- pmin(mu, del_censored - 1e-6)
  draws <- matrix(NA_real_, nrow=reps, ncol=k, dimnames=list(NULL, colnames(R)[censored]))
  for (it in seq_len(burn + reps * thin)){
    for (j in seq_len(k)){
      mu_j <- mu[j] - as.numeric(Q[j, -j, drop=FALSE] %*% (z[-j] - mu[-j])) / Q[j, j]
      z[j] <- rtnorm_upper(mu_j, sd_j[j], del_censored[j])
    }
    if (it > burn && (it - burn) %% thin == 0) draws[(it - burn) %/% thin, ] <- z
  }
  draws
}

nonzero_normal_scores <- function(X, winsor){
  # qnorm of the winsorised empirical CDF for the nonzero entries; zeros are left NA
  Z <- matrix(NA_real_, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))
  for (j in seq_len(ncol(X))){
    nonzero <- X[, j] != 0
    Z[nonzero, j] <- qnorm(pmin(ecdf(X[, j])(X[nonzero, j]), 1 - winsor))
  }
  Z
}

latent_estimation <- function(X, R=latent_corr(X), winsor=1/(2*nrow(X)),
                              reps=200, burn=100, thin=1){
  is_trunc <- (X == 0)
  del <- qnorm(colMeans(is_trunc))
  Z <- nonzero_normal_scores(X, winsor)

  samples <- vector("list", nrow(X))
  for (i in seq_len(nrow(X))){
    censored <- which(is_trunc[i, ])
    observed <- which(!is_trunc[i, ])
    samples[[i]] <- gibbs_trunc_block(R, observed, censored, Z[i, observed], del[censored],
                                      reps, burn, thin)
    Z[i, censored] <- colMeans(samples[[i]])
  }
  list(Z=Z, del=del, is_trunc=is_trunc, samples=samples, R=R)
}


# ---- TLGC-PCA --------------------------------------------------------------

TLGC_PCA <- function(X, r=FALSE){
  # r = FALSE selects the rank by profile likelihood on the eigenvalues of R
  R <- latent_corr(X)
  e <- eigen(R, symmetric=TRUE)
  scree <- e$values
  if (!r) r <- igraph::dim_select(scree)
  V <- e$vectors[, 1:r]
  Z <- latent_estimation(X, R)$Z
  list(R=R, V=V, U=Z %*% V, Z=Z, r=r, scree=scree)
}
