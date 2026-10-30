# Helpers shared by the simulation scripts.

chordal_dist <- function(A, B){
  # Equation 6, through the projections onto the column spaces of A and B
  P_A <- A %*% MASS::ginv(A)
  P_B <- B %*% MASS::ginv(B)
  norm(P_A - P_B, type="F") / sqrt(2)
}

transform_latent <- function(Z, transf){
  # The marginal transformations of Section 3.1
  switch(transf,
         scaling = 10*Z,
         cubic   = 10*Z^3,
         exp     = 10*exp(Z),
         mixed   = 10 + 10*Z + 10*Z^3 + 10*exp(Z),
         Z)
}

parse_command_line <- function(args){
  # "--name=value" or "--name value"; comma-separated values become vectors
  args <- trimws(strsplit(paste(args, collapse=" "), "\\s*--")[[1]])
  args <- args[nzchar(args)]
  options <- list()
  for (name_value in strsplit(args, "[= ]+")){
    options[[name_value[1]]] <- type.convert(strsplit(name_value[2], ",")[[1]], as.is=TRUE)
  }
  options
}

l_ecuyer_streams <- function(n_tasks, seed){
  # One random number stream per task, so results do not depend on the number of cores
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- list(.Random.seed)
  for (i in seq_len(n_tasks - 1)) streams[[i + 1]] <- parallel::nextRNGStream(streams[[i]])
  streams
}


# ---- Latent mixture data (Section 4.3) -------------------------------------

spiked_covariance_truth <- function(s2, r, sig){
  # W W' + sig^2 I with diagonal s2
  W <- matrix(rnorm(length(s2)*r), ncol=r)
  W <- sqrt(s2 - sig^2) * W / sqrt(rowSums(W^2))
  tcrossprod(W) + sig^2 * diag(length(s2))
}

simulate_mixture <- function(n, p, K, r, sig, delta, phi, zero_prop, s_range=c(0.75, 1.25)){
  # Group means ~ N(0, delta^2), centred. The within-group variances are rescaled
  # within each variable so every variable has the same total variance v, and
  # dividing by v gives each margin G_j mean 0 and variance 1. X = exp(Z), truncated
  # at Delta_j = G_j^{-1}(zero_prop).
  mu <- matrix(rnorm(K*p, sd=delta), nrow=K, ncol=p)
  mu <- sweep(mu, 2, colSums(phi * mu))
  between <- colSums(phi * mu^2)
  v <- 1 + max(between)
  s2 <- matrix(runif(K*p, s_range[1], s_range[2])^2, nrow=K, ncol=p)
  s2 <- sweep(s2, 2, (v - between) / colSums(phi * s2), "*")
  stopifnot(all(s2 > sig^2))
  Sigma <- lapply(1:K, function(k) spiked_covariance_truth(s2[k, ], r, sig) / v)
  mu <- mu / sqrt(v)
  del <- sapply(1:p, function(j) mixture_quantile(zero_prop, j, phi, mu, Sigma))

  c <- sample(1:K, n, replace=TRUE, prob=phi)
  Z <- matrix(NA_real_, nrow=n, ncol=p)
  for (k in 1:K){
    in_k <- which(c == k)
    E <- matrix(rnorm(length(in_k)*p), ncol=p)
    Z[in_k, ] <- sweep(E %*% chol(Sigma[[k]]), 2, mu[k, ], "+")
  }
  X <- exp(Z) * (Z > rep(del, each=n))

  list(X=X, Z=Z, c=c, phi=phi, mu=mu, Sigma=Sigma, del=del,
       V=lapply(Sigma, function(S) eigen(S, symmetric=TRUE)$vectors[, 1:r]))
}

match_groups <- function(labels, c, K){
  # m[k] is the fitted group matched to true group k, maximising agreement
  perms <- as.matrix(expand.grid(rep(list(1:K), K)))
  perms <- perms[apply(perms, 1, function(m) length(unique(m)) == K), , drop=FALSE]
  agreement <- apply(perms, 1, function(m) mean(labels == m[c]))
  perms[which.max(agreement), ]
}

bhattacharyya <- function(Y, c){
  # Gaussian Bhattacharyya distance between the true groups, averaged over pairs
  groups <- lapply(seq_len(max(c)), function(k)
    list(mean=colMeans(Y[c == k, , drop=FALSE]), cov=cov(Y[c == k, , drop=FALSE])))
  log_det <- function(M) as.numeric(determinant(M)$modulus)
  mean(apply(combn(length(groups), 2), 2, function(ab){
    a <- groups[[ab[1]]]; b <- groups[[ab[2]]]
    S <- (a$cov + b$cov) / 2
    d <- a$mean - b$mean
    sum(d * solve(S, d)) / 8 + (log_det(S) - (log_det(a$cov) + log_det(b$cov)) / 2) / 2
  }))
}
