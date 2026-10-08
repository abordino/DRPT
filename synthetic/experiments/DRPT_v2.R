## ============================================================
## DRPT v2
## Equation (12) from arXiv.v2
##
## Uses DRPT::starSampler() for the permutation sampler,
## but computes the statistic in R.
## ============================================================

.eval_r_rows = function(Z, r) {
  
  Z = as.matrix(Z)
  N = nrow(Z)
  d = ncol(Z)
  
  args = lapply(
    seq_len(d),
    function(j) Z[, j]
  )
  
  out = tryCatch(
    do.call(r, args),
    error = function(e) NULL
  )
  
  if (is.null(out) || length(out) != N) {
    
    out = vapply(
      seq_len(N),
      function(i) {
        as.numeric(
          do.call(
            r,
            as.list(Z[i, ])
          )
        )
      },
      numeric(1)
    )
  }
  
  out = as.numeric(out)
  
  if (any(!is.finite(out))) {
    stop("r(Z) contains non-finite values.")
  }
  
  if (any(out <= 0)) {
    stop(
      paste(
        "The DRPT assumes r(z) > 0.",
        "At least one pooled observation has r(z) <= 0."
      )
    )
  }
  
  out
}

.lambda_hat_v2 = function(rZ, n, m, tol = 1e-12) {
  
  objective = function(lambda) {
    
    sum(
      1 / (
        n +
          m * lambda * rZ
      )
    ) - 1
  }
  
  lower = 0
  upper = 1
  
  while (objective(upper) > 0) {
    
    upper = 2 * upper
    
    if (!is.finite(upper) || upper > 1e16) {
      stop("Could not bracket the root for lambda_hat.")
    }
  }
  
  uniroot(
    objective,
    interval = c(lower, upper),
    tol = tol
  )$root
}

.kernel_matrix = function(Z, kernel) {
  
  Z = as.matrix(Z)
  
  N = nrow(Z)
  
  K = matrix(
    0,
    nrow = N,
    ncol = N
  )
  
  for (i in seq_len(N)) {
    
    K[i, i] = kernel(
      Z[i, ],
      Z[i, ]
    )
    
    if (i < N) {
      
      for (j in seq.int(i + 1L, N)) {
        
        kij = kernel(
          Z[i, ],
          Z[j, ]
        )
        
        K[i, j] = kij
        K[j, i] = kij
      }
    }
  }
  
  if (any(!is.finite(K))) {
    stop("Kernel matrix contains non-finite values.")
  }
  
  K
}


.shiftedMMD_v2_indices = function(
    idx_X,
    idx_Y,
    rZ,
    K,
    lambda_hat,
    n,
    m
) {
  
  tau_nm = as.numeric(n) / as.numeric(m)
  
  rX = rZ[idx_X]
  rY = rZ[idx_Y]
  
  wX =
    lambda_hat * rX /
    (
      tau_nm +
        lambda_hat * rX
    )
  
  wY =
    1 /
    (
      tau_nm +
        lambda_hat * rY
    )
  
  
  
  KXX = K[
    idx_X,
    idx_X,
    drop = FALSE
  ]
  
  KYY = K[
    idx_Y,
    idx_Y,
    drop = FALSE
  ]
  
  KXY = K[
    idx_X,
    idx_Y,
    drop = FALSE
  ]
  
  
  ## ------------------------------------
  ## First term
  ## ------------------------------------
  
  first_all =
    drop(
      crossprod(
        wX,
        KXX %*% wX
      )
    )
  
  first_diag =
    sum(
      diag(KXX) *
        wX^2
    )
  
  first_term =
    (first_all - first_diag) /
    n^2
  
  
  ## ------------------------------------
  ## Second term
  ## ------------------------------------
  
  second_all =
    drop(
      crossprod(
        wY,
        KYY %*% wY
      )
    )
  
  second_diag =
    sum(
      diag(KYY) *
        wY^2
    )
  
  second_term =
    (second_all - second_diag) /
    m^2
  
  
  ## ------------------------------------
  ## Mixed X-Y term
  ## ------------------------------------
  
  mixed_term =
    drop(
      crossprod(
        wX,
        KXY %*% wY
      )
    ) /
    (n * m)
  
  
  scale =
    (
      as.numeric(n + m) /
        as.numeric(m)
    )^2
  
  
  U =
    scale *
    (
      first_term +
        second_term -
        2 * mixed_term
    )
  
  
  return(U)
}

shiftedMMD_v2 = function(
    X,
    Y,
    r,
    kernel
) {
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  n = nrow(X)
  m = nrow(Y)
  
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.")
  }
  
  Z = rbind(
    X,
    Y
  )

  rZ = .eval_r_rows(
    Z,
    r
  )
  
  lambda_hat = .lambda_hat_v2(
    rZ = rZ,
    n = n,
    m = m
  )
  
  K = .kernel_matrix(
    Z,
    kernel
  )
  
  idx_X = seq_len(n)
  
  idx_Y = n + seq_len(m)
  
  .shiftedMMD_v2_indices(
    idx_X = idx_X,
    idx_Y = idx_Y,
    rZ = rZ,
    K = K,
    lambda_hat = lambda_hat,
    n = n,
    m = m
  )
}

.row_keys = function(Z) {
  
  Z = as.matrix(Z)
  
  apply(
    Z,
    1,
    function(z) {
      
      paste(
        sprintf("%.17g", z),
        collapse = "\034"
      )
    }
  )
}


## ============================================================
## Density Ratio Permutation Test
## ============================================================

DRPT_v2 = function(
    X,
    Y,
    r,
    kernel,
    H = 99,
    S = 50,
    details = FALSE
) {
  
  if (!requireNamespace(
    "DRPT",
    quietly = TRUE
  )) {
    stop("Package 'DRPT' is required.")
  }
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  n = nrow(X)
  m = nrow(Y)
  
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.")
  }
  
  if (n < 1 || m < 1) {
    stop("Both X and Y must contain observations.")
  }
  
  Z0 = rbind(
    X,
    Y
  )
  
  data = DRPT::starSampler(
    X = X,
    Y = Y,
    r = r,
    H = H,
    S = S
  )

  
  rZ = .eval_r_rows(
    Z0,
    r
  )
  
  
  lambda_hat = .lambda_hat_v2(
    rZ = rZ,
    n = n,
    m = m
  )
  

  
  K = .kernel_matrix(
    Z0,
    kernel
  )
  
  
  keys0 = .row_keys(Z0)
  
  compute_one = function(Zh) {
    
    Zh = as.matrix(Zh)
    
    keysh = .row_keys(Zh)
    
    idx = match(
      keysh,
      keys0
    )
    
    if (anyNA(idx)) {
      stop(
        paste(
          "Could not match a permuted observation",
          "to the original pooled dataset."
        )
      )
    }
    
    idx_X = idx[
      seq_len(n)
    ]
    
    idx_Y = idx[
      n + seq_len(m)
    ]
    
    .shiftedMMD_v2_indices(
      idx_X = idx_X,
      idx_Y = idx_Y,
      rZ = rZ,
      K = K,
      lambda_hat = lambda_hat,
      n = n,
      m = m
    )
  }
  

  
  T_hats = vapply(
    data,
    compute_one,
    numeric(1)
  )
  
  
  T_obs = T_hats[1]
  
  T_perm = T_hats[-1]
  
  
  p_hat =
    (
      1 +
        sum(T_perm >= T_obs)
    ) /
    (H + 1)
  
  
  if (!details) {
    return(p_hat)
  }
  
  
  return(
    list(
      p.value = p_hat,
      statistic = T_obs,
      permutation.statistics = T_perm,
      lambda_hat = lambda_hat,
      n = n,
      m = m,
      tau_nm = as.numeric(n) / as.numeric(m),
      scale_v2 =
        (
          as.numeric(n + m) /
            as.numeric(m)
        )^2
    )
  )
}