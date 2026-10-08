## ============================================================
## DRPT v3
##
## Pure-R star sampler with a user-controlled number of
## proposed X-Y pairs per Markov step.
##
## When number_pairs_switched = min(n, m)
## gives the same Markov transition mechanism as the
## original DRPT::starSampler().
## ============================================================

.check_positive_integer_v3 = function(x, name) {
  
  if (
    length(x) != 1L ||
    !is.numeric(x) ||
    !is.finite(x) ||
    x < 1 ||
    x != floor(x)
  ) {
    stop(
      paste0(
        "'", name,
        "' must be a positive integer."
      )
    )
  }
  
  as.integer(x)
}

.starSampler_v3_indices = function(
    rZ,
    n,
    m,
    H = 99,
    S = 50,
    number_pairs_switched = min(n, m)
) {
  
  H = .check_positive_integer_v3(
    H,
    "H"
  )
  
  S = .check_positive_integer_v3(
    S,
    "S"
  )
  
  number_pairs_switched =
    .check_positive_integer_v3(
      number_pairs_switched,
      "number_pairs_switched"
    )
  
  
  K = min(n, m)
  
  if (number_pairs_switched > K) {
    stop(
      paste0(
        "'number_pairs_switched' cannot exceed min(n, m) = ",
        K,
        "."
      )
    )
  }
  
  
  N = n + m
  
  if (length(rZ) != N) {
    stop("length(rZ) must equal n + m.")
  }
  
  if (
    any(!is.finite(rZ)) ||
    any(rZ <= 0)
  ) {
    stop("All density-ratio values must be finite and positive.")
  }
  
  log_rZ = log(rZ)
  
  q = number_pairs_switched
  
  
  ## ------------------------------------
  ## One step of Algorithm 1
  ## ------------------------------------
  
  one_step = function(state) {
    
    pos_X = sample.int(
      n,
      size = q,
      replace = FALSE
    )
    
    pos_Y =
      n +
      sample.int(
        m,
        size = q,
        replace = FALSE
      )
    
    obs_X = state[pos_X]
    obs_Y = state[pos_Y]
    
    p_swap = plogis(
      log_rZ[obs_X] -
        log_rZ[obs_Y]
    )
    
    
    do_swap =
      runif(q) <
      p_swap

    
    if (any(do_swap)) {
      
      swap_X = pos_X[do_swap]
      swap_Y = pos_Y[do_swap]
      
      tmp = state[swap_X]
      
      state[swap_X] =
        state[swap_Y]
      
      state[swap_Y] =
        tmp
    }
    
    
    return(state)
  }
  
  
  ## ------------------------------------
  ## Run S Markov steps
  ## ------------------------------------
  
  run_S_steps = function(state) {
    
    for (s in seq_len(S)) {
      
      state =
        one_step(state)
    }
    
    return(state)
  }
  
  
  idx_original =
    seq_len(N)

  
  idx_star =
    run_S_steps(
      idx_original
    )

  
  out =
    vector(
      "list",
      H + 1L
    )
  
  
  out[[1L]] =
    idx_original
  
  
  for (h in seq_len(H)) {
    
    out[[h + 1L]] =
      run_S_steps(
        idx_star
      )
  }
  
  
  return(out)
}


## ============================================================
## Pure-R starSampler
## ============================================================

starSampler_v3 = function(
    X,
    Y,
    r,
    H = 99,
    S = 50,
    number_pairs_switched = NULL
) {
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  n = nrow(X)
  m = nrow(Y)
  
  
  if (ncol(X) != ncol(Y)) {
    stop(
      "X and Y must have the same number of columns."
    )
  }
  
  if (n < 1 || m < 1) {
    stop(
      "Both X and Y must contain observations."
    )
  }
  
  
  if (is.null(number_pairs_switched)) {
    
    number_pairs_switched =
      min(n, m)
  }
  
  
  Z0 =
    rbind(
      X,
      Y
    )
  
  rZ =
    .eval_r_rows(
      Z0,
      r
    )
  
  
  idx_list =
    .starSampler_v3_indices(
      rZ = rZ,
      n = n,
      m = m,
      H = H,
      S = S,
      number_pairs_switched =
        number_pairs_switched
    )
  
  
  data =
    lapply(
      idx_list,
      function(idx) {
        
        Z0[
          idx,
          ,
          drop = FALSE
        ]
      }
    )
  
  
  return(data)
}


## ============================================================
## Density Ratio Permutation Test v3
## ============================================================

DRPT_v3 = function(
    X,
    Y,
    r,
    kernel,
    H = 99,
    S = 50,
    number_pairs_switched = NULL,
    details = FALSE
) {
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  n = nrow(X)
  m = nrow(Y)
  
  
  if (ncol(X) != ncol(Y)) {
    stop(
      "X and Y must have the same number of columns."
    )
  }
  
  if (n < 1 || m < 1) {
    stop(
      "Both X and Y must contain observations."
    )
  }
  
  
  if (is.null(number_pairs_switched)) {
    
    number_pairs_switched =
      min(n, m)
  }
  
  
  number_pairs_switched =
    .check_positive_integer_v3(
      number_pairs_switched,
      "number_pairs_switched"
    )
  
  
  if (
    number_pairs_switched >
    min(n, m)
  ) {
    
    stop(
      paste0(
        "'number_pairs_switched' cannot exceed min(n, m) = ",
        min(n, m),
        "."
      )
    )
  }
  
  
  Z0 =
    rbind(
      X,
      Y
    )
  
  
  rZ =
    .eval_r_rows(
      Z0,
      r
    )

  
  lambda_hat =
    .lambda_hat_v2(
      rZ = rZ,
      n = n,
      m = m
    )
  
  
  K =
    .kernel_matrix(
      Z0,
      kernel
    )

  
  idx_data =
    .starSampler_v3_indices(
      rZ = rZ,
      n = n,
      m = m,
      H = H,
      S = S,
      number_pairs_switched =
        number_pairs_switched
    )
  
  compute_one = function(idx) {
    
    idx_X =
      idx[
        seq_len(n)
      ]
    
    idx_Y =
      idx[
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
  
  T_hats =
    vapply(
      idx_data,
      compute_one,
      numeric(1)
    )
  
  
  T_obs =
    T_hats[1L]
  
  T_perm =
    T_hats[-1L]
  
  
  p_hat =
    (
      1 +
        sum(
          T_perm >= T_obs
        )
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
      H = H,
      S = S,
      number_pairs_switched =
        number_pairs_switched,
      tau_nm =
        as.numeric(n) /
        as.numeric(m),
      scale_v2 =
        (
          as.numeric(n + m) /
            as.numeric(m)
        )^2
    )
  )
}