# ==============================================================================
# Simulation study: ZIDM-LVM, ZIMN-LVM, and ZINB-LVM generators
# ==============================================================================
# This file consolidates three original simulation programs.
# Shared model-fitting and evaluation functions are defined once below.
# Generator-specific implementations are scoped to each simulation runner.
# This preserves generator-specific behavior without overwriting functions.
#
# Usage on a SLURM array or interactively:
#   Rscript simulation_main.R ZIDM
#   Rscript simulation_main.R ZIMN
#   Rscript simulation_main.R ZINB
# SLURM_ARRAY_TASK_ID selects the scenario as in the original programs.
# Sourcing this script defines functions but does not start a simulation.
#
# Important: the original numerical algorithms and scenario configurations
# have been retained. Check the three scenario settings before publication.
# ==============================================================================

# ==============================================================================
# FUSED CI VERSION
# - Keeps the newly added Koslovsky ZIDM paper method + posterior credible intervals
# - Adds CI evaluation for ZINB_gllvm, ZIMN_soft/hard, ZIDM_soft/hard
# - gllvm ZINB uses coefficient-level SE validity: invalid beta coefficients are NA
#   and excluded individually from coverage; other valid beta coefficients remain.
# ==============================================================================
# ==============================================================================
# ✅ 1. HPC user library setup and package dependencies
# ==============================================================================
dir.create("~/R/library_452", recursive = TRUE, showWarnings = FALSE)
.libPaths("~/R/library_452")

if(!"ZIPPCAlnm" %in% rownames(installed.packages(lib.loc="~/R/library_452"))){
  if(!"remotes" %in% rownames(installed.packages(lib.loc="~/R/library_452"))){
    install.packages("remotes", lib="~/R/library_452")
  }
  remotes::install_github("YanyZeng/ZIPPCAlnm", lib="~/R/library_452")
}

############################################################
## ZIDM
############################################################
.libPaths(c(
  path.expand("~/R/library_452"),
  .libPaths()
))

if(!"ZIDM" %in% rownames(installed.packages(lib.loc="~/R/library_452"))){
  if(!"remotes" %in% rownames(installed.packages(lib.loc="~/R/library_452"))){
    install.packages("remotes", lib="~/R/library_452")
  }
  remotes::install_github("mkoslovsky/ZIDM", lib="~/R/library_452")
}

library(ZIDM)
library(MASS)
library(compositions)
library(robCompositions)
library(TMB)
library(gllvm)
library(ZIPPCAlnm)
library(vegan)   # Procrustes residuals for latent ordination recovery



# ==================== Shared simulation utilities ====================

generate_ZINB_data <- function(
    n,
    m,
    k,
    p = 2,
    alpha = NULL,
    beta = NULL,
    Lambda = NULL,
    zeta = NULL,
    zeta_range = c(0.05, 0.15),
    X = NULL,
    U_sampler = NULL,
    alpha_i = NULL,
    phi = 1,
    zero_inflation = TRUE,
    seed = NULL
){
  
  ############################################################
  ## 0. Seed
  ############################################################
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  
  ############################################################
  ## 1. Covariates X
  ############################################################
  
  if(is.null(X)){
    
    if(p == 2){
      
      X <- cbind(
        x1 = rnorm(n, 0, 1),
        x2 = rexp(n, rate = 1)
      )
      
    } else {
      
      X <- matrix(
        rnorm(n * p),
        nrow = n,
        ncol = p
      )
      
      colnames(X) <- paste0("x", seq_len(p))
    }
    
  } else {
    
    X <- as.matrix(X)
    
    if(nrow(X) != n){
      stop("X must have n rows.")
    }
    
    p <- ncol(X)
  }
  
  
  ############################################################
  ## 2. Species-specific intercept beta_0j
  ############################################################
  
  if(is.null(alpha)){
    
    alpha <- runif(
      m,
      min = -1,
      max = 1
    )
    
  } else {
    
    alpha <- as.numeric(alpha)
    
    if(length(alpha) != m){
      stop("alpha must have length m.")
    }
  }
  
  beta0 <- alpha
  
  
  ############################################################
  ## 3. Environmental coefficients beta_j
  ##
  ## dimension: m x p
  ############################################################
  
  if(is.null(beta)){
    
    beta <- matrix(
      runif(
        m * p,
        min = -1,
        max = 1
      ),
      nrow = m,
      ncol = p
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if(!all(dim(beta) == c(m, p))){
      stop("beta must have dimension m x p.")
    }
  }
  
  
  ############################################################
  ## 4. Latent loadings gamma_j
  ##
  ## dimension: m x k
  ############################################################
  
  if(is.null(Lambda)){
    
    Lambda <- matrix(
      runif(
        m * k,
        min = -2,
        max = 2
      ),
      nrow = m,
      ncol = k
    )
    
  } else {
    
    Lambda <- as.matrix(Lambda)
    
    if(!all(dim(Lambda) == c(m, k))){
      stop("Lambda must have dimension m x k.")
    }
  }
  
  
  ############################################################
  ## 5. Latent variables U
  ##
  ## u_i ~ N(0, I_k)
  ##
  ## dimension: n x k
  ############################################################
  
  if(is.null(U_sampler)){
    
    U <- matrix(
      rnorm(n * k),
      nrow = n,
      ncol = k
    )
    
  } else {
    
    U <- U_sampler(n, k)
    
    U <- as.matrix(U)
    
    if(!all(dim(U) == c(n, k))){
      stop("U_sampler must return an n x k matrix.")
    }
  }
  
  
  ############################################################
  ## 6. Sample-specific row effect alpha_i
  ##
  ## Used by ZINB to control sequencing depth.
  ##
  ## If not supplied, use zero row effects.
  ############################################################
  
  if(is.null(alpha_i)){
    
    alpha_i <- rep(0, n)
    
  } else {
    
    alpha_i <- as.numeric(alpha_i)
    
    if(length(alpha_i) == 1){
      
      alpha_i <- rep(
        alpha_i,
        n
      )
      
    } else if(length(alpha_i) != n){
      
      stop("alpha_i must be either scalar or length n.")
    }
  }
  
  
  ############################################################
  ## 7. Fixed-effect contribution X beta
  ##
  ## X:    n x p
  ## beta: m x p
  ##
  ## result: n x m
  ############################################################
  
  XB <- X %*% t(beta)
  
  
  ############################################################
  ## 8. Latent contribution U Gamma
  ##
  ## U:      n x k
  ## Lambda: m x k
  ##
  ## result: n x m
  ############################################################
  
  UG <- U %*% t(Lambda)
  
  
  ############################################################
  ## 9. Linear predictor
  ##
  ## eta_ij =
  ## alpha_i
  ## + beta_0j
  ## + x_i^T beta_j
  ## + u_i^T gamma_j
  ############################################################
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  eta <- eta +
    matrix(
      alpha_i,
      nrow = n,
      ncol = m,
      byrow = FALSE
    )
  
  eta <- eta +
    XB +
    UG
  
  
  ############################################################
  ## 10. Mean parameter
  ##
  ## mu_ij = exp(eta_ij)
  ############################################################
  
  mu <- exp(eta)
  
  
  ############################################################
  ## 11. TRUE composition
  ##
  ## alpha_i cancels after normalization:
  ##
  ## rho_ij =
  ## exp(beta0_j + x_i beta_j + u_i gamma_j)
  ## -------------------------------------------------
  ## sum_l exp(beta0_l + x_i beta_l + u_i gamma_l)
  ############################################################
  
  eta_comp <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  ) +
    XB +
    UG
  
  
  ############################################################
  ## stable softmax
  ############################################################
  
  rho_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    row_max <- max(
      eta_comp[i, ]
    )
    
    tmp <- exp(
      eta_comp[i, ] - row_max
    )
    
    rho_true[i, ] <-
      tmp / sum(tmp)
  }
  
  
  ############################################################
  ## 12. Taxon-specific zero-inflation probabilities
  ##
  ## IMPORTANT:
  ##
  ## zeta_j = structural-zero probability for taxon j
  ##
  ## zeta = (zeta_1, ..., zeta_m)
  ##
  ## dimension: m
  ############################################################
  
  if(is.null(zeta)){
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if(length(zeta) != m){
      
      stop(
        paste0(
          "For taxon-specific ZINB zero inflation, ",
          "zeta must have length m = ",
          m,
          "."
        )
      )
    }
  }
  
  
  ############################################################
  ## check zeta
  ############################################################
  
  if(any(!is.finite(zeta))){
    stop("All zeta values must be finite.")
  }
  
  if(any(zeta < 0 | zeta > 1)){
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  ############################################################
  ## 13. Structural-zero indicators
  ##
  ## z_ij ~ Bernoulli(zeta_j)
  ##
  ## Here:
  ##
  ## structural_zero_ij = 1
  ## means Y_ij is forced to zero.
  ############################################################
  
  if(zero_inflation){
    
    ##########################################################
    ## Create n x m matrix:
    ##
    ## row 1: zeta_1 ... zeta_m
    ## row 2: zeta_1 ... zeta_m
    ## ...
    ##########################################################
    
    zeta_mat <- matrix(
      rep(
        zeta,
        each = n
      ),
      nrow = n,
      ncol = m
    )
    
    
    structural_zero <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(zeta_mat)
      ),
      
      nrow = n,
      ncol = m
    )
    
  } else {
    
    structural_zero <- matrix(
      0L,
      nrow = n,
      ncol = m
    )
  }
  
  
  ############################################################
  ## 14. Generate negative-binomial counts
  ##
  ## NB parameterization:
  ##
  ## E(Y_ij) = mu_ij
  ##
  ## Var(Y_ij) =
  ## mu_ij + mu_ij^2 / phi
  ############################################################
  
  if(length(phi) == 1){
    
    phi_mat <- matrix(
      phi,
      nrow = n,
      ncol = m
    )
    
  } else if(length(phi) == m){
    
    phi_mat <- matrix(
      rep(
        phi,
        each = n
      ),
      nrow = n,
      ncol = m
    )
    
  } else {
    
    stop("phi must be scalar or length m.")
  }
  
  
  Y_sampling <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  
  for(i in seq_len(n)){
    
    for(j in seq_len(m)){
      
      Y_sampling[i, j] <- rnbinom(
        1,
        mu = mu[i, j],
        size = phi_mat[i, j]
      )
    }
  }
  
  
  ############################################################
  ## 15. Apply structural zeros
  ############################################################
  
  Y <- Y_sampling
  
  if(zero_inflation){
    
    Y[
      structural_zero == 1
    ] <- 0L
  }
  
  
  ############################################################
  ## 16. Zero diagnostics
  ############################################################
  
  observed_zero <-
    Y == 0
  
  sampling_zero <-
    (Y_sampling == 0) &
    (structural_zero == 0)
  
  
  observed_zero_prop <-
    mean(observed_zero)
  
  structural_zero_prop <-
    mean(structural_zero == 1)
  
  sampling_zero_prop <-
    mean(sampling_zero)
  
  
  ############################################################
  ## 17. Return
  ############################################################
  
  return(list(
    
    ##########################################################
    ## observed data
    ##########################################################
    
    Y = Y,
    
    X = X,
    
    
    ##########################################################
    ## latent quantities
    ##########################################################
    
    U = U,
    
    eta = eta,
    
    mu = mu,
    
    rho_true = rho_true,
    
    
    ##########################################################
    ## zero-inflation quantities
    ##########################################################
    
    zeta = zeta,
    
    zeta_mat = if(zero_inflation){
      zeta_mat
    } else {
      matrix(
        0,
        nrow = n,
        ncol = m
      )
    },
    
    structural_zero =
      structural_zero,
    
    
    ##########################################################
    ## counts before structural zero inflation
    ##########################################################
    
    Y_sampling =
      Y_sampling,
    
    
    ##########################################################
    ## zero proportions
    ##########################################################
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    
    ##########################################################
    ## true parameters
    ##########################################################
    
    truth = list(
      
      alpha_i =
        alpha_i,
      
      beta0 =
        beta0,
      
      beta =
        beta,
      
      gamma =
        Lambda,
      
      U =
        U,
      
      zeta =
        zeta,
      
      phi =
        phi,
      
      rho =
        rho_true
    )
  ))
}

fit_ZIDM_TMB <- function(
    Y,
    X = NULL,
    k = 2,
    family = 2,
    add_intercept = FALSE,
    xr = NULL,
    offset = NULL,
    sizes = NULL,
    extra = c(0, 0),
    
    init_r0 = NULL,
    init_b = NULL,
    init_theta = NULL,
    init_u = NULL,
    init_Au = NULL,
    init_lg_Si = NULL,
    init_lgt_zeta = NULL,
    
    fix_r0 = TRUE,
    optimizer = c("optim", "nlminb"),
    optim_method = "L-BFGS-B",
    optim_control = list(maxit = 30000, factr = 1e5, reltol = 1e-10),
    nlminb_control = list(iter.max = 30000, eval.max = 30000, rel.tol = 1e-10),
    
    DLL = "zimnVA",
    silent = TRUE
){
  optimizer <- match.arg(optimizer)
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  
  # =========================================================
  # ⭐ FIXED X handling (ALWAYS include intercept)
  # =========================================================
  if (is.null(X)) {
    
    # ⭐ 关键修复：即使没有 covariate，也加 intercept
    X <- matrix(1, n, 1)
    
  } else {
    
    X <- as.matrix(X)
    if (nrow(X) != n) stop("X must have the same number of rows as Y.")
    
    if (add_intercept) {
      X <- cbind(1, X)
    }
  }
  
  q <- ncol(X)
  
  # -----------------------------
  # xr handling
  # -----------------------------
  if (is.null(xr)) {
    xr <- matrix(0, 1, p)
  } else {
    xr <- as.matrix(xr)
  }
  
  # -----------------------------
  # offset
  # -----------------------------
  if (is.null(offset)) {
    offset <- matrix(0, n, p)
  } else {
    offset <- as.matrix(offset)
    if (!all(dim(offset) == c(n, p))) {
      stop("offset must have dimension n x p.")
    }
  }
  
  # -----------------------------
  # sizes
  # -----------------------------
  if (is.null(sizes)) {
    sizes <- rowSums(Y)
  }
  sizes <- as.integer(sizes)
  
  if (ncol(xr) != p) {
    stop("xr must have p columns so that r0 %*% xr has dimension n x p.")
  }
  
  # -----------------------------
  # theta length
  # -----------------------------
  theta_len <- if (k > 0) {
    p * k - (k * (k - 1)) / 2
  } else 0
  
  # -----------------------------
  # Au length
  # -----------------------------
  Au_len <- if (k > 0) {
    k * n + n * (k * (k - 1) / 2)
  } else 0
  
  # -----------------------------
  # initial values
  # -----------------------------
  if (is.null(init_r0)) {
    init_r0 <- matrix(0, n, nrow(xr))
  } else {
    init_r0 <- as.matrix(init_r0)
  }
  
  if (!all(dim(init_r0) == c(n, nrow(xr)))) {
    stop(sprintf("init_r0 must have dimension %d x %d.", n, nrow(xr)))
  }
  
  if (is.null(init_b)) {
    init_b <- matrix(0, q, p)
  } else {
    init_b <- as.matrix(init_b)
  }
  
  if (!all(dim(init_b) == c(q, p))) {
    stop(sprintf("init_b must have dimension %d x %d.", q, p))
  }
  
  if (is.null(init_theta)) {
    init_theta <- if (theta_len > 0) rnorm(theta_len, 0, 0.1) else numeric(0)
  } else {
    init_theta <- as.numeric(init_theta)
  }
  
  if (length(init_theta) != theta_len) {
    stop(sprintf("init_theta must have length %d.", theta_len))
  }
  
  if (is.null(init_u)) {
    init_u <- if (k > 0) matrix(rnorm(n * k, 0, 0.1), n, k) else matrix(0, n, 0)
  } else {
    init_u <- as.matrix(init_u)
  }
  
  if (!all(dim(init_u) == c(n, k))) {
    stop(sprintf("init_u must have dimension %d x %d.", n, k))
  }
  
  if (is.null(init_Au)) {
    if (k > 0) {
      init_Au <- c(
        rep(log(0.2), k * n),
        rep(0, n * (k * (k - 1) / 2))
      )
    } else {
      init_Au <- numeric(0)
    }
  } else {
    init_Au <- as.numeric(init_Au)
  }
  
  if (length(init_Au) != Au_len) {
    stop(sprintf("init_Au must have length %d.", Au_len))
  }
  
  if (is.null(init_lg_Si)) {
    init_lg_Si <- log(pmax(sizes, 1))
  } else {
    init_lg_Si <- as.numeric(init_lg_Si)
  }
  
  if (length(init_lg_Si) != n) {
    stop(sprintf("init_lg_Si must have length %d.", n))
  }
  
  if (is.null(init_lgt_zeta)) {
    init_lgt_zeta <- rep(0, p)
  } else {
    init_lgt_zeta <- as.numeric(init_lgt_zeta)
  }
  
  if (length(init_lgt_zeta) != p) {
    stop(sprintf("init_lgt_zeta must have length %d.", p))
  }
  
  # -----------------------------
  # data
  # -----------------------------
  data <- list(
    y = Y,
    x = X,
    xr = xr,
    offset = offset,
    num_lv = k,
    family = family,
    extra = extra,
    sizes = sizes
  )
  
  parameters <- list(
    r0 = init_r0,
    b = init_b,
    theta = init_theta,
    u = init_u,
    Au = init_Au,
    lg_Si = init_lg_Si,
    lgt_zeta = init_lgt_zeta
  )
  
  # -----------------------------
  # map
  # -----------------------------
  map <- NULL
  if (fix_r0) {
    map <- list(
      r0 = factor(rep(NA, length(parameters$r0)))
    )
  }
  
  # -----------------------------
  # TMB object
  # -----------------------------
  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    map = map,
    DLL = DLL,
    silent = silent
  )
  
  # -----------------------------
  # optimize
  # -----------------------------
  if (optimizer == "optim") {
    opt <- optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      method = optim_method,
      control = optim_control
    )
    convergence <- opt$convergence
    objective_value <- opt$value
    par_hat <- opt$par
  } else {
    opt <- nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = nlminb_control
    )
    convergence <- opt$convergence
    objective_value <- opt$objective
    par_hat <- opt$par
  }
  
  obj$env$last.par.best <- par_hat
  rep <- obj$report()
  par_list <- obj$env$parList(par_hat)
  out <- list(
    
    obj = obj,
    opt = opt,
    
    # all raw estimated parameters
    par_hat = par_list,
    
    # regression parameters
    b_hat = par_list$b,
    
    beta0_hat = par_list$b[1, ],
    
    beta_hat =
      if (nrow(par_list$b) > 1) {
        par_list$b[-1, , drop = FALSE]
      } else {
        matrix(numeric(0), 0, p)
      },
    
    # latent structure
    u_hat = par_list$u,
    theta_hat = par_list$theta,
    Au_hat = par_list$Au,
    
    # ZIDM augmentation parameters
    lg_Si_hat = par_list$lg_Si,
    Si_hat = rep$Si,
    
    # zero-inflation parameters
    lgt_zeta_hat = par_list$lgt_zeta,
    zeta_hat = plogis(par_list$lgt_zeta),
    
    # actual loading matrix
    gamma_hat = rep$newlam,
    
    # fitted quantities
    rho0_hat = rep$rho0,
    rhoZI_hat = rep$rhoZI,
    rho_hat = rep$rho0,
    
    eta_hat = rep$eta,
    cQ_hat = rep$cQ,
    
    logl = rep$logl,
    nll = rep$nll,
    
    convergence = convergence,
    objective = objective_value
  )
  
  if (k > 0) {
    out$ordination <- rep$u
    out$loadings <- rep$newlam
  }
  
  return(out)
}

fit_TMB <- function(
    Y,
    X = NULL,
    k = 2,
    model = c("ZIDM","ZIMN"),
    extra = c(0,0),
    
    init_method = c(
      "FA",
      "PCA",
      "RANDOM",
      "DS"
    ),
    
    init_u = NULL,
    init_lgt_zeta = NULL,
    init_lg_Si = NULL,
    init_b = NULL
){
  
  model <- match.arg(model)
  
  init_method <- match.arg(init_method)
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  
  m <- ncol(Y)
  
  ##########################################################
  # X
  ##########################################################
  
  if (!is.null(X)) {
    
    X <- as.matrix(X)
    
    if (nrow(X) != n) {
      stop("X must have the same number of rows as Y.")
    }
  }
  
  # number of rows needed in b:
  # 1 intercept + number of covariates
  q_fit <- if (is.null(X)) {
    1
  } else {
    ncol(X) + 1
  }
  
  ##########################################################
  # family
  ##########################################################
  
  family <- switch(
    model,
    "ZIMN" = 1,
    "ZIDM" = 2
  )
  
  ##########################################################
  # fit0
  ##########################################################
  
  fit0_init <- initialize_fit0(
    
    Y = Y,
    
    X = X,
    
    family = family,
    
    extra = c(0,0)
  )
  
  ##########################################################
  # shared nuisance
  ##########################################################
  
  if(!is.null(fit0_init)){
    
    if(is.null(init_b)){
      
      init_b <- fit0_init$b
    }
    
    if(is.null(init_lgt_zeta)){
      
      init_lgt_zeta <- fit0_init$lgt_zeta
    }
    
    if(is.null(init_lg_Si)){
      
      init_lg_Si <- fit0_init$lg_Si
    }
  }
  
  ##########################################################
  # fallback
  ##########################################################
  
  if(is.null(init_b)){
    
    init_b <- matrix(
      0,
      q_fit,
      m
    )
  }
  
  if(is.null(init_lgt_zeta)){
    
    init_lgt_zeta <- rnorm(m, -3, 0.5)
  }
  
  if(is.null(init_lg_Si)){
    
    init_lg_Si <- log(rowSums(Y) + 1)
  }
  
  ##########################################################
  # latent variables
  ##########################################################
  
  if(is.null(init_u)){
    
    ########################################################
    # FA initialization
    ########################################################
    
    if(init_method == "FA"){
      
      if(is.null(fit0_init)){
        
        clr_y <- cenLR(Y + 1)$x.clr
        
        pca <- prcomp(
          clr_y,
          center = TRUE,
          scale. = FALSE
        )
        
        init_u <- scale(
          pca$x[,1:k,drop=FALSE]
        )
        
      } else {
        
        rho_fit <- fit0_init$rho_fit
        
        rho_obs <- sweep(
          Y + 0.5,
          1,
          rowSums(Y + 0.5),
          "/"
        )
        
        clr_obs <- cenLR(rho_obs)$x.clr
        
        clr_fit <- cenLR(rho_fit)$x.clr
        
        res_clr <- clr_obs - clr_fit
        
        res_clr[!is.finite(res_clr)] <- 0
        
        pca_fa <- prcomp(
          res_clr,
          center = TRUE,
          scale. = FALSE
        )
        
        init_u <- scale(
          pca_fa$x[,1:k,drop=FALSE]
        )
      }
    }
    
    ########################################################
    # PCA initialization
    ########################################################
    
    if(init_method == "PCA"){
      
      clr_y <- cenLR(Y + 1)$x.clr
      
      pca <- prcomp(
        clr_y,
        center = TRUE,
        scale. = FALSE
      )
      
      init_u <- scale(
        pca$x[,1:k,drop=FALSE]
      )
    }
    
    ########################################################
    # RANDOM initialization
    ########################################################
    
    if(init_method == "RANDOM"){
      
      init_u <- matrix(
        rnorm(n*k,0,0.1),
        n,
        k
      )
      
      init_u <- scale(init_u)
    }
    
    ########################################################
    # DS initialization
    ########################################################
    
    if(init_method == "DS"){
      
      if (is.null(fit0_init)) {
        stop("DS initialization requires a successful fit0.")
      }
      
      DS_init <- DS_FA_initialize(
        
        Y = Y,
        
        fit0_init = fit0_init,
        
        k = k,
        
        model = model,
        
        route = "B",
        
        jitter = FALSE,
        
        seed = 123
      )
      
      init_u <- DS_init$init_u
    }
  }
  
  ##########################################################
  # final fit
  ##########################################################
  fit <- fit_ZIDM_TMB(
    
    Y = Y,
    
    X = X,
    
    
    add_intercept = !is.null(X),
    
    k = k,
    
    family = family,
    
    extra = extra,
    
    init_b = init_b,
    
    init_u = init_u,
    
    init_lgt_zeta = init_lgt_zeta,
    
    init_lg_Si = init_lg_Si,
    
    fix_r0 = TRUE,
    
    optimizer = "optim",
    
    optim_method = "L-BFGS-B",
    
    optim_control = list(
      maxit = 50000,
      factr = 1e4,
      pgtol = 1e-8,
      trace = 0
    )
  )
  
  ##########################################################
  # recover all parameters
  ##########################################################
  
  par_hat <- fit$obj$env$parList(
    fit$opt$par
  )
  
  b_hat <- par_hat$b
  
  # first row = beta_0j
  beta0_hat <- b_hat[1, ]
  
  # remaining rows = environmental beta
  beta_hat <- if(nrow(b_hat) > 1){
    
    b_hat[-1, , drop = FALSE]
    
  } else {
    
    matrix(
      numeric(0),
      nrow = 0,
      ncol = m
    )
  }
  
  ##########################################################
  # outputs
  ##########################################################
  
  return(list(
    
    # regression
    beta0_hat = beta0_hat,
    
    beta_hat = beta_hat,
    
    b_hat = b_hat,
    
    # latent variables
    u_hat = fit$u_hat,
    
    gamma_hat = fit$gamma_hat,
    
    theta_hat = fit$theta_hat,
    
    Au_hat = fit$Au_hat,
    
    # zero inflation
    lgt_zeta_hat = fit$lgt_zeta_hat,
    
    zeta_hat = fit$zeta_hat,
    
    # ZIDM auxiliary parameters
    lg_Si_hat = fit$lg_Si_hat,
    
    Si_hat = fit$Si_hat,
    
    # compositions
    rho_hat = fit$rho0_hat,
    
    rhoZI_hat = fit$rhoZI_hat,
    
    # fitted predictor
    eta_hat = fit$eta_hat,
    
    cQ_hat = fit$cQ_hat,
    
    # old names, if existing simulation code uses them
    ordination = fit$u_hat,
    
    loadings = fit$gamma_hat,
    
    lgt_zeta = fit$lgt_zeta_hat,
    
    # optimization information
    convergence = fit$convergence,
    
    objective = fit$objective,
    
    logl = fit$logl,
    
    nll = fit$nll,
    
    # model information
    model = model,
    
    extra = extra,
    
    init_method = init_method,
    
    # everything
    all_parameters = par_hat,
    
    fit0_init = fit0_init,
    
    fit = fit
  ))
}

initialize_fit0 <- function(
    Y,
    X = NULL,
    family = 2,
    extra = c(0,0),
    DLL = "zimnVA",
    optimizer = "optim",
    optim_method = "L-BFGS-B"
){
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  ##########################################################
  # X dimension
  ##########################################################
  
  if (is.null(X)) {
    
    # only species-specific intercept
    q_fit <- 1
    
  } else {
    
    X <- as.matrix(X)
    
    if (nrow(X) != n) {
      stop("X must have the same number of rows as Y.")
    }
    
    # intercept + environmental covariates
    q_fit <- ncol(X) + 1
  }
  
  ##########################################################
  # initial b
  ##########################################################
  
  init_b0 <- matrix(
    rnorm(q_fit * p, 0, 0.1),
    q_fit,
    p
  )
  
  ##########################################################
  # ZIPPCAlnm-style ZI init
  ##########################################################
  
  pzero_col <- colMeans(Y == 0)
  
  pi_init <- t(
    ifelse(
      t(Y) == 0,
      pzero_col,
      0
    )
  )
  
  eta_init <- colMeans(pi_init)
  
  eta_init <- pmin(
    pmax(eta_init, 1e-4),
    1 - 1e-4
  )
  
  init_lgt_zeta0 <- qlogis(eta_init)
  
  ##########################################################
  # stable library size
  ##########################################################
  
  init_lg_Si0 <- log(rowSums(Y) + 1)
  
  ##########################################################
  # fit0
  ##########################################################
  
  fit0 <- try(
    
    fit_ZIDM_TMB(
      Y = Y,
      
      X = X,
      
      # intercept
      add_intercept = !is.null(X),
      
      k = 0,
      
      family = family,
      
      DLL = DLL,
      
      extra = extra,
      
      init_b = init_b0,
      
      init_lgt_zeta = init_lgt_zeta0,
      
      init_lg_Si = init_lg_Si0,
      
      fix_r0 = TRUE,
      
      optimizer = optimizer,
      
      optim_method = optim_method,
      
      optim_control = list(
        maxit = 30000,
        factr = 1e7,
        pgtol = 1e-4
      )
    ),
    
    silent = TRUE
  )
  
  ##########################################################
  # failed fit
  ##########################################################
  
  if(inherits(fit0, "try-error")){
    return(NULL)
  }
  
  ##########################################################
  # convergence check
  ##########################################################
  
  if(
    fit0$convergence != 0 ||
    any(is.na(fit0$rho_hat)) ||
    any(!is.finite(fit0$rho_hat))
  ){
    return(NULL)
  }
  
  ##########################################################
  # extract nuisance
  ##########################################################
  
  par0 <- fit0$obj$env$parList(
    fit0$opt$par
  )
  
  shared_b <- par0$b
  
  shared_zeta <- par0$lgt_zeta
  
  shared_lgSi <- par0$lg_Si
  
  ##########################################################
  # fitted composition
  ##########################################################
  
  rho_fit <- pmax(
    fit0$rho_hat,
    1e-8
  )
  
  rho_fit <- rho_fit / rowSums(rho_fit)
  
  ##########################################################
  # return
  ##########################################################
  
  return(list(
    
    fit0 = fit0,
    
    b = shared_b,
    
    lgt_zeta = shared_zeta,
    
    lg_Si = shared_lgSi,
    
    rho_fit = rho_fit
  ))
}

DS_FA_initialize <- function(
    Y,
    fit0_init,
    k,
    model = c("ZIDM", "ZIMN"),
    route = c("A", "B"),
    jitter = FALSE,
    jitter_sd = 0.1,
    seed = 123
){
  
  ##########################################################
  ## route / model selection
  ##########################################################
  
  route <- match.arg(route)
  model <- match.arg(model)
  
  set.seed(seed)
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  m <- ncol(Y)
  
  N <- rowSums(Y)
  
  
  ##########################################################
  ## ======================================================
  ## Step 1:
  ## recover fitted quantities
  ## ======================================================
  ##########################################################
  
  rho_fit <- fit0_init$rho_fit
  
  if(model == "ZIDM"){
    
    if(route == "A"){
      
      ######################################################
      ## gamma_ij = S_i * rho_ij
      ######################################################
      
      lg_Si <- fit0_init$lg_Si
      
      S_fit <- exp(lg_Si)
      
      gamma_fit <- sweep(
        rho_fit,
        1,
        S_fit,
        "*"
      )
      
    } else {
      
      ######################################################
      ## gamma_ij = exp(eta_ij)
      ######################################################
      
      eta_fit <- fit0_init$fit0$eta
      
      gamma_fit <- exp(eta_fit)
    }
    
    gamma_plus <- rowSums(gamma_fit)
  }
  
  
  ##########################################################
  ## ======================================================
  ## Step 2:
  ## Dunn–Smyth residuals
  ## ======================================================
  ##########################################################
  
  R_DS <- matrix(0, n, m)
  
  for(i in 1:n){
    
    for(j in 1:m){
      
      ######################################################
      ## observed count
      ######################################################
      
      y_ij <- Y[i, j]
      
      ######################################################
      ## total count
      ######################################################
      
      N_i <- N[i]
      
      
      ######################################################
      ## ==================================================
      ## ZIDM:
      ## marginal beta-binomial DS residual
      ## ==================================================
      ######################################################
      
      if(model == "ZIDM"){
        
        a_ij <- gamma_fit[i, j]
        
        b_ij <- gamma_plus[i] -
          gamma_fit[i, j]
        
        ####################################################
        ## numerical stability
        ####################################################
        
        a_ij <- max(a_ij, 1e-10)
        
        b_ij <- max(b_ij, 1e-10)
        
        ####################################################
        ## lower CDF: F(y-1)
        ####################################################
        
        if(y_ij == 0){
          
          F_low <- 0
          
        } else {
          
          x_low <- 0:(y_ij - 1)
          
          logpmf_low <- lchoose(
            N_i,
            x_low
          ) +
            lbeta(
              x_low + a_ij,
              N_i - x_low + b_ij
            ) -
            lbeta(a_ij, b_ij)
          
          F_low <- sum(exp(logpmf_low))
        }
        
        ####################################################
        ## upper CDF: F(y)
        ####################################################
        
        x_high <- 0:y_ij
        
        logpmf_high <- lchoose(
          N_i,
          x_high
        ) +
          lbeta(
            x_high + a_ij,
            N_i - x_high + b_ij
          ) -
          lbeta(a_ij, b_ij)
        
        F_high <- sum(exp(logpmf_high))
        
        
        ######################################################
        ## ==================================================
        ## ZIMN:
        ## marginal binomial DS residual
        ## ==================================================
        ######################################################
        
      } else if(model == "ZIMN"){
        
        ####################################################
        ## fitted multinomial probability
        ####################################################
        
        p_ij <- rho_fit[i, j]
        
        p_ij <- pmin(
          pmax(p_ij, 1e-10),
          1 - 1e-10
        )
        
        ####################################################
        ## lower CDF: F(y-1)
        ####################################################
        
        if(y_ij == 0){
          
          F_low <- 0
          
        } else {
          
          F_low <- pbinom(
            y_ij - 1,
            size = N_i,
            prob = p_ij
          )
        }
        
        ####################################################
        ## upper CDF: F(y)
        ####################################################
        
        F_high <- pbinom(
          y_ij,
          size = N_i,
          prob = p_ij
        )
      }
      
      
      ######################################################
      ## clipping
      ######################################################
      
      F_low <- pmin(
        pmax(F_low, 1e-10),
        1 - 1e-10
      )
      
      F_high <- pmin(
        pmax(F_high, 1e-10),
        1 - 1e-10
      )
      
      ######################################################
      ## Dunn–Smyth randomization
      ######################################################
      
      u_ij <- runif(
        1,
        F_low,
        F_high
      )
      
      ######################################################
      ## DS residual
      ######################################################
      
      R_DS[i, j] <- qnorm(u_ij)
    }
  }
  
  R_DS[!is.finite(R_DS)] <- 0
  
  
  ##########################################################
  ## ======================================================
  ## Step 3:
  ## Factor analysis
  ## ======================================================
  ##########################################################
  
  if(n > m){
    
    fa <- factanal(
      R_DS,
      factors = k,
      scores = "regression"
    )
    
    init_u <- as.matrix(fa$scores)
    
    init_gamma <- as.matrix(
      fa$loadings
    )[1:m,1:k,drop=FALSE]
    
  } else if(n < m){
    
    fa <- factanal(
      t(R_DS),
      factors = k,
      scores = "regression"
    )
    
    init_u <- as.matrix(
      fa$loadings
    )[1:n,1:k,drop=FALSE]
    
    init_gamma <- as.matrix(
      fa$scores
    )[1:m,1:k,drop=FALSE]
    
  } else {
    
    fa <- factanal(
      rbind(
        R_DS,
        rnorm(m,0,0.01)
      ),
      factors = k,
      scores = "regression"
    )
    
    init_u <- as.matrix(
      fa$scores
    )[1:n,1:k,drop=FALSE]
    
    init_gamma <- as.matrix(
      fa$loadings
    )[1:m,1:k,drop=FALSE]
  }
  
  
  ##########################################################
  ## ======================================================
  ## Step 4:
  ## QR rotation
  ## ======================================================
  ##########################################################
  
  if(k > 1 && m > 2){
    
    qr_gamma <- qr(
      t(init_gamma)
    )
    
    gamma_new <- t(
      qr.R(qr_gamma)
    )
    
    Qmat <- qr.Q(qr_gamma)
    
    sig <- sign(
      diag(gamma_new)
    )
    
    sig[sig == 0] <- 1
    
    init_gamma <- t(
      t(gamma_new) * sig
    )
    
    init_u <- init_u %*% Qmat
    
    init_u <- t(
      t(init_u) * sig
    )
    
  } else {
    
    sig <- sign(init_gamma[1,1])
    
    if(sig == 0) sig <- 1
    
    init_gamma <- init_gamma * sig
    
    init_u <- init_u * sig
  }
  
  
  ##########################################################
  ## ======================================================
  ## Step 5:
  ## scaling
  ## ======================================================
  ##########################################################
  
  sdu <- sqrt(
    diag(cov(init_u))
  )
  
  sdg <- sqrt(
    diag(cov(init_gamma))
  )
  
  sdu[is.na(sdu) | sdu == 0] <- 1
  
  sdg[is.na(sdg) | sdg == 0] <- 1
  
  init_u <- init_u %*%
    diag(1/sdu, nrow=k)
  
  init_gamma <- init_gamma %*%
    diag(0.8/sdg, nrow=k)
  
  
  ##########################################################
  ## ======================================================
  ## Step 6:
  ## optional jitter
  ## ======================================================
  ##########################################################
  
  if(jitter){
    
    init_u <- init_u +
      MASS::mvrnorm(
        n = n,
        mu = rep(0,k),
        Sigma = diag(k) * jitter_sd^2
      )
  }
  
  
  ##########################################################
  ## return
  ##########################################################
  
  out <- list(
    
    model = model,
    
    route = route,
    
    R_DS = R_DS,
    
    init_u = init_u,
    
    init_gamma = init_gamma,
    
    fa = fa
  )
  
  if(model == "ZIDM"){
    
    out$gamma_fit <- gamma_fit
    out$gamma_plus <- gamma_plus
  }
  
  return(out)
}

fit_ZR <- function(Y){
  
  Y <- as.matrix(Y)
  
  
  ##########################################################
  ## zero replacement
  ##########################################################
  
  Y_zr <- Y
  
  Y_zr[Y_zr == 0] <- 0.5
  
  
  ##########################################################
  ## composition
  ##########################################################
  
  rho_zr <- Y_zr /
    rowSums(Y_zr)
  
  
  ##########################################################
  ## return unified format
  ##########################################################
  
  return(list(
    
    rho_hat = rho_zr,
    
    ord_hat = NULL,
    
    beta0_hat = NULL,
    
    beta_hat = NULL,
    
    gamma_hat = NULL,
    
    zeta_hat = NULL,
    
    convergence = TRUE,
    
    raw_fit = NULL
  ))
}

evaluate_composition <- function(rho_hat, rho_true, eps = 1e-12){
  
  # -----------------------------
  # basic checks
  # -----------------------------
  
  if(is.null(rho_hat) || !is.matrix(rho_hat)){
    return(list(
      Frobenius = NA,
      KL = NA,
      Shannon_MSE = NA,
      Simpson_MSE = NA
    ))
  }
  
  if(any(is.na(rho_hat))){
    return(list(
      Frobenius = NA,
      KL = NA,
      Shannon_MSE = NA,
      Simpson_MSE = NA
    ))
  }
  
  if(!all(dim(rho_hat) == dim(rho_true))){
    return(list(
      Frobenius = NA,
      KL = NA,
      Shannon_MSE = NA,
      Simpson_MSE = NA
    ))
  }
  
  # -----------------------------
  # avoid log(0)
  # -----------------------------
  
  rho_hat  <- pmax(rho_hat, eps)
  rho_true <- pmax(rho_true, eps)
  
  # renormalize (important for safety)
  rho_hat  <- rho_hat / rowSums(rho_hat)
  rho_true <- rho_true / rowSums(rho_true)
  
  # -----------------------------
  # 1 Frobenius norm
  # -----------------------------
  
  frob <- sqrt(sum((rho_hat - rho_true)^2))
  
  # -----------------------------
  # 2 KL divergence
  # -----------------------------
  
  KL <- mean(rowSums(rho_true * log(rho_true / rho_hat)))
  
  # -----------------------------
  # 3 Shannon index MSE
  # -----------------------------
  
  shannon_true <- -rowSums(rho_true * log(rho_true))
  shannon_hat  <- -rowSums(rho_hat  * log(rho_hat))
  
  shannon_mse <- mean((shannon_hat - shannon_true)^2)
  
  # -----------------------------
  # 4 Simpson index MSE
  # -----------------------------
  
  simpson_true <- 1 - rowSums(rho_true^2)
  simpson_hat  <- 1 - rowSums(rho_hat^2)
  
  simpson_mse <- mean((simpson_hat - simpson_true)^2)
  
  return(list(
    Frobenius = frob,
    KL = KL,
    Shannon_MSE = shannon_mse,
    Simpson_MSE = simpson_mse
  ))
}

compute_procrustes <- function(true_ord, est_ord){
  
  # -----------------------------
  # basic checks
  # -----------------------------
  
  if(is.null(est_ord) || is.null(true_ord)) return(NA)
  
  if(!is.matrix(est_ord) || !is.matrix(true_ord)) return(NA)
  
  if(nrow(est_ord) != nrow(true_ord)) return(NA)
  
  if(any(is.na(est_ord))) return(NA)
  
  # -----------------------------
  # scale ordinations
  # -----------------------------
  
  true_ord <- scale(true_ord)
  est_ord  <- scale(est_ord)
  
  # -----------------------------
  # procrustes
  # -----------------------------
  
  p <- try(vegan::procrustes(true_ord, est_ord), silent = TRUE)
  
  if(inherits(p, "try-error")) return(NA)
  
  return(p$ss)
}

evaluate_parameter <- function(est, true){
  
  ##########################################################
  ## basic checks
  ##########################################################
  
  if(is.null(est) || is.null(true)){
    return(list(
      Bias = NA,
      AbsBias = NA,
      RMSE = NA
    ))
  }
  
  est  <- as.numeric(est)
  true <- as.numeric(true)
  
  if(length(est) != length(true)){
    return(list(
      Bias = NA,
      AbsBias = NA,
      RMSE = NA
    ))
  }
  
  if(
    any(!is.finite(est)) ||
    any(!is.finite(true))
  ){
    return(list(
      Bias = NA,
      AbsBias = NA,
      RMSE = NA
    ))
  }
  
  
  ##########################################################
  ## estimation error
  ##########################################################
  
  error <- est - true
  
  
  ##########################################################
  ## Bias
  ##########################################################
  
  bias <- mean(error)
  
  
  ##########################################################
  ## Absolute Bias
  ##########################################################
  
  abs_bias <- mean(abs(error))
  
  
  ##########################################################
  ## RMSE
  ##########################################################
  
  rmse <- sqrt(
    mean(error^2)
  )
  
  
  ##########################################################
  ## return
  ##########################################################
  
  return(list(
    Bias = bias,
    AbsBias = abs_bias,
    RMSE = rmse
  ))
}

fit_ZIDM_paper <- function(
    Y,
    X = NULL,
    X_theta = NULL,
    iterations = 10000,
    thin = 10,
    burnin_prop = 0.5,
    sigma2_beta_gamma = sqrt(5),
    sigma2_beta_theta = sqrt(5),
    seed = 123
){
  
  set.seed(seed)
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  m <- ncol(Y)
  
  
  # ============================================================================
  # 1. Abundance design matrix
  #
  # Add intercept ONCE
  # ============================================================================
  
  if(is.null(X)){
    
    X_design <- matrix(
      1,
      nrow = n,
      ncol = 1
    )
    
    colnames(X_design) <- "Intercept"
    
  } else {
    
    X <- as.matrix(X)
    
    if(nrow(X) != n){
      stop("X must have the same number of rows as Y.")
    }
    
    X_design <- cbind(
      Intercept = 1,
      X
    )
  }
  
  
  # ============================================================================
  # 2. Zero-inflation design matrix
  #
  # Your current simulation:
  # zeta_j does NOT depend on x1/x2
  #
  # therefore intercept only
  # ============================================================================
  
  if(is.null(X_theta)){
    
    X_theta_design <- matrix(
      1,
      nrow = n,
      ncol = 1
    )
    
    colnames(X_theta_design) <- "Intercept"
    
  } else {
    
    X_theta <- as.matrix(X_theta)
    
    if(nrow(X_theta) != n){
      stop("X_theta must have same number of rows as Y.")
    }
    
    X_theta_design <- cbind(
      Intercept = 1,
      X_theta
    )
  }
  
  
  q <- ncol(X_design)
  
  q_theta <- ncol(X_theta_design)
  
  
  # ============================================================================
  # 3. Number of stored MCMC samples
  # ============================================================================
  
  samples <- floor(
    iterations / thin
  )
  
  if(samples < 2){
    stop("Too few stored MCMC samples.")
  }
  
  
  # ============================================================================
  # 4. Initialize beta_gamma
  #
  # taxa x covariates x samples
  # ============================================================================
  
  beta_gamma <- array(
    0,
    dim = c(
      m,
      q,
      samples
    )
  )
  
  beta_gamma[, , 1] <- matrix(
    rnorm(m * q),
    nrow = m,
    ncol = q
  )
  
  
  # ============================================================================
  # 5. Initialize latent gamma variables c_ij
  # ============================================================================
  
  cc <- array(
    0,
    dim = c(
      n,
      m,
      samples
    )
  )
  
  cc[, , 1] <- Y
  
  
  # ============================================================================
  # 6. Auxiliary u_i
  #
  # follows official ZIDM_R source
  # ============================================================================
  
  temp_uu <- rep(
    0,
    n
  )
  
  for(i in seq_len(n)){
    
    sum_Y <- sum(
      Y[i, ]
    )
    
    if(sum_Y > 0){
      
      temp_uu[i] <- rgamma(
        1,
        shape = sum_Y,
        rate = sum_Y
      )
      
    } else {
      
      temp_uu[i] <- 1
    }
  }
  
  
  # ============================================================================
  # 7. Structural-zero indicators eta
  # ============================================================================
  
  eta <- array(
    0,
    dim = c(
      n,
      m,
      samples
    )
  )
  
  eta[, , 1] <- 1
  
  
  eta_accept <- array(
    0,
    dim = c(
      n,
      m,
      samples
    )
  )
  
  
  eta_index <- vector(
    "list",
    m
  )
  
  
  for(j in seq_len(m)){
    
    zero_id <- which(
      Y[, j] == 0
    )
    
    # C++ uses zero-based indexing
    eta_index[[j]] <-
      zero_id - 1
    
    
    if(length(zero_id) > 0){
      
      eta[
        zero_id,
        j,
        1
      ] <- rbinom(
        length(zero_id),
        size = 1,
        prob = 0.5
      )
    }
  }
  
  
  # ============================================================================
  # 8. beta_theta
  #
  # taxa x ZI covariates x samples
  # ============================================================================
  
  beta_theta <- array(
    0,
    dim = c(
      m,
      q_theta,
      samples
    )
  )
  
  beta_theta[, , 1] <- matrix(
    rnorm(m * q_theta),
    nrow = m,
    ncol = q_theta
  )
  
  
  # ============================================================================
  # 9. Polya-Gamma auxiliary omega
  # ============================================================================
  
  omega <- array(
    0,
    dim = c(
      n,
      m,
      samples
    )
  )
  
  omega[, , 1] <- 1
  
  
  # ============================================================================
  # 10. Run OFFICIAL ZIDM C++ sampler
  # ============================================================================
  
  t0 <- Sys.time()
  
  
  MCMC <- ZIDM::zidm(
    
    iterations = iterations,
    
    thin = thin,
    
    z = Y,
    
    x = X_design,
    
    x_theta = X_theta_design,
    
    beta_gamma = beta_gamma,
    
    cc = cc,
    
    temp_uu = temp_uu,
    
    eta = eta,
    
    eta_accept = eta_accept,
    
    eta_index = eta_index,
    
    omega = omega,
    
    beta_theta = beta_theta,
    
    zeta = array(
      1,
      dim = c(
        m,
        q_theta,
        samples
      )
    ),
    
    varphi = array(
      1,
      dim = c(
        m,
        q,
        samples
      )
    ),
    
    sigma2_beta_gamma =
      sigma2_beta_gamma,
    
    sigma2_beta_theta =
      sigma2_beta_theta
  )
  
  
  runtime <- as.numeric(
    Sys.time() - t0,
    units = "secs"
  )
  
  
  # ============================================================================
  # 11. Official output ordering
  # ============================================================================
  
  raw_fit <- MCMC[1:8]
  
  names(raw_fit) <- c(
    "varphi",
    "beta_gamma",
    "eta",
    "beta_theta",
    "zeta",
    "omega",
    "eta_accept",
    "cc"
  )
  
  
  # ============================================================================
  # 12. Posterior summary
  # ============================================================================
  
  summary_fit <- summarize_ZIDM_paper(
    zidm_obj = raw_fit,
    X = X_design,
    burnin_prop = burnin_prop
  )
  
  
  # ============================================================================
  # 13. Unified output
  # ============================================================================
  
  return(
    c(
      summary_fit,
      
      list(
        
        raw_fit =
          raw_fit,
        
        X_design =
          X_design,
        
        X_theta_design =
          X_theta_design,
        
        runtime =
          runtime,
        
        iterations =
          iterations,
        
        thin =
          thin,
        
        burnin_prop =
          burnin_prop,
        
        convergence =
          TRUE,
        ######################################################
        ## no latent ordination in Koslovsky ZIDM
        ######################################################
        
        ordination = NULL,
        
        loadings = NULL,
        model =
          "ZIDM_paper"
      )
    )
  )
}

summarize_ZIDM_paper <- function(
    zidm_obj,
    X,
    burnin_prop = 0.5,
    CI_level = 0.95
){
  
  beta_gamma <- zidm_obj$beta_gamma
  beta_theta <- zidm_obj$beta_theta
  cc <- zidm_obj$cc
  
  
  # ============================================================================
  # dimensions
  # ============================================================================
  
  m <- dim(beta_gamma)[1]
  
  q <- dim(beta_gamma)[2]
  
  S <- dim(beta_gamma)[3]
  
  n <- dim(cc)[1]
  
  
  # ============================================================================
  # burn-in
  #
  # IMPORTANT:
  # use STORED sample index, not 0.5 directly
  # ============================================================================
  
  burn <- floor(
    S * burnin_prop
  )
  
  keep <- seq.int(
    burn + 1,
    S
  )
  
  
  cat(
    "Stored samples =",
    S,
    "| burn-in =",
    burn,
    "| retained =",
    length(keep),
    "\n"
  )
  
  
  # ============================================================================
  # 1. Regression coefficient posterior means
  # ============================================================================
  
  beta_gamma_hat <- apply(
    beta_gamma[, , keep, drop = FALSE],
    c(1, 2),
    mean,
    na.rm = TRUE
  )
  
  
  beta0_hat <-
    beta_gamma_hat[, 1]
  
  
  beta_hat <- if(q > 1){
    
    beta_gamma_hat[
      ,
      2:q,
      drop = FALSE
    ]
    
  } else {
    
    matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================================
  # 2. Posterior CIs for regression coefficients
  # ============================================================================
  
  alpha <- (1 - CI_level) / 2
  
  
  beta_gamma_lower <- apply(
    beta_gamma[, , keep, drop = FALSE],
    c(1,2),
    quantile,
    probs = alpha,
    na.rm = TRUE
  )
  
  
  beta_gamma_upper <- apply(
    beta_gamma[, , keep, drop = FALSE],
    c(1,2),
    quantile,
    probs = 1 - alpha,
    na.rm = TRUE
  )
  
  
  # ============================================================================
  # 3. Structural-zero probability
  #
  # Source package's estimates_ZIDM uses:
  #
  # 1 / (1 + exp(beta_theta))
  #
  # This is structural-zero probability.
  #
  # Current simulation uses intercept-only zero inflation.
  # ============================================================================
  
  q_theta <- dim(beta_theta)[2]
  
  
  if(q_theta == 1){
    
    zeta_draws <-
      1 / (
        1 +
          exp(
            beta_theta[
              ,
              1,
              keep,
              drop = FALSE
            ]
          )
      )
    
    
    zeta_hat <- apply(
      zeta_draws,
      1,
      mean,
      na.rm = TRUE
    )
    
    
    zeta_lower <- apply(
      zeta_draws,
      1,
      quantile,
      probs = alpha,
      na.rm = TRUE
    )
    
    
    zeta_upper <- apply(
      zeta_draws,
      1,
      quantile,
      probs = 1 - alpha,
      na.rm = TRUE
    )
    
  } else {
    
    zeta_hat <- NULL
    zeta_lower <- NULL
    zeta_upper <- NULL
  }
  
  
  # ============================================================================
  # 4. Baseline composition rho0
  #
  # For every retained MCMC draw:
  #
  # eta_ij = X_i beta_gamma_j
  #
  # rho0 = softmax(eta)
  #
  # This is BEFORE structural-zero realization.
  # ============================================================================
  
  rho0_sum <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  
  for(s in keep){
    
    B_s <- beta_gamma[
      ,
      ,
      s,
      drop = FALSE
    ]
    
    B_s <- matrix(
      B_s,
      nrow = m,
      ncol = q
    )
    
    
    eta_s <-
      X %*%
      t(B_s)
    
    
    rho_s <- matrix(
      0,
      nrow = n,
      ncol = m
    )
    
    
    for(i in seq_len(n)){
      
      temp <- exp(
        eta_s[i, ] -
          max(
            eta_s[i, ]
          )
      )
      
      rho_s[i, ] <-
        temp / sum(temp)
    }
    
    
    rho0_sum <-
      rho0_sum +
      rho_s
  }
  
  
  rho0_hat <-
    rho0_sum /
    length(keep)
  
  
  # ============================================================================
  # 5. Individual zero-inflated composition psi
  #
  # psi_ij = c_ij / sum_j c_ij
  # ============================================================================
  
  psi_sum <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  
  valid_draws <- 0
  
  
  for(s in keep){
    
    cc_s <- cc[
      ,
      ,
      s,
      drop = FALSE
    ]
    
    cc_s <- matrix(
      cc_s,
      nrow = n,
      ncol = m
    )
    
    
    rs <- rowSums(
      cc_s
    )
    
    
    valid_rows <- rs > 0
    
    
    psi_s <- matrix(
      0,
      nrow = n,
      ncol = m
    )
    
    
    psi_s[
      valid_rows,
    ] <-
      cc_s[
        valid_rows,
        ,
        drop = FALSE
      ] /
      rs[
        valid_rows
      ]
    
    
    psi_sum <-
      psi_sum +
      psi_s
    
    
    valid_draws <-
      valid_draws + 1
  }
  
  
  psi_hat <-
    psi_sum /
    valid_draws
  
  
  # ============================================================================
  # return
  # ============================================================================
  
  return(
    list(
      
      # ------------------------------------------------------------------------
      # parameter recovery
      # ------------------------------------------------------------------------
      
      beta0_hat =
        beta0_hat,
      
      beta_hat =
        beta_hat,
      
      beta_gamma_hat =
        beta_gamma_hat,
      
      beta_gamma_lower =
        beta_gamma_lower,
      
      beta_gamma_upper =
        beta_gamma_upper,
      
      
      # ------------------------------------------------------------------------
      # zero inflation
      # ------------------------------------------------------------------------
      
      zeta_hat =
        zeta_hat,
      
      zeta_lower =
        zeta_lower,
      
      zeta_upper =
        zeta_upper,
      
      
      # ------------------------------------------------------------------------
      # composition
      #
      # rho_hat = baseline composition
      # ------------------------------------------------------------------------
      
      rho_hat =
        rho0_hat,
      
      rho0_hat =
        rho0_hat,
      
      rhoZI_hat =
        psi_hat,
      
      psi_hat =
        psi_hat,
      
      
      # ------------------------------------------------------------------------
      # diagnostics
      # ------------------------------------------------------------------------
      
      burnin =
        burn,
      
      retained_samples =
        length(keep)
    )
  )
}

compute_TMB_beta_CI <- function(
    fit,
    beta0_true = NULL,
    beta1_true = NULL,
    beta2_true = NULL,
    level = 0.95
){
  
  ##########################################################
  ## safety
  ##########################################################
  
  if(is.null(fit)){
    return(NULL)
  }
  
  if(is.null(fit$fit)){
    return(NULL)
  }
  
  if(is.null(fit$fit$obj) ||
     is.null(fit$fit$opt)){
    return(NULL)
  }
  
  
  ##########################################################
  ## TMB object and optimum
  ##########################################################
  
  obj <- fit$fit$obj
  
  par_opt <- fit$fit$opt$par
  
  
  ##########################################################
  ## Hessian at optimum
  ##########################################################
  
  H_full <- obj$he(par_opt)
  
  H_full <- as.matrix(H_full)
  
  
  if(
    any(!is.finite(H_full))
  ){
    return(NULL)
  }
  
  
  ##########################################################
  ## locate b parameters
  ##
  ## b contains:
  ## row 1 = beta0
  ## row 2 = beta1
  ## row 3 = beta2
  ##########################################################
  
  par_names <- names(par_opt)
  
  idx_beta <- which(
    par_names == "b"
  )
  
  
  if(length(idx_beta) == 0){
    stop("Could not identify parameter 'b' in TMB parameter vector.")
  }
  
  
  idx_other <- setdiff(
    seq_along(par_opt),
    idx_beta
  )
  
  
  ##########################################################
  ## Hessian blocks
  ##########################################################
  
  H_bb <- H_full[
    idx_beta,
    idx_beta,
    drop = FALSE
  ]
  
  
  ##########################################################
  ## profile Hessian
  ##########################################################
  
  if(length(idx_other) > 0){
    
    H_bo <- H_full[
      idx_beta,
      idx_other,
      drop = FALSE
    ]
    
    H_oo <- H_full[
      idx_other,
      idx_other,
      drop = FALSE
    ]
    
    H_profile <-
      H_bb -
      H_bo %*%
      MASS::ginv(H_oo) %*%
      t(H_bo)
    
  } else {
    
    H_profile <- H_bb
  }
  
  
  ##########################################################
  ## covariance matrix for b
  ##########################################################
  
  V_beta <- MASS::ginv(
    H_profile
  )
  
  
  if(
    any(!is.finite(V_beta))
  ){
    return(NULL)
  }
  
  
  ##########################################################
  ## standard errors
  ##########################################################
  
  variance_beta <- diag(
    V_beta
  )
  
  if(any(variance_beta < 0)){
    warning(
      "Negative variance detected in TMB beta covariance."
    )
  }
  
  se_vec <- sqrt(
    pmax(
      variance_beta,
      0
    )
  )
  
  
  ##########################################################
  ## recover fitted b matrix
  ##########################################################
  
  b_hat <- fit$b_hat
  
  nr <- nrow(b_hat)
  nc <- ncol(b_hat)
  
  
  ##########################################################
  ## restore matrix shape
  ##########################################################
  
  se_mat <- matrix(
    se_vec,
    nrow = nr,
    ncol = nc
  )
  
  
  ##########################################################
  ## confidence intervals
  ##########################################################
  
  zcrit <- qnorm(
    1 - (1 - level) / 2
  )
  
  lower <- b_hat -
    zcrit * se_mat
  
  upper <- b_hat +
    zcrit * se_mat
  
  
  ##########################################################
  ## create taxon-level table
  ##########################################################
  
  out <- data.frame(
    
    taxon = seq_len(nc),
    
    beta0_hat = b_hat[1, ],
    beta0_SE = se_mat[1, ],
    beta0_lower = lower[1, ],
    beta0_upper = upper[1, ],
    
    stringsAsFactors = FALSE
  )
  
  
  ##########################################################
  ## beta1
  ##########################################################
  
  if(nr >= 2){
    
    out$beta1_hat <- b_hat[2, ]
    
    out$beta1_SE <- se_mat[2, ]
    
    out$beta1_lower <- lower[2, ]
    
    out$beta1_upper <- upper[2, ]
    
  } else {
    
    out$beta1_hat <- NA_real_
    out$beta1_SE <- NA_real_
    out$beta1_lower <- NA_real_
    out$beta1_upper <- NA_real_
  }
  
  
  ##########################################################
  ## beta2
  ##########################################################
  
  if(nr >= 3){
    
    out$beta2_hat <- b_hat[3, ]
    
    out$beta2_SE <- se_mat[3, ]
    
    out$beta2_lower <- lower[3, ]
    
    out$beta2_upper <- upper[3, ]
    
  } else {
    
    out$beta2_hat <- NA_real_
    out$beta2_SE <- NA_real_
    out$beta2_lower <- NA_real_
    out$beta2_upper <- NA_real_
  }
  
  
  ##########################################################
  ## TRUE values + coverage
  ##########################################################
  
  if(!is.null(beta0_true)){
    
    out$beta0_true <- beta0_true
    
    out$beta0_cover <-
      beta0_true >= out$beta0_lower &
      beta0_true <= out$beta0_upper
  }
  
  
  if(!is.null(beta1_true)){
    
    out$beta1_true <- beta1_true
    
    out$beta1_cover <-
      beta1_true >= out$beta1_lower &
      beta1_true <= out$beta1_upper
  }
  
  
  if(!is.null(beta2_true)){
    
    out$beta2_true <- beta2_true
    
    out$beta2_cover <-
      beta2_true >= out$beta2_lower &
      beta2_true <= out$beta2_upper
  }
  
  
  ##########################################################
  ## return
  ##########################################################
  
  return(list(
    
    table = out,
    
    b_hat = b_hat,
    
    SE = se_mat,
    
    lower = lower,
    
    upper = upper,
    
    V_beta = V_beta,
    
    H_profile = H_profile,
    
    H_full = H_full
  ))
}

compute_gllvm_ZINB_beta_CI <- function(
    fit,
    beta0_true = NULL,
    beta1_true = NULL,
    beta2_true = NULL,
    level = 0.95
){

  if(is.null(fit) || is.null(fit$model)){
    return(NULL)
  }

  model <- fit$model
  ci_warning <- NULL

  ##########################################################
  ## 1. confint.gllvm() needs stored standard errors.
  ##    Keep the original fit unchanged and calculate these
  ##    only in this local copy of the model.
  ##########################################################

  se_obj <- tryCatch(
    withCallingHandlers(
      gllvm::se(model),
      warning = function(w){
        msg <- conditionMessage(w)
        ci_warning <<- paste(c(ci_warning, msg), collapse = " | ")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if(is.null(se_obj) || is.null(se_obj$sd)){
    return(NULL)
  }

  model$sd <- se_obj$sd
  model$Hess <- se_obj$Hess

  ##########################################################
  ## 2. Built-in gllvm confidence intervals
  ##########################################################

  ci_beta0 <- tryCatch(
    withCallingHandlers(
      stats::confint(
        model,
        parm = "beta0",
        level = level
      ),
      warning = function(w){
        msg <- conditionMessage(w)
        ci_warning <<- paste(c(ci_warning, msg), collapse = " | ")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  ci_Xcoef <- tryCatch(
    withCallingHandlers(
      stats::confint(
        model,
        parm = "Xcoef",
        level = level
      ),
      warning = function(w){
        msg <- conditionMessage(w)
        ci_warning <<- paste(c(ci_warning, msg), collapse = " | ")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if(is.null(ci_beta0) || is.null(ci_Xcoef)){
    return(NULL)
  }

  ci_beta0 <- as.matrix(ci_beta0)
  ci_Xcoef <- as.matrix(ci_Xcoef)

  ##########################################################
  ## 3. Point estimates and dimensions
  ##########################################################

  beta0_hat <- as.numeric(fit$beta0_hat)
  beta_hat <- as.matrix(fit$beta_hat)

  m <- length(beta0_hat)

  if(
    length(dim(beta_hat)) != 2 ||
    nrow(beta_hat) != m ||
    ncol(beta_hat) < 2
  ){
    return(NULL)
  }

  q <- ncol(beta_hat)

  if(
    ncol(ci_beta0) < 2 ||
    nrow(ci_beta0) != m ||
    ncol(ci_Xcoef) < 2 ||
    nrow(ci_Xcoef) != m * q
  ){
    return(NULL)
  }

  ##########################################################
  ## 4. Restore Xcoef CI to the same m x q orientation as
  ##    fit$params$Xcoef. confint.gllvm() flattens Xcoef in
  ##    R column-major order, so matrix(..., nrow = m) restores
  ##    the original species x covariate layout.
  ##########################################################

  beta0_lower <- as.numeric(ci_beta0[, 1])
  beta0_upper <- as.numeric(ci_beta0[, 2])

  beta_lower <- matrix(
    as.numeric(ci_Xcoef[, 1]),
    nrow = m,
    ncol = q
  )

  beta_upper <- matrix(
    as.numeric(ci_Xcoef[, 2]),
    nrow = m,
    ncol = q
  )

  ##########################################################
  ## 5. Coverage validity: CI bounds only
  ##########################################################

  valid_beta0 <-
    is.finite(beta0_lower) &
    is.finite(beta0_upper) &
    beta0_lower < beta0_upper

  valid_beta1 <-
    is.finite(beta_lower[, 1]) &
    is.finite(beta_upper[, 1]) &
    beta_lower[, 1] < beta_upper[, 1]

  valid_beta2 <-
    is.finite(beta_lower[, 2]) &
    is.finite(beta_upper[, 2]) &
    beta_lower[, 2] < beta_upper[, 2]

  ##########################################################
  ## 6. Compatibility SE columns
  ##
  ## summarize_beta_CI() in the rest of this file is left
  ## unchanged. These values are reconstructed from the CI
  ## width only so that its existing interface still works.
  ## They are NOT used to construct the gllvm intervals.
  ##########################################################

  zcrit <- qnorm(1 - (1 - level) / 2)

  beta0_se <-
    (beta0_upper - beta0_lower) /
    (2 * zcrit)

  beta_se <-
    (beta_upper - beta_lower) /
    (2 * zcrit)

  out <- data.frame(
    taxon = seq_len(m),

    beta0_hat = beta0_hat,
    beta0_SE = beta0_se,
    beta0_lower = beta0_lower,
    beta0_upper = beta0_upper,
    beta0_valid = valid_beta0,

    beta1_hat = beta_hat[, 1],
    beta1_SE = beta_se[, 1],
    beta1_lower = beta_lower[, 1],
    beta1_upper = beta_upper[, 1],
    beta1_valid = valid_beta1,

    beta2_hat = beta_hat[, 2],
    beta2_SE = beta_se[, 2],
    beta2_lower = beta_lower[, 2],
    beta2_upper = beta_upper[, 2],
    beta2_valid = valid_beta2,

    stringsAsFactors = FALSE
  )

  ##########################################################
  ## 7. Coverage
  ##########################################################

  if(!is.null(beta0_true)){
    beta0_true <- as.numeric(beta0_true)
    out$beta0_true <- beta0_true
    out$beta0_cover <- NA
    out$beta0_cover[valid_beta0] <-
      beta0_true[valid_beta0] >= beta0_lower[valid_beta0] &
      beta0_true[valid_beta0] <= beta0_upper[valid_beta0]
  }

  if(!is.null(beta1_true)){
    beta1_true <- as.numeric(beta1_true)
    out$beta1_true <- beta1_true
    out$beta1_cover <- NA
    out$beta1_cover[valid_beta1] <-
      beta1_true[valid_beta1] >= beta_lower[valid_beta1, 1] &
      beta1_true[valid_beta1] <= beta_upper[valid_beta1, 1]
  }

  if(!is.null(beta2_true)){
    beta2_true <- as.numeric(beta2_true)
    out$beta2_true <- beta2_true
    out$beta2_cover <- NA
    out$beta2_cover[valid_beta2] <-
      beta2_true[valid_beta2] >= beta_lower[valid_beta2, 2] &
      beta2_true[valid_beta2] <= beta_upper[valid_beta2, 2]
  }

  return(list(
    table = out,
    SE_beta0 = beta0_se,
    SE_beta = beta_se,
    lower_beta0 = beta0_lower,
    upper_beta0 = beta0_upper,
    lower_beta = beta_lower,
    upper_beta = beta_upper,
    Hessian = se_obj$Hess,
    ci_warning = ci_warning,
    confint_beta0 = ci_beta0,
    confint_Xcoef = ci_Xcoef,
    beta0_nValid = sum(valid_beta0),
    beta1_nValid = sum(valid_beta1),
    beta2_nValid = sum(valid_beta2),
    beta0_validProp = mean(valid_beta0),
    beta1_validProp = mean(valid_beta1),
    beta2_validProp = mean(valid_beta2)
  ))
}

compute_ZIDM_paper_beta_CI <- function(
    fit,
    beta0_true = NULL,
    beta1_true = NULL,
    beta2_true = NULL,
    level = 0.95
){
  
  if(
    is.null(fit) ||
    is.null(fit$raw_fit) ||
    is.null(fit$raw_fit$beta_gamma)
  ){
    return(NULL)
  }
  
  
  beta_draws <- fit$raw_fit$beta_gamma
  
  # taxa x coefficients x stored MCMC samples
  m <- dim(beta_draws)[1]
  q <- dim(beta_draws)[2]
  S <- dim(beta_draws)[3]
  
  
  ############################################################
  ## burn-in
  ############################################################
  
  burn <- floor(
    S * fit$burnin_prop
  )
  
  keep <- seq.int(
    burn + 1,
    S
  )
  
  
  ############################################################
  ## posterior probabilities
  ############################################################
  
  alpha <- (1 - level) / 2
  
  
  ############################################################
  ## posterior means
  ############################################################
  
  beta_mean <- apply(
    beta_draws[, , keep, drop = FALSE],
    c(1, 2),
    mean,
    na.rm = TRUE
  )
  
  
  ############################################################
  ## posterior SD
  ############################################################
  
  beta_sd <- apply(
    beta_draws[, , keep, drop = FALSE],
    c(1, 2),
    sd,
    na.rm = TRUE
  )
  
  
  ############################################################
  ## credible intervals
  ############################################################
  
  beta_lower <- apply(
    beta_draws[, , keep, drop = FALSE],
    c(1, 2),
    quantile,
    probs = alpha,
    na.rm = TRUE
  )
  
  
  beta_upper <- apply(
    beta_draws[, , keep, drop = FALSE],
    c(1, 2),
    quantile,
    probs = 1 - alpha,
    na.rm = TRUE
  )
  
  
  ############################################################
  ## need intercept + x1 + x2
  ############################################################
  
  if(q < 3){
    
    stop(
      "ZIDM_paper requires intercept + two covariates ",
      "for the current simulation."
    )
  }
  
  
  ############################################################
  ## output
  ############################################################
  
  out <- data.frame(
    
    taxon = seq_len(m),
    
    # beta0
    beta0_hat = beta_mean[, 1],
    beta0_SE = beta_sd[, 1],
    beta0_lower = beta_lower[, 1],
    beta0_upper = beta_upper[, 1],
    
    # beta1
    beta1_hat = beta_mean[, 2],
    beta1_SE = beta_sd[, 2],
    beta1_lower = beta_lower[, 2],
    beta1_upper = beta_upper[, 2],
    
    # beta2
    beta2_hat = beta_mean[, 3],
    beta2_SE = beta_sd[, 3],
    beta2_lower = beta_lower[, 3],
    beta2_upper = beta_upper[, 3],
    
    stringsAsFactors = FALSE
  )
  
  
  ############################################################
  ## coverage
  ############################################################
  
  if(!is.null(beta0_true)){
    
    beta0_true <- as.numeric(beta0_true)
    
    out$beta0_true <- beta0_true
    
    out$beta0_cover <-
      beta0_true >= out$beta0_lower &
      beta0_true <= out$beta0_upper
  }
  
  
  if(!is.null(beta1_true)){
    
    beta1_true <- as.numeric(beta1_true)
    
    out$beta1_true <- beta1_true
    
    out$beta1_cover <-
      beta1_true >= out$beta1_lower &
      beta1_true <= out$beta1_upper
  }
  
  
  if(!is.null(beta2_true)){
    
    beta2_true <- as.numeric(beta2_true)
    
    out$beta2_true <- beta2_true
    
    out$beta2_cover <-
      beta2_true >= out$beta2_lower &
      beta2_true <= out$beta2_upper
  }
  
  
  return(
    list(
      table = out,
      posterior_mean = beta_mean,
      posterior_sd = beta_sd,
      lower = beta_lower,
      upper = beta_upper
    )
  )
}

summarize_beta_CI <- function(ci){

  fail <- list(
    ci_success = FALSE,
    beta0_Coverage = NA_real_,
    beta1_Coverage = NA_real_,
    beta2_Coverage = NA_real_,
    beta0_meanSE = NA_real_,
    beta1_meanSE = NA_real_,
    beta2_meanSE = NA_real_,
    beta0_meanWidth = NA_real_,
    beta1_meanWidth = NA_real_,
    beta2_meanWidth = NA_real_,
    beta0_nValid = NA_integer_,
    beta1_nValid = NA_integer_,
    beta2_nValid = NA_integer_,
    beta0_validProp = NA_real_,
    beta1_validProp = NA_real_,
    beta2_validProp = NA_real_
  )

  if(is.null(ci) || is.null(ci$table)){
    return(fail)
  }

  tab <- ci$table

  get_valid <- function(prefix){
    se_col <- paste0(prefix, "_SE")
    lo_col <- paste0(prefix, "_lower")
    up_col <- paste0(prefix, "_upper")
    valid_col <- paste0(prefix, "_valid")

    if(!all(c(se_col, lo_col, up_col) %in% names(tab))){
      return(rep(FALSE, nrow(tab)))
    }

    vv <- is.finite(tab[[se_col]]) &
      tab[[se_col]] > 0 &
      is.finite(tab[[lo_col]]) &
      is.finite(tab[[up_col]]) &
      tab[[lo_col]] < tab[[up_col]]

    if(valid_col %in% names(tab)){
      vv <- vv & (tab[[valid_col]] %in% TRUE)
    }

    vv
  }

  summarize_one <- function(prefix){
    valid <- get_valid(prefix)
    cover_col <- paste0(prefix, "_cover")
    se_col <- paste0(prefix, "_SE")
    lo_col <- paste0(prefix, "_lower")
    up_col <- paste0(prefix, "_upper")

    n_valid <- sum(valid)
    prop_valid <- if(nrow(tab) > 0) mean(valid) else NA_real_

    if(n_valid == 0){
      return(list(
        Coverage = NA_real_, meanSE = NA_real_, meanWidth = NA_real_,
        nValid = 0L, validProp = prop_valid
      ))
    }

    cov_val <- NA_real_
    if(cover_col %in% names(tab)){
      cc <- tab[[cover_col]][valid]
      if(any(!is.na(cc))){
        cov_val <- mean(cc, na.rm = TRUE)
      }
    }

    list(
      Coverage = cov_val,
      meanSE = mean(tab[[se_col]][valid]),
      meanWidth = mean(tab[[up_col]][valid] - tab[[lo_col]][valid]),
      nValid = as.integer(n_valid),
      validProp = prop_valid
    )
  }

  s0 <- summarize_one("beta0")
  s1 <- summarize_one("beta1")
  s2 <- summarize_one("beta2")

  list(
    ci_success = any(c(s0$nValid, s1$nValid, s2$nValid) > 0),

    beta0_Coverage = s0$Coverage,
    beta1_Coverage = s1$Coverage,
    beta2_Coverage = s2$Coverage,

    beta0_meanSE = s0$meanSE,
    beta1_meanSE = s1$meanSE,
    beta2_meanSE = s2$meanSE,

    beta0_meanWidth = s0$meanWidth,
    beta1_meanWidth = s1$meanWidth,
    beta2_meanWidth = s2$meanWidth,

    beta0_nValid = s0$nValid,
    beta1_nValid = s1$nValid,
    beta2_nValid = s2$nValid,

    beta0_validProp = s0$validProp,
    beta1_validProp = s1$validProp,
    beta2_validProp = s2$validProp
  )
}


# ==================== ZIDM generator ====================
run_simulation_ZIDM <- local({
  # Generator-specific implementation: generate_ZIDM_data
generate_ZIDM_data <- function(
    
  # ============================================================
  # Basic simulation dimensions
  # ============================================================
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.30)
  # high   : U(0.35, 0.45)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional user-specified true values
  # ============================================================
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter distributions from simulation table
  # ============================================================
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenario settings
  # If supplied, override zi_level/depth_level
  # ============================================================
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Structural-zero handling
  # ============================================================
  zero_inflation = TRUE,
  
  # "zero_row":
  # if all taxa structurally absent, Y_i = 0
  #
  # "resample":
  # redraw Z until at least one taxon survives
  all_zero_action = c("zero_row", "resample"),
  
  max_resample = 1000,
  
  # ============================================================
  # Output control
  # ============================================================
  return_Z = TRUE,
  return_lambda = TRUE,
  return_eta = TRUE,
  return_alpha_ij = TRUE,
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  all_zero_action <- match.arg(all_zero_action)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Scenario definitions
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.30),
      
      high = c(0.35, 0.45)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be a positive number.")
  }
  
  
  # ============================================================
  # 2. Generate observed covariates X
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 3. Generate latent variables
  #
  # u_i ~ N(0, I_k)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 4. Generate beta0
  #
  # beta_0j ~ U(-1,1)
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Generate covariate effects
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Generate factor loadings
  #
  # gamma_j ~ U(-2,2)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Zero-inflation probabilities
  #
  # zeta_j ~ Uniform(...)
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 8. Linear predictor
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 9. Dirichlet/Gamma shape
  #
  # alpha_ij = exp(eta_ij)
  # ============================================================
  
  alpha_ij <- exp(eta)
  
  
  if (
    any(!is.finite(alpha_ij)) ||
    any(alpha_ij <= 0)
  ) {
    stop(
      "Non-finite or non-positive exp(eta) encountered."
    )
  }
  
  
  # ============================================================
  # 10. Gamma latent variables
  #
  # lambda0_ij ~ Gamma(exp(eta_ij), rate = 1)
  # ============================================================
  
  lambda0 <- matrix(
    
    rgamma(
      n * m,
      shape = as.vector(alpha_ij),
      rate = 1
    ),
    
    nrow = n,
    ncol = m
  )
  
  
  # ============================================================
  # 11. Composition BEFORE zero inflation
  #
  # theta0_ij =
  # lambda0_ij / sum_l lambda0_il
  # ============================================================
  
  lambda0_rowsum <- rowSums(lambda0)
  
  rho0_true <- lambda0 /
    lambda0_rowsum
  
  
  # ============================================================
  # 12. Structural-zero process
  #
  # zeta_j = structural-zero probability
  #
  # Therefore:
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  #
  # Z = 1 : taxon survives
  # Z = 0 : structural zero
  # ============================================================
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      rep(
        1 - zeta,
        each = n
      ),
      nrow = n,
      ncol = m
    )
    
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
    
  } else {
    
    Z <- matrix(
      1L,
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 13. Optional handling of rows with all structural zeros
  # ============================================================
  
  if (
    zero_inflation &&
    all_zero_action == "resample"
  ) {
    
    zero_rows <- which(
      rowSums(Z) == 0
    )
    
    iter <- 0
    
    while (length(zero_rows) > 0) {
      
      iter <- iter + 1
      
      if (iter > max_resample) {
        stop(
          "Exceeded max_resample for all-zero structural rows."
        )
      }
      
      for (i in zero_rows) {
        
        Z[i, ] <- rbinom(
          m,
          size = 1,
          prob = 1 - zeta
        )
      }
      
      zero_rows <- which(
        rowSums(Z) == 0
      )
    }
  }
  
  
  # ============================================================
  # 14. Apply structural zeros
  # ============================================================
  
  lambdaZI <- lambda0 * Z
  
  row_mass <- rowSums(lambdaZI)
  
  
  # ============================================================
  # 15. ZI composition
  # ============================================================
  
  rhoZI_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  
  positive_rows <- row_mass > 0
  
  
  rhoZI_true[positive_rows, ] <-
    lambdaZI[positive_rows, , drop = FALSE] /
    row_mass[positive_rows]
  
  
  # ============================================================
  # 16. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.integer(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  
  N <- as.integer(N)
  
  
  # ============================================================
  # 17. Multinomial observation
  #
  # If all taxa are structurally zero:
  #
  # Y_i = 0
  #
  # otherwise:
  #
  # Y_i ~ Multinomial(N_i, rhoZI_i)
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  
  for (i in seq_len(n)) {
    
    if (row_mass[i] > 0) {
      
      Y[i, ] <- as.integer(
        
        rmultinom(
          n = 1,
          size = N[i],
          prob = rhoZI_true[i, ]
        )
      )
    }
  }
  
  
  # ============================================================
  # 18. Useful simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  all_structural_zero_prop <-
    mean(rowSums(Z) == 0)
  
  
  # ============================================================
  # 19. Return EVERYTHING needed for parameter recovery
  # ============================================================
  
  out <- list(
    
    # --------------------------------------------------------
    # observed data
    # --------------------------------------------------------
    
    Y = Y,
    X = X,
    N = N,
    
    
    # --------------------------------------------------------
    # TRUE composition
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for model recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # TRUE latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # TRUE zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # true linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    alpha_ij = alpha_ij,
    
    
    # --------------------------------------------------------
    # scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    all_structural_zero_prop =
      all_structural_zero_prop
  )
  
  
  # ============================================================
  # Optional latent outputs
  # ============================================================
  
  if (return_Z) {
    
    out$Z <- Z
  }
  
  
  if (return_lambda) {
    
    out$lambda0 <- lambda0
    
    out$lambdaZI <- lambdaZI
  }
  
  
  if (!return_eta) {
    
    out$eta <- NULL
  }
  
  
  if (!return_alpha_ij) {
    
    out$alpha_ij <- NULL
  }
  
  
  return(out)
}

  # Generator-specific implementation: generate_ZIMN_data
generate_ZIMN_data <- function(
    
  # ============================================================
  # Basic dimensions
  # ============================================================
  
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.40)
  # high   : U(0.50, 0.70)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional true parameter values
  #
  # IMPORTANT:
  # parameters are now defined directly in TAXON SPACE
  # ============================================================
  
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter ranges
  # ============================================================
  
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenarios
  # ============================================================
  
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Zero inflation
  # ============================================================
  
  zero_inflation = TRUE,
  
  # ============================================================
  # Random seed
  # ============================================================
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 1 || m != as.integer(m)) {
    stop("m must be an integer greater than 1.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Simulation scenario settings
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.40),
      
      high = c(0.50, 0.70)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be positive.")
  }
  
  
  # ============================================================
  # 2. Latent variables
  #
  # u_i ~ N(0, I_k)
  #
  # U: n x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 3. Observed covariates
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  #
  # X: n x 2
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 4. Species-specific intercepts
  #
  # beta_0j ~ U(-1,1)
  #
  # beta0: length m
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Covariate coefficients
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # beta: m x 2
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Latent loadings
  #
  # gamma_j ~ U(-2,2)
  #
  # Gamma: m x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Linear predictor DIRECTLY in taxon space
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  #
  # eta: n x m
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  #
  # X:        n x 2
  # t(beta):  2 x m
  #
  # result:   n x m
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  #
  # U:         n x k
  # t(Gamma):  k x m
  #
  # result:    n x m
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 8. TRUE baseline composition
  #
  # rho_ij =
  #
  # exp(eta_ij)
  # ------------------
  # sum_l exp(eta_il)
  #
  # stable softmax
  # ============================================================
  
  rho0_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    row_max <- max(
      eta[i, ]
    )
    
    tmp <- exp(
      eta[i, ] - row_max
    )
    
    rho0_true[i, ] <-
      tmp / sum(tmp)
  }
  
  
  # ============================================================
  # 9. Zero-inflation probabilities
  #
  # zeta_j = P(structural zero)
  #
  # zeta_j ~ Uniform(...)
  #
  # zeta: length m
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 10. Structural-zero indicators
  #
  # zeta_j = P(structural zero)
  #
  # Define:
  #
  # Z_ij = 1 : taxon is PRESENT
  # Z_ij = 0 : structural zero
  #
  # therefore
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  # ============================================================
  
  Z <- matrix(
    1L,
    nrow = n,
    ncol = m
  )
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      
      rep(
        1 - zeta,
        each = n
      ),
      
      nrow = n,
      ncol = m
    )
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 11. Apply structural zeros
  #
  # rhoZI_ij proportional to
  #
  # Z_ij * rho0_ij
  # ============================================================
  
  rhoZI_raw <-
    rho0_true * Z
  
  rs <- rowSums(
    rhoZI_raw
  )
  
  
  # ------------------------------------------------------------
  # Extremely rare case:
  # all taxa removed in one sample
  #
  # Activate one taxon so multinomial probabilities are valid
  # ------------------------------------------------------------
  
  bad <- which(
    rs == 0
  )
  
  if (length(bad) > 0) {
    
    for (i in bad) {
      
      j <- sample(
        seq_len(m),
        size = 1
      )
      
      Z[i, j] <- 1L
      
      rhoZI_raw[i, j] <-
        rho0_true[i, j]
    }
    
    rs <- rowSums(
      rhoZI_raw
    )
  }
  
  
  # ============================================================
  # 12. TRUE zero-inflated composition
  #
  # rhoZI_ij =
  #
  # Z_ij rho0_ij
  # --------------------
  # sum_l Z_il rho0_il
  # ============================================================
  
  rhoZI_true <-
    rhoZI_raw / rs
  
  
  # ============================================================
  # 13. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.numeric(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  N <- as.integer(
    round(N)
  )
  
  if (any(N < 0)) {
    stop("N must be nonnegative.")
  }
  
  
  # ============================================================
  # 14. Multinomial sampling
  #
  # Y_i ~ Multinomial(
  #     N_i,
  #     rhoZI_true_i
  # )
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    Y[i, ] <- as.integer(
      
      rmultinom(
        n = 1,
        size = N[i],
        prob = rhoZI_true[i, ]
      )
    )
  }
  
  
  # ============================================================
  # 15. Simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  
  # ============================================================
  # 16. Return EVERYTHING needed for simulation evaluation
  # ============================================================
  
  return(list(
    
    # --------------------------------------------------------
    # Observed data
    # --------------------------------------------------------
    
    Y = Y,
    
    X = X,
    
    N = N,
    
    
    # --------------------------------------------------------
    # True compositions
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for composition recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    P_true = rho0_true,
    
    P_ZI = rhoZI_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    #
    # NOW ALL PARAMETERS ARE IN TAXON SPACE
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # True latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    Beta = t(beta),
    
    Lambda = Gamma,
    
    
    # --------------------------------------------------------
    # True zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # Linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    eta_taxon = eta,
    
    
    # --------------------------------------------------------
    # Scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    
    # --------------------------------------------------------
    # Structural-zero indicator
    # --------------------------------------------------------
    
    Z = Z
  ))
}

  # Generator-specific implementation: fit_gllvm_ZINB
fit_gllvm_ZINB <- function(
    Y,
    X = NULL,
    k = 2,
    seed = 123
){
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  m <- ncol(Y)
  
  if(!is.null(X)){
    
    X <- as.matrix(X)
    
    if(nrow(X) != n){
      stop("X must have the same number of rows as Y.")
    }
  }
  
  
  ##########################################################
  ## 1. Fit ZINB-GLLVM
  ##
  ## Improved convergence:
  ## - more iterations
  ## - multiple initial values
  ## - residual-based starting values
  ##########################################################
  
  fit <- gllvm::gllvm(
    y = Y,
    X = X,
    family = "ZINB",
    num.lv = k,
    row.eff = "fixed",
    method = "VA",
    sd.errors = FALSE,
    seed = seed,
    
    control = list(
      maxit = 20000,
      max.iter = 20000,
      reltol = 1e-10
    ),
    
    control.start = list(
      starting.val = "res",
      n.init = 3
    )
  )
  
  
  ##########################################################
  ## 2. Extract fitted model parameters
  ##########################################################
  
  # species-specific intercepts beta_0j
  beta0_hat <- fit$params$beta0
  
  # environmental coefficients beta_j
  beta_hat <- fit$params$Xcoef
  
  # latent scores u_i
  u_hat <- as.matrix(
    gllvm::getLV(fit)
  )
  
  # properly scaled loading vectors gamma_j
  gamma_hat <- as.matrix(
    gllvm::getLoadings(fit)
  )
  
  # fixed sample / row effects alpha_i
  # sample 1 is the reference, alpha_1 = 0
  alpha_hat <- c(
    0,
    fit$params$row.params.fixed
  )
  
  
  ##########################################################
  ## 3. Zero-inflation and dispersion parameters
  ##########################################################
  
  # ZI probability p_j
  # This corresponds to zeta_j in our notation
  zeta_hat <- fit$params$phi
  
  # NB dispersion parameter
  dispersion_hat <- fit$params$ZINB.phi
  
  # inverse NB dispersion
  inv_dispersion_hat <- fit$params$ZINB.inv.phi
  
  
  ##########################################################
  ## 4. Variational covariance matrices A_i
  ##
  ## A_hat[i,,] = A_i
  ## dimension n x k x k
  ##########################################################
  
  A_hat <- fit$A
  
  
  ##########################################################
  ## 5. Construct fitted linear predictor eta
  ##
  ## eta_ij =
  ## alpha_i
  ## + beta_0j
  ## + x_i^T beta_j
  ## + u_i^T gamma_j
  ##########################################################
  
  eta_hat <- matrix(
    beta0_hat,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  ##########################################################
  ## add row effect alpha_i
  ##########################################################
  
  eta_hat <- eta_hat +
    matrix(
      alpha_hat,
      nrow = n,
      ncol = m,
      byrow = FALSE
    )
  
  
  ##########################################################
  ## add environmental contribution X beta
  ##########################################################
  
  if(!is.null(X) && !is.null(beta_hat)){
    
    beta_mat <- as.matrix(beta_hat)
    
    # usual gllvm orientation: m x q
    if(
      nrow(beta_mat) == m &&
      ncol(beta_mat) == ncol(X)
    ){
      
      eta_hat <- eta_hat +
        X %*% t(beta_mat)
      
    } else if(
      nrow(beta_mat) == ncol(X) &&
      ncol(beta_mat) == m
    ){
      
      eta_hat <- eta_hat +
        X %*% beta_mat
      
    } else {
      
      stop(
        "Dimension of fit$params$Xcoef is not compatible with X."
      )
    }
  }
  
  
  ##########################################################
  ## add latent contribution U Gamma
  ##########################################################
  
  # convert loadings to k x m
  if(
    nrow(gamma_hat) == m &&
    ncol(gamma_hat) == k
  ){
    
    Gamma_mat <- t(gamma_hat)
    
  } else if(
    nrow(gamma_hat) == k &&
    ncol(gamma_hat) == m
  ){
    
    Gamma_mat <- gamma_hat
    
  } else {
    
    stop(
      "Dimension of getLoadings(fit) is not compatible with k and m."
    )
  }
  
  eta_hat <- eta_hat +
    u_hat %*% Gamma_mat
  
  
  ##########################################################
  ## 6. Variational quadratic correction
  ##
  ## cQ_ij =
  ## 0.5 * gamma_j^T A_i gamma_j
  ##########################################################
  
  cQ_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    Ai <- A_hat[i, , ]
    
    for(j in seq_len(m)){
      
      gamma_j <- Gamma_mat[, j]
      
      cQ_hat[i, j] <-
        0.5 *
        as.numeric(
          t(gamma_j) %*%
            Ai %*%
            gamma_j
        )
    }
  }
  
  
  ##########################################################
  ## 7. Variationally corrected predictor
  ##########################################################
  
  eta_corrected_hat <-
    eta_hat + cQ_hat
  
  
  ##########################################################
  ## 8. Corrected fitted composition rho
  ##
  ## rho_ij =
  ## exp(eta_ij + cQ_ij)
  ## -------------------------
  ## sum_l exp(eta_il + cQ_il)
  ##########################################################
  
  rho_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    row_max <- max(
      eta_corrected_hat[i, ]
    )
    
    tmp <- exp(
      eta_corrected_hat[i, ] -
        row_max
    )
    
    rho_hat[i, ] <-
      tmp / sum(tmp)
  }
  
  
  ##########################################################
  ## 9. Return everything needed for simulation
  ##########################################################
  
  return(list(
    
    ########################################################
    ## regression / fitted model parameters
    ########################################################
    
    alpha_hat = alpha_hat,
    
    beta0_hat = beta0_hat,
    
    beta_hat = beta_hat,
    
    
    ########################################################
    ## latent variables
    ########################################################
    
    u_hat = u_hat,
    
    gamma_hat = gamma_hat,
    
    
    ########################################################
    ## zero-inflation parameters
    ########################################################
    
    zeta_hat = zeta_hat,
    
    dispersion_hat = dispersion_hat,
    
    inv_dispersion_hat = inv_dispersion_hat,
    
    
    ########################################################
    ## variational quantities
    ########################################################
    
    A_hat = A_hat,
    
    cQ_hat = cQ_hat,
    
    
    ########################################################
    ## predictors
    ########################################################
    
    eta_hat = eta_hat,
    
    eta_corrected_hat = eta_corrected_hat,
    
    
    ########################################################
    ## fitted composition
    ########################################################
    
    rho_hat = rho_hat,
    
    
    ########################################################
    ## raw loading parameters
    ########################################################
    
    theta_raw = fit$params$theta,
    
    sigma_lv = fit$params$sigma.lv,
    
    
    ########################################################
    ## convenient aliases
    ########################################################
    
    ordination = u_hat,
    
    loadings = gamma_hat,
    
    
    ########################################################
    ## complete fitted object
    ########################################################
    
    model = fit
  ))
}

run_simulation_ZIDM <- function(
    rep_k = 500,  # rep_k = 100
    n = 50,
    m = 100,
    k = 2,
    zi_level = c("low", "medium", "high"),
    depth_level = c("low", "medium", "high"),
    init_method = "DS",
    seed_base = 2026,
    compute_CI = TRUE,
    CI_level = 0.95,
    include_ZIDM_paper = TRUE,
    ZIDM_paper_iterations = 20000,
    ZIDM_paper_thin = 10,
    ZIDM_paper_burnin = 0.5
){

  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)

  ## FINAL calibrated sequencing-depth distribution.
  ## No dataset filtering is applied.
  if(n == 50 && m == 100){
    calibrated_depth_range <- switch(
      zi_level,
      low = c(320L, 480L),
      medium = c(800L, 1200L),
      high = c(10000L, 15000L)
    )
  } else if(n == 100 && m == 50){
    calibrated_depth_range <- switch(
      zi_level,
      low = c(140L, 210L),
      medium = c(320L, 480L),
      high = c(4000L, 6000L)
    )
  } else {
    stop("No calibrated ZIDM depth range for this (n,m) setting.")
  }

  tmb_model_list <- list(
    list(name = "ZIMN_soft", model = "ZIMN", extra = c(0, 0)),
    list(name = "ZIMN_hard", model = "ZIMN", extra = c(0, 1)),
    list(name = "ZIDM_soft", model = "ZIDM", extra = c(0, 0)),
    list(name = "ZIDM_hard", model = "ZIDM", extra = c(0, 1))
  )

  methods_all <- c(
    "ZINB_gllvm",
    "ZIMN_soft",
    "ZIMN_hard",
    "ZIDM_soft",
    "ZIDM_hard"
  )
  if(include_ZIDM_paper){
    methods_all <- c(methods_all, "ZIDM_paper")
  }
  methods_all <- c(methods_all, "ZR")

  ## ------------------------------------------------------------
  ## Helpers
  ## ------------------------------------------------------------

  valid_TMB <- function(fit){
    isTRUE(
      !is.null(fit) &&
      !is.null(fit$convergence) &&
      length(fit$convergence) == 1L &&
      is.finite(fit$convergence) &&
      fit$convergence == 0 &&
      !is.null(fit$rho_hat) &&
      all(is.finite(fit$rho_hat))
    )
  }

  valid_gllvm <- function(fit){
    isTRUE(
      !is.null(fit) &&
      !is.null(fit$model) &&
      !is.null(fit$model$convergence) &&
      isTRUE(fit$model$convergence) &&
      !is.null(fit$rho_hat) &&
      all(is.finite(fit$rho_hat))
    )
  }

  valid_ZIDM_paper <- function(fit){
    isTRUE(
      !is.null(fit) &&
      !is.null(fit$rho_hat) &&
      is.matrix(fit$rho_hat) &&
      all(dim(fit$rho_hat) == c(n, m)) &&
      all(is.finite(fit$rho_hat)) &&
      !is.null(fit$beta0_hat) &&
      length(fit$beta0_hat) == m &&
      all(is.finite(fit$beta0_hat)) &&
      !is.null(fit$beta_hat) &&
      all(dim(fit$beta_hat) == c(m, 2)) &&
      all(is.finite(fit$beta_hat)) &&
      !is.null(fit$zeta_hat) &&
      length(fit$zeta_hat) == m &&
      all(is.finite(fit$zeta_hat))
    )
  }

  empty_ci_summary <- function(requested = compute_CI){
    list(
      ci_success = if(requested) FALSE else NA,
      beta0_Coverage = NA_real_,
      beta1_Coverage = NA_real_,
      beta2_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta1_meanSE = NA_real_,
      beta2_meanSE = NA_real_,
      beta0_meanWidth = NA_real_,
      beta1_meanWidth = NA_real_,
      beta2_meanWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta1_nValid = NA_integer_,
      beta2_nValid = NA_integer_,
      beta0_validProp = NA_real_,
      beta1_validProp = NA_real_,
      beta2_validProp = NA_real_
    )
  }

  make_na_row <- function(
      r,
      method,
      runtime = NA_real_,
      convergence = NA,
      ci_runtime = NA_real_,
      observed_zero_prop = NA_real_,
      structural_zero_prop = NA_real_,
      sampling_zero_prop = NA_real_,
      generation_runtime = NA_real_
  ){
    data.frame(
      rep = r,
      method = method,
      success = FALSE,
      convergence = convergence,
      runtime = runtime,
      ci_runtime = ci_runtime,
      runtime_total =
        if(is.finite(runtime) && is.finite(ci_runtime)) runtime + ci_runtime
        else if(is.finite(runtime)) runtime
        else NA_real_,
      ci_success = if(compute_CI && method != "ZR") FALSE else NA,

      Frobenius = NA_real_,
      KL = NA_real_,
      Shannon_MSE = NA_real_,
      Simpson_MSE = NA_real_,
      Procrustes = NA_real_,

      beta0_Bias = NA_real_,
      beta0_AbsBias = NA_real_,
      beta0_RMSE = NA_real_,
      beta1_Bias = NA_real_,
      beta1_AbsBias = NA_real_,
      beta1_RMSE = NA_real_,
      beta2_Bias = NA_real_,
      beta2_AbsBias = NA_real_,
      beta2_RMSE = NA_real_,
      zeta_Bias = NA_real_,
      zeta_AbsBias = NA_real_,
      zeta_RMSE = NA_real_,

      ## IMPORTANT: failed rows contain EXACTLY the same CI columns
      ## as successful rows. This prevents rbind column mismatch.
      beta0_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta0_meanCIWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta0_validProp = NA_real_,

      beta1_Coverage = NA_real_,
      beta1_meanSE = NA_real_,
      beta1_meanCIWidth = NA_real_,
      beta1_nValid = NA_integer_,
      beta1_validProp = NA_real_,

      beta2_Coverage = NA_real_,
      beta2_meanSE = NA_real_,
      beta2_meanCIWidth = NA_real_,
      beta2_nValid = NA_integer_,
      beta2_validProp = NA_real_,

      observed_zero_prop = observed_zero_prop,
      structural_zero_prop = structural_zero_prop,
      sampling_zero_prop = sampling_zero_prop,
      generation_runtime = generation_runtime,
      stringsAsFactors = FALSE
    )
  }

  ## Extra safety: even if a future method row gains/loses a column,
  ## bind by column name instead of crashing the whole replicate.
  rbind_fill <- function(x){
    x <- x[!vapply(x, is.null, logical(1))]
    if(length(x) == 0L) return(NULL)

    all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))

    x <- lapply(x, function(z){
      miss <- setdiff(all_cols, names(z))
      if(length(miss) > 0L){
        for(nm in miss) z[[nm]] <- NA
      }
      z[, all_cols, drop = FALSE]
    })

    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
  }

  timed_fit <- function(expr, validator, type, convergence_fun){
    t0 <- Sys.time()
    err <- NULL

    fit <- tryCatch(
      force(expr),
      error = function(e){
        err <<- conditionMessage(e)
        NULL
      }
    )

    runtime <- as.numeric(Sys.time() - t0, units = "secs")

    ok <- tryCatch(
      isTRUE(validator(fit)),
      error = function(e){
        err <<- paste0("validation: ", conditionMessage(e))
        FALSE
      }
    )

    conv <- tryCatch(
      convergence_fun(fit),
      error = function(e) NA
    )

    list(
      fit = if(ok) fit else NULL,
      runtime = runtime,
      success = ok,
      convergence = conv,
      type = type,
      error = err
    )
  }

  results <- vector("list", rep_k)

  ## ------------------------------------------------------------
  ## Replicates
  ## ------------------------------------------------------------

  for(r in seq_len(rep_k)){

    cat(
      "\n========================================\n",
      "ZIDM simulation replicate: ", r, " / ", rep_k, "\n",
      "n = ", n, " | m = ", m,
      " | zi = ", zi_level,
      " | depth = ", depth_level, "\n",
      "========================================\n",
      sep = ""
    )

    rep_seed <- seed_base + r

    ## ==========================================================
    ## 1. DATA GENERATION
    ## Only a data-generation failure is allowed to make all
    ## methods NA, because no method can be fitted without data.
    ## ==========================================================

    t0_generate <- Sys.time()

    ## N_i ~ discrete Uniform(L,U), independently for each sample.
    set.seed(rep_seed + 9000000L)
    N_now <- sample(
      seq.int(calibrated_depth_range[1], calibrated_depth_range[2]),
      size = n,
      replace = TRUE
    )

    dat <- tryCatch(
      generate_ZIDM_data(
        n = n,
        m = m,
        k = k,
        zi_level = zi_level,
        depth_level = depth_level,
        N = N_now,
        all_zero_action = "resample",
        zero_inflation = TRUE,
        seed = rep_seed
      ),
      error = function(e){
        cat(
          "  DATA GENERATION FAILED: ",
          conditionMessage(e),
          "\n",
          sep = ""
        )
        NULL
      }
    )

    generation_runtime <- as.numeric(
      Sys.time() - t0_generate,
      units = "secs"
    )

    if(is.null(dat)){
      results[[r]] <- rbind_fill(
        lapply(
          methods_all,
          function(mm){
            make_na_row(
              r = r,
              method = mm,
              generation_runtime = generation_runtime
            )
          }
        )
      )
      next
    }

    Y <- dat$Y
    X <- dat$X
    rho_true <- dat$rho_true
    U_true <- dat$U
    beta0_true <- dat$beta0
    beta1_true <- dat$beta1
    beta2_true <- dat$beta2
    zeta_true <- dat$zeta

    observed_zero_prop <- dat$observed_zero_prop
    structural_zero_prop <- dat$structural_zero_prop
    sampling_zero_prop <- dat$sampling_zero_prop

    fit_results <- list()

    ## ==========================================================
    ## 2. FITS
    ## Every method is protected independently.
    ## ==========================================================

    cat("  fitting ZINB_gllvm ...\n")

    fit_results$ZINB_gllvm <- timed_fit(
      fit_gllvm_ZINB(
        Y = Y,
        X = X,
        k = k,
        seed = rep_seed
      ),
      validator = valid_gllvm,
      type = "gllvm",
      convergence_fun = function(fit){
        if(is.null(fit) || is.null(fit$model$convergence)) return(FALSE)
        isTRUE(fit$model$convergence)
      }
    )

    for(mod in tmb_model_list){

      cat("  fitting ", mod$name, " ...\n", sep = "")

      fit_results[[mod$name]] <- timed_fit(
        fit_TMB(
          Y = Y,
          X = X,
          k = k,
          model = mod$model,
          extra = mod$extra,
          init_method = init_method
        ),
        validator = valid_TMB,
        type = "TMB",
        convergence_fun = function(fit){
          if(is.null(fit) || is.null(fit$convergence)) return(NA)
          fit$convergence
        }
      )
    }

    if(include_ZIDM_paper){

      cat("  fitting ZIDM_paper ...\n")

      fit_results$ZIDM_paper <- timed_fit(
        fit_ZIDM_paper(
          Y = Y,
          X = X,
          X_theta = NULL,
          iterations = ZIDM_paper_iterations,
          thin = ZIDM_paper_thin,
          burnin_prop = ZIDM_paper_burnin,
          seed = rep_seed
        ),
        validator = valid_ZIDM_paper,
        type = "ZIDM_paper",
        convergence_fun = function(fit){
          !is.null(fit)
        }
      )
    }

    cat("  fitting ZR ...\n")

    fit_results$ZR <- timed_fit(
      fit_ZR(Y),
      validator = function(fit){
        isTRUE(
          !is.null(fit) &&
          !is.null(fit$rho_hat) &&
          all(is.finite(fit$rho_hat))
        )
      },
      type = "ZR",
      convergence_fun = function(fit){
        !is.null(fit)
      }
    )

    ## ==========================================================
    ## 3. CI
    ## CI failure NEVER changes fit success.
    ## ==========================================================

    ci_results <- list()

    for(method in methods_all){

      requested_ci <- compute_CI && method != "ZR"

      ci_results[[method]] <- list(
        ci = NULL,
        summary = empty_ci_summary(requested_ci),
        runtime = NA_real_
      )

      if(!requested_ci) next

      fit_obj <- fit_results[[method]]

      if(
        is.null(fit_obj) ||
        !isTRUE(fit_obj$success) ||
        is.null(fit_obj$fit)
      ){
        next
      }

      t0_ci <- Sys.time()

      ci_tmp <- tryCatch(
        {
          if(fit_obj$type == "gllvm"){

            compute_gllvm_ZINB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )

          } else if(fit_obj$type == "TMB"){

            compute_TMB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )

          } else if(fit_obj$type == "ZIDM_paper"){

            compute_ZIDM_paper_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )

          } else {
            NULL
          }
        },
        error = function(e){
          cat(
            "  CI FAILED [", method, "]: ",
            conditionMessage(e),
            "\n",
            sep = ""
          )
          NULL
        }
      )

      runtime_ci <- as.numeric(
        Sys.time() - t0_ci,
        units = "secs"
      )

      ci_sum <- if(is.null(ci_tmp)){
        empty_ci_summary(TRUE)
      } else {
        tryCatch(
          summarize_beta_CI(ci_tmp),
          error = function(e){
            cat(
              "  CI SUMMARY FAILED [", method, "]: ",
              conditionMessage(e),
              "\n",
              sep = ""
            )
            empty_ci_summary(TRUE)
          }
        )
      }

      ci_results[[method]] <- list(
        ci = ci_tmp,
        summary = ci_sum,
        runtime = runtime_ci
      )
    }

    ## ==========================================================
    ## 4. EVALUATION
    ## Every method row is built independently. If anything in
    ## one method's evaluation fails, ONLY that method is NA.
    ## ==========================================================

    out_rep <- vector("list", length(methods_all))
    names(out_rep) <- methods_all

    for(method in methods_all){

      obj <- fit_results[[method]]

      if(
        is.null(obj) ||
        !isTRUE(obj$success) ||
        is.null(obj$fit)
      ){
        out_rep[[method]] <- make_na_row(
          r = r,
          method = method,
          runtime =
            if(!is.null(obj) && !is.null(obj$runtime)) obj$runtime
            else NA_real_,
          convergence =
            if(!is.null(obj) && !is.null(obj$convergence)) obj$convergence
            else NA,
          ci_runtime =
            if(!is.null(ci_results[[method]])) ci_results[[method]]$runtime
            else NA_real_,
          observed_zero_prop = observed_zero_prop,
          structural_zero_prop = structural_zero_prop,
          sampling_zero_prop = sampling_zero_prop,
          generation_runtime = generation_runtime
        )

        if(!is.null(obj$error) && !is.null(obj$error)){
          cat(
            "  METHOD FAILED [", method, "]",
            if(!is.null(obj$error)) paste0(": ", obj$error) else "",
            "\n",
            sep = ""
          )
        }

        next
      }

      out_rep[[method]] <- tryCatch(
        {
          fit <- obj$fit

          ci_obj <- ci_results[[method]]
          if(is.null(ci_obj)){
            ci_obj <- list(
              ci = NULL,
              summary = empty_ci_summary(compute_CI && method != "ZR"),
              runtime = NA_real_
            )
          }

          ci_sum <- ci_obj$summary
          if(is.null(ci_sum)){
            ci_sum <- empty_ci_summary(compute_CI && method != "ZR")
          }

          ci_runtime <- ci_obj$runtime

          runtime_total <- if(
            is.finite(obj$runtime) &&
            is.finite(ci_runtime)
          ){
            obj$runtime + ci_runtime
          } else {
            obj$runtime
          }

          ## Composition
          comp <- tryCatch(
            evaluate_composition(
              rho_hat = fit$rho_hat,
              rho_true = rho_true
            ),
            error = function(e) NULL
          )

          if(is.null(comp)){
            comp <- list(
              Frobenius = NA_real_,
              KL = NA_real_,
              Shannon_MSE = NA_real_,
              Simpson_MSE = NA_real_
            )
          }

          ## Ordination
          proc <- NA_real_

          if(method != "ZR" && method != "ZIDM_paper"){
            proc_tmp <- tryCatch(
              compute_procrustes(
                true_ord = U_true,
                est_ord = fit$ordination
              ),
              error = function(e) NA_real_
            )

            if(
              length(proc_tmp) == 1L &&
              is.finite(proc_tmp)
            ){
              proc <- proc_tmp
            }
          }

          ## Parameter metrics
          empty_met <- list(
            Bias = NA_real_,
            AbsBias = NA_real_,
            RMSE = NA_real_
          )

          beta0_met <- empty_met
          beta1_met <- empty_met
          beta2_met <- empty_met
          zeta_met <- empty_met

          if(method != "ZR"){

            tmp <- tryCatch(
              evaluate_parameter(
                est = fit$beta0_hat,
                true = beta0_true
              ),
              error = function(e) NULL
            )
            if(!is.null(tmp)) beta0_met <- tmp

            if(obj$type %in% c("gllvm", "ZIDM_paper")){

              tmp <- tryCatch(
                evaluate_parameter(
                  est = fit$beta_hat[, 1],
                  true = beta1_true
                ),
                error = function(e) NULL
              )
              if(!is.null(tmp)) beta1_met <- tmp

              tmp <- tryCatch(
                evaluate_parameter(
                  est = fit$beta_hat[, 2],
                  true = beta2_true
                ),
                error = function(e) NULL
              )
              if(!is.null(tmp)) beta2_met <- tmp

            } else if(obj$type == "TMB"){

              tmp <- tryCatch(
                evaluate_parameter(
                  est = fit$beta_hat[1, ],
                  true = beta1_true
                ),
                error = function(e) NULL
              )
              if(!is.null(tmp)) beta1_met <- tmp

              tmp <- tryCatch(
                evaluate_parameter(
                  est = fit$beta_hat[2, ],
                  true = beta2_true
                ),
                error = function(e) NULL
              )
              if(!is.null(tmp)) beta2_met <- tmp
            }

            tmp <- tryCatch(
              evaluate_parameter(
                est = fit$zeta_hat,
                true = zeta_true
              ),
              error = function(e) NULL
            )
            if(!is.null(tmp)) zeta_met <- tmp
          }

          ## Successful method row
          data.frame(
            rep = r,
            method = method,
            success = TRUE,
            convergence = obj$convergence,
            runtime = obj$runtime,
            ci_runtime = ci_runtime,
            runtime_total = runtime_total,
            ci_success = ci_sum$ci_success,

            Frobenius = comp$Frobenius,
            KL = comp$KL,
            Shannon_MSE = comp$Shannon_MSE,
            Simpson_MSE = comp$Simpson_MSE,
            Procrustes = proc,

            beta0_Bias = beta0_met$Bias,
            beta0_AbsBias = beta0_met$AbsBias,
            beta0_RMSE = beta0_met$RMSE,

            beta1_Bias = beta1_met$Bias,
            beta1_AbsBias = beta1_met$AbsBias,
            beta1_RMSE = beta1_met$RMSE,

            beta2_Bias = beta2_met$Bias,
            beta2_AbsBias = beta2_met$AbsBias,
            beta2_RMSE = beta2_met$RMSE,

            zeta_Bias = zeta_met$Bias,
            zeta_AbsBias = zeta_met$AbsBias,
            zeta_RMSE = zeta_met$RMSE,

            beta0_Coverage = ci_sum$beta0_Coverage,
            beta0_meanSE = ci_sum$beta0_meanSE,
            beta0_meanCIWidth = ci_sum$beta0_meanWidth,
            beta0_nValid = ci_sum$beta0_nValid,
            beta0_validProp = ci_sum$beta0_validProp,

            beta1_Coverage = ci_sum$beta1_Coverage,
            beta1_meanSE = ci_sum$beta1_meanSE,
            beta1_meanCIWidth = ci_sum$beta1_meanWidth,
            beta1_nValid = ci_sum$beta1_nValid,
            beta1_validProp = ci_sum$beta1_validProp,

            beta2_Coverage = ci_sum$beta2_Coverage,
            beta2_meanSE = ci_sum$beta2_meanSE,
            beta2_meanCIWidth = ci_sum$beta2_meanWidth,
            beta2_nValid = ci_sum$beta2_nValid,
            beta2_validProp = ci_sum$beta2_validProp,

            observed_zero_prop = observed_zero_prop,
            structural_zero_prop = structural_zero_prop,
            sampling_zero_prop = sampling_zero_prop,
            generation_runtime = generation_runtime,
            stringsAsFactors = FALSE
          )
        },
        error = function(e){

          cat(
            "  EVALUATION FAILED [", method, "]: ",
            conditionMessage(e),
            " -> only this method set to NA\n",
            sep = ""
          )

          make_na_row(
            r = r,
            method = method,
            runtime = obj$runtime,
            convergence = obj$convergence,
            ci_runtime =
              if(!is.null(ci_results[[method]])) ci_results[[method]]$runtime
              else NA_real_,
            observed_zero_prop = observed_zero_prop,
            structural_zero_prop = structural_zero_prop,
            sampling_zero_prop = sampling_zero_prop,
            generation_runtime = generation_runtime
          )
        }
      )
    }

    ## This binder cannot fail merely because method rows have
    ## different columns.
    results[[r]] <- rbind_fill(out_rep)
  }

  ## ------------------------------------------------------------
  ## Combine replicates safely
  ## ------------------------------------------------------------

  results_df <- rbind_fill(results)

  if(is.null(results_df)){
    stop("No simulation results were produced.")
  }

  results_df$DGP <- "ZIDM"
  results_df$n <- n
  results_df$m <- m
  results_df$k <- k
  results_df$zi_level <- zi_level
  results_df$depth_level <- depth_level

  results_df$ZIDM_paper_iterations <-
    if(include_ZIDM_paper) ZIDM_paper_iterations else NA_integer_

  results_df$ZIDM_paper_thin <-
    if(include_ZIDM_paper) ZIDM_paper_thin else NA_integer_

  results_df$ZIDM_paper_burnin <-
    if(include_ZIDM_paper) ZIDM_paper_burnin else NA_real_

  results_df$compute_CI <- compute_CI
  results_df$CI_level <- if(compute_CI) CI_level else NA_real_

  id_cols <- c(
    "DGP", "n", "m", "k",
    "zi_level", "depth_level",
    "rep", "method"
  )

  results_df <- results_df[
    ,
    c(id_cols, setdiff(names(results_df), id_cols)),
    drop = FALSE
  ]

  cat("\n========================================\n")
  cat("Finished ZIDM simulation\n")
  cat("Methods:", paste(methods_all, collapse = ", "), "\n")
  cat("Expected rows:", rep_k * length(methods_all), "\n")
  cat("Actual rows:", nrow(results_df), "\n")
  cat("Successful fits:", sum(results_df$success %in% TRUE, na.rm = TRUE), "\n")
  cat("Failed fits:", sum(results_df$success %in% FALSE, na.rm = TRUE), "\n")
  cat("========================================\n")

  results_df
}
})

# Original ZIDM scenario table and SLURM output handling.
run_scenario_ZIDM <- function() {








############################################################
## Puhti scenario control
## DGP = ZIDM
## FINAL calibrated design
## NO dataset filtering
############################################################

scenario_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(scenario_id)) scenario_id <- 1
cat("SLURM_ARRAY_TASK_ID =", scenario_id, "\n")

scenario_table <- data.frame(
  scenario = 1:6,
  n = c(50,50,50,100,100,100),
  m = c(100,100,100,50,50,50),
  zi_level = c("low","medium","high","low","medium","high"),
  depth_level = c("low","medium","high","low","medium","high"),
  zeta_low = c(0.05,0.20,0.35,0.05,0.20,0.35),
  zeta_high = c(0.15,0.30,0.45,0.15,0.30,0.45),
  depth_low = c(320,800,10000,140,320,4000),
  depth_high = c(480,1200,15000,210,480,6000),
  stringsAsFactors = FALSE
)

scenario_table$rep_k <- 100

scenario_table$label <- paste0(
  "ZIDM_n", scenario_table$n,
  "_m", scenario_table$m,
  "_ZI", scenario_table$zi_level,
  "_N", scenario_table$depth_low, "-", scenario_table$depth_high,
  "_R", scenario_table$rep_k
)

print(scenario_table)

if(scenario_id < 1 || scenario_id > nrow(scenario_table)){
  stop("scenario_id out of range.")
}

row <- scenario_table[scenario_id, , drop = FALSE]

cat("\n========================================\n")
cat("Running FINAL ZIDM scenario\n")
cat("========================================\n")
cat("scenario =", row$scenario, "\n")
cat("n =", row$n, "\n")
cat("m =", row$m, "\n")
cat("ZI =", row$zi_level, "\n")
cat("zeta_j ~ U(", row$zeta_low, ", ", row$zeta_high, ")\n", sep = "")
cat("N_i ~ discrete U{", row$depth_low, ", ", row$depth_high, "}\n", sep = "")
cat("replications =", row$rep_k, "\n")
cat("NO dataset filtering\n")
cat("label =", row$label, "\n")

output_dir <- "results_ZIDM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outfile_rds <- file.path(output_dir, paste0(row$label, ".rds"))
outfile_csv <- file.path(output_dir, paste0(row$label, "_full.csv"))

if(file.exists(outfile_rds)){
  cat("Skipping existing file:", outfile_rds, "\n")
  quit(save = "no")
}

res <- try(
  run_simulation_ZIDM(
    rep_k = row$rep_k,
    n = row$n,
    m = row$m,
    k = 2,
    zi_level = row$zi_level,
    depth_level = row$depth_level,
    init_method = "DS",
    seed_base = 2026000 + scenario_id * 10000
  ),
  silent = TRUE
)

if(inherits(res, "try-error")){
  cat("ZIDM scenario failed.\n")
  failed_file <- file.path(output_dir, paste0("FAILED_", row$label, ".rds"))
  saveRDS(list(scenario = as.list(row), error = as.character(res)), failed_file)
  stop(res)
}

## Metadata only; no replicate is screened or rejected.
res$design_zeta_low <- row$zeta_low
res$design_zeta_high <- row$zeta_high
res$design_depth_low <- row$depth_low
res$design_depth_high <- row$depth_high
res$design_depth_distribution <- "discrete_uniform"

saveRDS(res, outfile_rds)
write.csv(res, outfile_csv, row.names = FALSE)

cat("\nSaved:\n")
cat(outfile_rds, "\n")
cat(outfile_csv, "\n")
cat("Finished ZIDM scenario:", row$label, "\n")

}


# ==================== ZIMN generator ====================
run_simulation_ZIMN <- local({
  # Generator-specific implementation: generate_ZIDM_data
generate_ZIDM_data <- function(
    
  # ============================================================
  # Basic simulation dimensions
  # ============================================================
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.40)
  # high   : U(0.50, 0.70)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional user-specified true values
  # ============================================================
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter distributions from simulation table
  # ============================================================
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenario settings
  # If supplied, override zi_level/depth_level
  # ============================================================
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Structural-zero handling
  # ============================================================
  zero_inflation = TRUE,
  
  # "zero_row":
  # if all taxa structurally absent, Y_i = 0
  #
  # "resample":
  # redraw Z until at least one taxon survives
  all_zero_action = c("zero_row", "resample"),
  
  max_resample = 1000,
  
  # ============================================================
  # Output control
  # ============================================================
  return_Z = TRUE,
  return_lambda = TRUE,
  return_eta = TRUE,
  return_alpha_ij = TRUE,
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  all_zero_action <- match.arg(all_zero_action)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Scenario definitions
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.40),
      
      high = c(0.50, 0.70)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be a positive number.")
  }
  
  
  # ============================================================
  # 2. Generate observed covariates X
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 3. Generate latent variables
  #
  # u_i ~ N(0, I_k)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 4. Generate beta0
  #
  # beta_0j ~ U(-1,1)
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Generate covariate effects
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Generate factor loadings
  #
  # gamma_j ~ U(-2,2)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Zero-inflation probabilities
  #
  # zeta_j ~ Uniform(...)
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 8. Linear predictor
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 9. Dirichlet/Gamma shape
  #
  # alpha_ij = exp(eta_ij)
  # ============================================================
  
  alpha_ij <- exp(eta)
  
  
  if (
    any(!is.finite(alpha_ij)) ||
    any(alpha_ij <= 0)
  ) {
    stop(
      "Non-finite or non-positive exp(eta) encountered."
    )
  }
  
  
  # ============================================================
  # 10. Gamma latent variables
  #
  # lambda0_ij ~ Gamma(exp(eta_ij), rate = 1)
  # ============================================================
  
  lambda0 <- matrix(
    
    rgamma(
      n * m,
      shape = as.vector(alpha_ij),
      rate = 1
    ),
    
    nrow = n,
    ncol = m
  )
  
  
  # ============================================================
  # 11. Composition BEFORE zero inflation
  #
  # theta0_ij =
  # lambda0_ij / sum_l lambda0_il
  # ============================================================
  
  lambda0_rowsum <- rowSums(lambda0)
  
  rho0_true <- lambda0 /
    lambda0_rowsum
  
  
  # ============================================================
  # 12. Structural-zero process
  #
  # zeta_j = structural-zero probability
  #
  # Therefore:
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  #
  # Z = 1 : taxon survives
  # Z = 0 : structural zero
  # ============================================================
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      rep(
        1 - zeta,
        each = n
      ),
      nrow = n,
      ncol = m
    )
    
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
    
  } else {
    
    Z <- matrix(
      1L,
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 13. Optional handling of rows with all structural zeros
  # ============================================================
  
  if (
    zero_inflation &&
    all_zero_action == "resample"
  ) {
    
    zero_rows <- which(
      rowSums(Z) == 0
    )
    
    iter <- 0
    
    while (length(zero_rows) > 0) {
      
      iter <- iter + 1
      
      if (iter > max_resample) {
        stop(
          "Exceeded max_resample for all-zero structural rows."
        )
      }
      
      for (i in zero_rows) {
        
        Z[i, ] <- rbinom(
          m,
          size = 1,
          prob = 1 - zeta
        )
      }
      
      zero_rows <- which(
        rowSums(Z) == 0
      )
    }
  }
  
  
  # ============================================================
  # 14. Apply structural zeros
  # ============================================================
  
  lambdaZI <- lambda0 * Z
  
  row_mass <- rowSums(lambdaZI)
  
  
  # ============================================================
  # 15. ZI composition
  # ============================================================
  
  rhoZI_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  
  positive_rows <- row_mass > 0
  
  
  rhoZI_true[positive_rows, ] <-
    lambdaZI[positive_rows, , drop = FALSE] /
    row_mass[positive_rows]
  
  
  # ============================================================
  # 16. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.integer(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  
  N <- as.integer(N)
  
  
  # ============================================================
  # 17. Multinomial observation
  #
  # If all taxa are structurally zero:
  #
  # Y_i = 0
  #
  # otherwise:
  #
  # Y_i ~ Multinomial(N_i, rhoZI_i)
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  
  for (i in seq_len(n)) {
    
    if (row_mass[i] > 0) {
      
      Y[i, ] <- as.integer(
        
        rmultinom(
          n = 1,
          size = N[i],
          prob = rhoZI_true[i, ]
        )
      )
    }
  }
  
  
  # ============================================================
  # 18. Useful simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  all_structural_zero_prop <-
    mean(rowSums(Z) == 0)
  
  
  # ============================================================
  # 19. Return EVERYTHING needed for parameter recovery
  # ============================================================
  
  out <- list(
    
    # --------------------------------------------------------
    # observed data
    # --------------------------------------------------------
    
    Y = Y,
    X = X,
    N = N,
    
    
    # --------------------------------------------------------
    # TRUE composition
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for model recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # TRUE latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # TRUE zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # true linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    alpha_ij = alpha_ij,
    
    
    # --------------------------------------------------------
    # scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    all_structural_zero_prop =
      all_structural_zero_prop
  )
  
  
  # ============================================================
  # Optional latent outputs
  # ============================================================
  
  if (return_Z) {
    
    out$Z <- Z
  }
  
  
  if (return_lambda) {
    
    out$lambda0 <- lambda0
    
    out$lambdaZI <- lambdaZI
  }
  
  
  if (!return_eta) {
    
    out$eta <- NULL
  }
  
  
  if (!return_alpha_ij) {
    
    out$alpha_ij <- NULL
  }
  
  
  return(out)
}

  # Generator-specific implementation: generate_ZIMN_data
generate_ZIMN_data <- function(
    
  # ============================================================
  # Basic dimensions
  # ============================================================
  
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.30)
  # high   : U(0.35, 0.45)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional true parameter values
  #
  # IMPORTANT:
  # parameters are now defined directly in TAXON SPACE
  # ============================================================
  
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter ranges
  # ============================================================
  
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenarios
  # ============================================================
  
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Zero inflation
  # ============================================================
  
  zero_inflation = TRUE,
  
  # ============================================================
  # Random seed
  # ============================================================
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 1 || m != as.integer(m)) {
    stop("m must be an integer greater than 1.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Simulation scenario settings
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.30),
      
      high = c(0.35, 0.45)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be positive.")
  }
  
  
  # ============================================================
  # 2. Latent variables
  #
  # u_i ~ N(0, I_k)
  #
  # U: n x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 3. Observed covariates
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  #
  # X: n x 2
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 4. Species-specific intercepts
  #
  # beta_0j ~ U(-1,1)
  #
  # beta0: length m
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Covariate coefficients
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # beta: m x 2
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Latent loadings
  #
  # gamma_j ~ U(-2,2)
  #
  # Gamma: m x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Linear predictor DIRECTLY in taxon space
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  #
  # eta: n x m
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  #
  # X:        n x 2
  # t(beta):  2 x m
  #
  # result:   n x m
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  #
  # U:         n x k
  # t(Gamma):  k x m
  #
  # result:    n x m
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 8. TRUE baseline composition
  #
  # rho_ij =
  #
  # exp(eta_ij)
  # ------------------
  # sum_l exp(eta_il)
  #
  # stable softmax
  # ============================================================
  
  rho0_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    row_max <- max(
      eta[i, ]
    )
    
    tmp <- exp(
      eta[i, ] - row_max
    )
    
    rho0_true[i, ] <-
      tmp / sum(tmp)
  }
  
  
  # ============================================================
  # 9. Zero-inflation probabilities
  #
  # zeta_j = P(structural zero)
  #
  # zeta_j ~ Uniform(...)
  #
  # zeta: length m
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 10. Structural-zero indicators
  #
  # zeta_j = P(structural zero)
  #
  # Define:
  #
  # Z_ij = 1 : taxon is PRESENT
  # Z_ij = 0 : structural zero
  #
  # therefore
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  # ============================================================
  
  Z <- matrix(
    1L,
    nrow = n,
    ncol = m
  )
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      
      rep(
        1 - zeta,
        each = n
      ),
      
      nrow = n,
      ncol = m
    )
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 11. Apply structural zeros
  #
  # rhoZI_ij proportional to
  #
  # Z_ij * rho0_ij
  # ============================================================
  
  rhoZI_raw <-
    rho0_true * Z
  
  rs <- rowSums(
    rhoZI_raw
  )
  
  
  # ------------------------------------------------------------
  # Extremely rare case:
  # all taxa removed in one sample
  #
  # Activate one taxon so multinomial probabilities are valid
  # ------------------------------------------------------------
  
  bad <- which(
    rs == 0
  )
  
  if (length(bad) > 0) {
    
    for (i in bad) {
      
      j <- sample(
        seq_len(m),
        size = 1
      )
      
      Z[i, j] <- 1L
      
      rhoZI_raw[i, j] <-
        rho0_true[i, j]
    }
    
    rs <- rowSums(
      rhoZI_raw
    )
  }
  
  
  # ============================================================
  # 12. TRUE zero-inflated composition
  #
  # rhoZI_ij =
  #
  # Z_ij rho0_ij
  # --------------------
  # sum_l Z_il rho0_il
  # ============================================================
  
  rhoZI_true <-
    rhoZI_raw / rs
  
  
  # ============================================================
  # 13. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.numeric(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  N <- as.integer(
    round(N)
  )
  
  if (any(N < 0)) {
    stop("N must be nonnegative.")
  }
  
  
  # ============================================================
  # 14. Multinomial sampling
  #
  # Y_i ~ Multinomial(
  #     N_i,
  #     rhoZI_true_i
  # )
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    Y[i, ] <- as.integer(
      
      rmultinom(
        n = 1,
        size = N[i],
        prob = rhoZI_true[i, ]
      )
    )
  }
  
  
  # ============================================================
  # 15. Simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  
  # ============================================================
  # 16. Return EVERYTHING needed for simulation evaluation
  # ============================================================
  
  return(list(
    
    # --------------------------------------------------------
    # Observed data
    # --------------------------------------------------------
    
    Y = Y,
    
    X = X,
    
    N = N,
    
    
    # --------------------------------------------------------
    # True compositions
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for composition recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    P_true = rho0_true,
    
    P_ZI = rhoZI_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    #
    # NOW ALL PARAMETERS ARE IN TAXON SPACE
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # True latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    Beta = t(beta),
    
    Lambda = Gamma,
    
    
    # --------------------------------------------------------
    # True zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # Linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    eta_taxon = eta,
    
    
    # --------------------------------------------------------
    # Scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    
    # --------------------------------------------------------
    # Structural-zero indicator
    # --------------------------------------------------------
    
    Z = Z
  ))
}

  # Generator-specific implementation: fit_gllvm_ZINB
fit_gllvm_ZINB <- function(
    Y,
    X = NULL,
    k = 2,
    seed = 123
){
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  m <- ncol(Y)
  
  if(!is.null(X)){
    
    X <- as.matrix(X)
    
    if(nrow(X) != n){
      stop("X must have the same number of rows as Y.")
    }
  }
  
  
  ##########################################################
  ## 1. Fit ZINB-GLLVM
  ##
  ## Improved convergence:
  ## - more iterations
  ## - multiple initial values
  ## - residual-based starting values
  ##########################################################
  
  fit <- gllvm::gllvm(
    y = Y,
    X = X,
    family = "ZINB",
    num.lv = k,
    row.eff = "fixed",
    method = "VA",
    sd.errors = FALSE,
    seed = seed,
    
    control = list(
      maxit = 20000,
      max.iter = 20000,
      reltol = 1e-10
    ),
    
    control.start = list(
      starting.val = "res",
      n.init = 3
    )
  )
  
  
  ##########################################################
  ## 2. Extract fitted model parameters
  ##########################################################
  
  # species-specific intercepts beta_0j
  beta0_hat <- fit$params$beta0
  
  # environmental coefficients beta_j
  beta_hat <- fit$params$Xcoef
  
  # latent scores u_i
  u_hat <- as.matrix(
    gllvm::getLV(fit)
  )
  
  # properly scaled loading vectors gamma_j
  gamma_hat <- as.matrix(
    gllvm::getLoadings(fit)
  )
  
  # fixed sample / row effects alpha_i
  # sample 1 is the reference, alpha_1 = 0
  alpha_hat <- c(
    0,
    fit$params$row.params.fixed
  )
  
  
  ##########################################################
  ## 3. Zero-inflation and dispersion parameters
  ##########################################################
  
  # ZI probability p_j
  # This corresponds to zeta_j in our notation
  zeta_hat <- fit$params$phi
  
  # NB dispersion parameter
  dispersion_hat <- fit$params$ZINB.phi
  
  # inverse NB dispersion
  inv_dispersion_hat <- fit$params$ZINB.inv.phi
  
  
  ##########################################################
  ## 4. Variational covariance matrices A_i
  ##
  ## A_hat[i,,] = A_i
  ## dimension n x k x k
  ##########################################################
  
  A_hat <- fit$A
  
  
  ##########################################################
  ## 5. Construct fitted linear predictor eta
  ##
  ## eta_ij =
  ## alpha_i
  ## + beta_0j
  ## + x_i^T beta_j
  ## + u_i^T gamma_j
  ##########################################################
  
  eta_hat <- matrix(
    beta0_hat,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  ##########################################################
  ## add row effect alpha_i
  ##########################################################
  
  eta_hat <- eta_hat +
    matrix(
      alpha_hat,
      nrow = n,
      ncol = m,
      byrow = FALSE
    )
  
  
  ##########################################################
  ## add environmental contribution X beta
  ##########################################################
  
  if(!is.null(X) && !is.null(beta_hat)){
    
    beta_mat <- as.matrix(beta_hat)
    
    # usual gllvm orientation: m x q
    if(
      nrow(beta_mat) == m &&
      ncol(beta_mat) == ncol(X)
    ){
      
      eta_hat <- eta_hat +
        X %*% t(beta_mat)
      
    } else if(
      nrow(beta_mat) == ncol(X) &&
      ncol(beta_mat) == m
    ){
      
      eta_hat <- eta_hat +
        X %*% beta_mat
      
    } else {
      
      stop(
        "Dimension of fit$params$Xcoef is not compatible with X."
      )
    }
  }
  
  
  ##########################################################
  ## add latent contribution U Gamma
  ##########################################################
  
  # convert loadings to k x m
  if(
    nrow(gamma_hat) == m &&
    ncol(gamma_hat) == k
  ){
    
    Gamma_mat <- t(gamma_hat)
    
  } else if(
    nrow(gamma_hat) == k &&
    ncol(gamma_hat) == m
  ){
    
    Gamma_mat <- gamma_hat
    
  } else {
    
    stop(
      "Dimension of getLoadings(fit) is not compatible with k and m."
    )
  }
  
  eta_hat <- eta_hat +
    u_hat %*% Gamma_mat
  
  
  ##########################################################
  ## 6. Variational quadratic correction
  ##
  ## cQ_ij =
  ## 0.5 * gamma_j^T A_i gamma_j
  ##########################################################
  
  cQ_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    Ai <- A_hat[i, , ]
    
    for(j in seq_len(m)){
      
      gamma_j <- Gamma_mat[, j]
      
      cQ_hat[i, j] <-
        0.5 *
        as.numeric(
          t(gamma_j) %*%
            Ai %*%
            gamma_j
        )
    }
  }
  
  
  ##########################################################
  ## 7. Variationally corrected predictor
  ##########################################################
  
  eta_corrected_hat <-
    eta_hat + cQ_hat
  
  
  ##########################################################
  ## 8. Corrected fitted composition rho
  ##
  ## rho_ij =
  ## exp(eta_ij + cQ_ij)
  ## -------------------------
  ## sum_l exp(eta_il + cQ_il)
  ##########################################################
  
  rho_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    row_max <- max(
      eta_corrected_hat[i, ]
    )
    
    tmp <- exp(
      eta_corrected_hat[i, ] -
        row_max
    )
    
    rho_hat[i, ] <-
      tmp / sum(tmp)
  }
  
  
  ##########################################################
  ## 9. Return everything needed for simulation
  ##########################################################
  
  return(list(
    
    ########################################################
    ## regression / fitted model parameters
    ########################################################
    
    alpha_hat = alpha_hat,
    
    beta0_hat = beta0_hat,
    
    beta_hat = beta_hat,
    
    
    ########################################################
    ## latent variables
    ########################################################
    
    u_hat = u_hat,
    
    gamma_hat = gamma_hat,
    
    
    ########################################################
    ## zero-inflation parameters
    ########################################################
    
    zeta_hat = zeta_hat,
    
    dispersion_hat = dispersion_hat,
    
    inv_dispersion_hat = inv_dispersion_hat,
    
    
    ########################################################
    ## variational quantities
    ########################################################
    
    A_hat = A_hat,
    
    cQ_hat = cQ_hat,
    
    
    ########################################################
    ## predictors
    ########################################################
    
    eta_hat = eta_hat,
    
    eta_corrected_hat = eta_corrected_hat,
    
    
    ########################################################
    ## fitted composition
    ########################################################
    
    rho_hat = rho_hat,
    
    
    ########################################################
    ## raw loading parameters
    ########################################################
    
    theta_raw = fit$params$theta,
    
    sigma_lv = fit$params$sigma.lv,
    
    
    ########################################################
    ## convenient aliases
    ########################################################
    
    ordination = u_hat,
    
    loadings = gamma_hat,
    
    
    ########################################################
    ## complete fitted object
    ########################################################
    
    model = fit
  ))
}

run_simulation_ZIMN <- function(
    rep_k = 500, #100
    n = 50,
    m = 100,
    k = 2,
    zi_level = c("low", "medium", "high"),
    depth_level = c("low", "medium", "high"),
    init_method = "DS",
    seed_base = 2026,
    compute_CI = TRUE,
    CI_level = 0.95,
    include_ZIDM_paper = TRUE,
    ZIDM_paper_iterations = 20000, #10000
    ZIDM_paper_thin = 10,
    ZIDM_paper_burnin = 0.5
){
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)

  ## FINAL calibrated sequencing-depth distribution.
  ## No dataset filtering is applied.
  if(n == 50 && m == 100){
    calibrated_depth_range <- switch(
      zi_level,
      low = c(200L, 300L),
      medium = c(320L, 480L),
      high = c(800L, 1200L)
    )
  } else if(n == 100 && m == 50){
    calibrated_depth_range <- switch(
      zi_level,
      low = c(100L, 150L),
      medium = c(140L, 210L),
      high = c(400L, 600L)
    )
  } else {
    stop("No calibrated ZIMN depth range for this (n,m) setting.")
  }
  
  
  ############################################################
  ## TMB models: NO STR
  ############################################################
  
  tmb_model_list <- list(
    
    list(
      name = "ZIMN_soft",
      model = "ZIMN",
      extra = c(0,0)
    ),
    
    list(
      name = "ZIMN_hard",
      model = "ZIMN",
      extra = c(0,1)
    ),
    
    list(
      name = "ZIDM_soft",
      model = "ZIDM",
      extra = c(0,0)
    ),
    
    list(
      name = "ZIDM_hard",
      model = "ZIDM",
      extra = c(0,1)
    )
  )
  
  
  methods_all <- c(
    "ZINB_gllvm",
    "ZIMN_soft",
    "ZIMN_hard",
    "ZIDM_soft",
    "ZIDM_hard"
  )
  
  if(include_ZIDM_paper){
    methods_all <- c(
      methods_all,
      "ZIDM_paper"
    )
  }
  
  methods_all <- c(
    methods_all,
    "ZR"
  )
  
  
  ############################################################
  ## SAFE
  ############################################################
  
  safe_fit <- function(expr){
    
    out <- try(
      expr,
      silent = TRUE
    )
    
    if(inherits(out, "try-error")){
      return(NULL)
    }
    
    out
  }
  
  
  safe_eval <- function(expr){
    
    out <- try(
      expr,
      silent = TRUE
    )
    
    if(inherits(out, "try-error")){
      return(NULL)
    }
    
    out
  }
  
  
  ############################################################
  ## convergence checks
  ############################################################
  
  valid_TMB <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    
    if(is.null(fit$convergence)) return(FALSE)
    
    if(fit$convergence != 0) return(FALSE)
    
    if(is.null(fit$rho_hat)) return(FALSE)
    
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    TRUE
  }
  
  
  valid_gllvm <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    
    if(is.null(fit$model)) return(FALSE)
    
    if(is.null(fit$model$convergence)) return(FALSE)
    
    if(!isTRUE(fit$model$convergence)) return(FALSE)
    
    if(is.null(fit$rho_hat)) return(FALSE)
    
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    TRUE
  }
  
  
  ############################################################
  ## paper ZIDM validity
  ############################################################
  
  valid_ZIDM_paper <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    if(is.null(fit$rho_hat)) return(FALSE)
    if(!is.matrix(fit$rho_hat)) return(FALSE)
    if(!all(dim(fit$rho_hat) == c(n, m))) return(FALSE)
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    if(is.null(fit$beta0_hat)) return(FALSE)
    if(length(fit$beta0_hat) != m) return(FALSE)
    if(any(!is.finite(fit$beta0_hat))) return(FALSE)
    
    if(is.null(fit$beta_hat)) return(FALSE)
    if(!all(dim(fit$beta_hat) == c(m, 2))) return(FALSE)
    if(any(!is.finite(fit$beta_hat))) return(FALSE)
    
    if(is.null(fit$zeta_hat)) return(FALSE)
    if(length(fit$zeta_hat) != m) return(FALSE)
    if(any(!is.finite(fit$zeta_hat))) return(FALSE)
    
    TRUE
  }
  
  
  ############################################################
  ## NA result
  ############################################################
  
    ############################################################
  ## CI summary helper
  ############################################################

  empty_ci_summary <- function(requested = compute_CI){
    list(
      ci_success = if(requested) FALSE else NA,
      beta0_Coverage = NA_real_,
      beta1_Coverage = NA_real_,
      beta2_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta1_meanSE = NA_real_,
      beta2_meanSE = NA_real_,
      beta0_meanWidth = NA_real_,
      beta1_meanWidth = NA_real_,
      beta2_meanWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta1_nValid = NA_integer_,
      beta2_nValid = NA_integer_,
      beta0_validProp = NA_real_,
      beta1_validProp = NA_real_,
      beta2_validProp = NA_real_
    )
  }


  ############################################################
  ## NA row
  ############################################################
  
  make_na_row <- function(
    r,
    method,
    runtime = NA_real_,
    convergence = NA,
    ci_runtime = NA_real_,
    observed_zero_prop = NA_real_,
    structural_zero_prop = NA_real_,
    sampling_zero_prop = NA_real_,
    generation_runtime = NA_real_
  ){

    runtime_total <- if(
      is.finite(runtime) && is.finite(ci_runtime)
    ){
      runtime + ci_runtime
    } else if(is.finite(runtime)){
      runtime
    } else {
      NA_real_
    }

    data.frame(
      rep = r,
      method = method,
      success = FALSE,
      convergence = convergence,

      runtime = runtime,
      ci_runtime = ci_runtime,
      runtime_total = runtime_total,
      ci_success = if(compute_CI && method != "ZR") FALSE else NA,

      Frobenius = NA_real_,
      KL = NA_real_,
      Shannon_MSE = NA_real_,
      Simpson_MSE = NA_real_,
      Procrustes = NA_real_,

      beta0_Bias = NA_real_,
      beta0_AbsBias = NA_real_,
      beta0_RMSE = NA_real_,

      beta1_Bias = NA_real_,
      beta1_AbsBias = NA_real_,
      beta1_RMSE = NA_real_,

      beta2_Bias = NA_real_,
      beta2_AbsBias = NA_real_,
      beta2_RMSE = NA_real_,

      zeta_Bias = NA_real_,
      zeta_AbsBias = NA_real_,
      zeta_RMSE = NA_real_,

      beta0_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta0_meanCIWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta0_validProp = NA_real_,

      beta1_Coverage = NA_real_,
      beta1_meanSE = NA_real_,
      beta1_meanCIWidth = NA_real_,
      beta1_nValid = NA_integer_,
      beta1_validProp = NA_real_,

      beta2_Coverage = NA_real_,
      beta2_meanSE = NA_real_,
      beta2_meanCIWidth = NA_real_,
      beta2_nValid = NA_integer_,
      beta2_validProp = NA_real_,

      observed_zero_prop = observed_zero_prop,
      structural_zero_prop = structural_zero_prop,
      sampling_zero_prop = sampling_zero_prop,
      generation_runtime = generation_runtime,

      stringsAsFactors = FALSE
    )
  }

  ############################################################
  ## safe bind by column name
  ## Different method rows are padded with NA before rbind.
  ############################################################

  rbind_fill <- function(x){

    x <- x[
      !vapply(
        x,
        is.null,
        logical(1)
      )
    ]

    if(length(x) == 0L){
      return(NULL)
    }

    all_cols <- unique(
      unlist(
        lapply(
          x,
          names
        ),
        use.names = FALSE
      )
    )

    x <- lapply(
      x,
      function(z){

        miss <- setdiff(
          all_cols,
          names(z)
        )

        if(length(miss) > 0L){
          for(nm in miss){
            z[[nm]] <- NA
          }
        }

        z[
          ,
          all_cols,
          drop = FALSE
        ]
      }
    )

    out <- do.call(
      rbind,
      x
    )

    rownames(out) <- NULL

    out
  }
  
  
  ############################################################
  ## storage
  ############################################################
  
  results <- vector(
    "list",
    rep_k
  )
  
  
  ############################################################
  ## repetitions
  ############################################################
  
  for(r in seq_len(rep_k)){
    
    cat(
      "\n========================================\n"
    )
    
    cat(
      "ZIMN replicate:",
      r,
      "/",
      rep_k,
      "\n"
    )
    
    cat(
      "n =", n,
      "| m =", m,
      "| zi =", zi_level,
      "| depth =", depth_level,
      "\n"
    )
    
    cat(
      "========================================\n"
    )
    
    
    rep_seed <- seed_base + r
    
    
    ##########################################################
    ## FULL REP PROTECTION
    ##########################################################
    
    rep_result <- try({
      
      ########################################################
      ## 1. Generate ZIMN data
      ########################################################
      
      t0_generate <- Sys.time()
      
      ## N_i ~ discrete Uniform(L,U), independently for each sample.
      set.seed(rep_seed + 9000000L)
      N_now <- sample(
        seq.int(calibrated_depth_range[1], calibrated_depth_range[2]),
        size = n,
        replace = TRUE
      )

      dat <- generate_ZIMN_data(
        
        n = n,
        
        m = m,
        
        k = k,
        
        zi_level = zi_level,
        
        depth_level = depth_level,
        
        N = N_now,
        
        zero_inflation = TRUE,
        
        seed = rep_seed
      )
      
      
      generation_runtime <- as.numeric(
        Sys.time() - t0_generate,
        units = "secs"
      )
      
      
      ########################################################
      ## true quantities
      ########################################################
      
      Y <- dat$Y
      
      X <- dat$X
      
      rho_true <- dat$rho_true
      
      U_true <- dat$U
      
      beta0_true <- dat$beta0
      
      beta1_true <- dat$beta1
      
      beta2_true <- dat$beta2
      
      zeta_true <- dat$zeta
      
      
      ########################################################
      ## fitted models
      ########################################################
      
      fit_results <- list()
      
      
      ########################################################
      ## 2. gllvm ZINB
      ########################################################
      
      cat("  fitting ZINB_gllvm ...\n")
      
      t0 <- Sys.time()
      
      fit_g <- safe_fit(
        
        fit_gllvm_ZINB(
          
          Y = Y,
          
          X = X,
          
          k = k,
          
          seed = rep_seed
        )
      )
      
      
      runtime_g <- as.numeric(
        Sys.time() - t0,
        units = "secs"
      )
      
      
      if(valid_gllvm(fit_g)){
        
        fit_results$ZINB_gllvm <- list(
          
          fit = fit_g,
          
          runtime = runtime_g,
          
          success = TRUE,
          
          convergence = TRUE,
          
          type = "gllvm"
        )
        
      } else {
        
        fit_results$ZINB_gllvm <- list(
          
          fit = NULL,
          
          runtime = runtime_g,
          
          success = FALSE,
          
          convergence = FALSE,
          
          type = "gllvm"
        )
      }
      
      
      ########################################################
      ## 3. Four TMB models
      ########################################################
      
      for(mod in tmb_model_list){
        
        cat(
          "  fitting",
          mod$name,
          "...\n"
        )
        
        
        t0 <- Sys.time()
        
        
        fit_tmp <- safe_fit(
          
          fit_TMB(
            
            Y = Y,
            
            X = X,
            
            k = k,
            
            model = mod$model,
            
            extra = mod$extra,
            
            init_method = init_method
          )
        )
        
        
        runtime_tmp <- as.numeric(
          Sys.time() - t0,
          units = "secs"
        )
        
        
        if(valid_TMB(fit_tmp)){
          
          fit_results[[mod$name]] <- list(
            
            fit = fit_tmp,
            
            runtime = runtime_tmp,
            
            success = TRUE,
            
            convergence = fit_tmp$convergence,
            
            type = "TMB"
          )
          
        } else {
          
          fit_results[[mod$name]] <- list(
            
            fit = NULL,
            
            runtime = runtime_tmp,
            
            success = FALSE,
            
            convergence =
              if(!is.null(fit_tmp$convergence)){
                fit_tmp$convergence
              } else {
                NA
              },
            
            type = "TMB"
          )
        }
      }
      
      
      ########################################################
      ## 4. Koslovsky paper ZIDM
      ## standard regression, no BVS, no CI
      ########################################################
      
      if(include_ZIDM_paper){
        
        cat("  fitting ZIDM_paper ...\n")
        
        t0 <- Sys.time()
        
        fit_paper <- safe_fit(
          fit_ZIDM_paper(
            Y = Y,
            X = X,
            X_theta = NULL,
            iterations = ZIDM_paper_iterations,
            thin = ZIDM_paper_thin,
            burnin_prop = ZIDM_paper_burnin,
            seed = rep_seed
          )
        )
        
        runtime_paper <- as.numeric(
          Sys.time() - t0,
          units = "secs"
        )
        
        if(valid_ZIDM_paper(fit_paper)){
          
          fit_results$ZIDM_paper <- list(
            fit = fit_paper,
            runtime = runtime_paper,
            success = TRUE,
            convergence = TRUE,
            type = "ZIDM_paper"
          )
          
        } else {
          
          fit_results$ZIDM_paper <- list(
            fit = NULL,
            runtime = runtime_paper,
            success = FALSE,
            convergence = FALSE,
            type = "ZIDM_paper"
          )
        }
      }
      
      
      ########################################################
      ## 5. ZR baseline
      ########################################################
      
      cat("  fitting ZR ...\n")
      
      t0 <- Sys.time()
      
      fit_zr <- safe_fit(
        fit_ZR(Y)
      )
      
      
      runtime_zr <- as.numeric(
        Sys.time() - t0,
        units = "secs"
      )
      
      
      if(
        !is.null(fit_zr) &&
        !is.null(fit_zr$rho_hat) &&
        all(is.finite(fit_zr$rho_hat))
      ){
        
        fit_results$ZR <- list(
          
          fit = fit_zr,
          
          runtime = runtime_zr,
          
          success = TRUE,
          
          convergence = TRUE,
          
          type = "ZR"
        )
        
      } else {
        
        fit_results$ZR <- list(
          
          fit = NULL,
          
          runtime = runtime_zr,
          
          success = FALSE,
          
          convergence = FALSE,
          
          type = "ZR"
        )
      }
      
      
      ########################################################
      ## 6. Evaluation
      ########################################################
      
      ########################################################
      ## Confidence intervals / credible intervals
      ########################################################

      ci_results <- list()

      for(method in methods_all){

        ci_results[[method]] <- list(
          ci = NULL,
          summary = empty_ci_summary(compute_CI),
          runtime = NA_real_
        )

        if(!compute_CI){
          next
        }

        if(method == "ZR"){
          ci_results[[method]] <- list(
            ci = NULL,
            summary = empty_ci_summary(FALSE),
            runtime = NA_real_
          )
          next
        }

        fit_obj <- fit_results[[method]]

        if(is.null(fit_obj) || !isTRUE(fit_obj$success) || is.null(fit_obj$fit)){
          next
        }

        t0_ci <- Sys.time()
        ci_tmp <- NULL

        if(fit_obj$type == "gllvm"){

          ci_tmp <- safe_eval(
            compute_gllvm_ZINB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )

        } else if(fit_obj$type == "TMB"){

          ci_tmp <- safe_eval(
            compute_TMB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )

        } else if(fit_obj$type == "ZIDM_paper"){

          ci_tmp <- safe_eval(
            compute_ZIDM_paper_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )
        }

        runtime_ci <- as.numeric(Sys.time() - t0_ci, units = "secs")

        if(is.null(ci_tmp)){
          ci_sum <- empty_ci_summary(TRUE)
        } else {
          ci_sum <- safe_eval(summarize_beta_CI(ci_tmp))
          if(is.null(ci_sum)){
            ci_sum <- empty_ci_summary(TRUE)
          }
        }

        ci_results[[method]] <- list(
          ci = ci_tmp,
          summary = ci_sum,
          runtime = runtime_ci
        )
      }


      out_rep <- list()
      
      
      for(method in methods_all){
        
        obj <- fit_results[[method]]
        
        
        ######################################################
        ## failed model
        ######################################################
        
        if(
          is.null(obj) ||
          !isTRUE(obj$success) ||
          is.null(obj$fit)
        ){
          
          out_rep[[method]] <- make_na_row(
            
            r = r,
            
            method = method,
            
            runtime =
              if(!is.null(obj$runtime)){
                obj$runtime
              } else {
                NA_real_
              }
          )
          
          next
        }
        
        
        method_eval <- try({

        fit <- obj$fit

        ci_obj <- ci_results[[method]]
        ci_sum <- ci_obj$summary
        ci_runtime <- ci_obj$runtime

        runtime_total <- if(is.finite(ci_runtime)){
          obj$runtime + ci_runtime
        } else {
          obj$runtime
        }
        
        
        ######################################################
        ## composition
        ######################################################
        
        comp <- safe_eval(
          
          evaluate_composition(
            
            rho_hat = fit$rho_hat,
            
            rho_true = rho_true
          )
        )
        
        
        if(is.null(comp)){
          
          comp <- list(
            Frobenius = NA,
            KL = NA,
            Shannon_MSE = NA,
            Simpson_MSE = NA
          )
        }
        
        
        ######################################################
        ## ordination
        ######################################################
        
        proc <- NA_real_
        
        
        if(method != "ZR" && method != "ZIDM_paper"){
          
          proc_tmp <- safe_eval(
            
            compute_procrustes(
              
              true_ord = U_true,
              
              est_ord = fit$ordination
            )
          )
          
          
          if(
            !is.null(proc_tmp) &&
            length(proc_tmp) == 1 &&
            is.finite(proc_tmp)
          ){
            
            proc <- proc_tmp
          }
        }
        
        
        ######################################################
        ## parameter defaults
        ######################################################
        
        empty_met <- list(
          Bias = NA,
          AbsBias = NA,
          RMSE = NA
        )
        
        
        beta0_met <- empty_met
        beta1_met <- empty_met
        beta2_met <- empty_met
        zeta_met <- empty_met
        
        
        ######################################################
        ## parameter evaluation
        ######################################################
        
        if(method != "ZR"){
          
          ####################################################
          ## beta0
          ####################################################
          
          tmp <- safe_eval(
            
            evaluate_parameter(
              
              est = fit$beta0_hat,
              
              true = beta0_true
            )
          )
          
          if(!is.null(tmp)){
            beta0_met <- tmp
          }
          
          
          ####################################################
          ## beta1 + beta2
          ####################################################
          
          if(obj$type %in% c("gllvm", "ZIDM_paper")){
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[,1],
                
                true = beta1_true
              )
            )
            
            if(!is.null(tmp)){
              beta1_met <- tmp
            }
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[,2],
                
                true = beta2_true
              )
            )
            
            if(!is.null(tmp)){
              beta2_met <- tmp
            }
            
            
          } else if(obj$type == "TMB"){
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[1,],
                
                true = beta1_true
              )
            )
            
            if(!is.null(tmp)){
              beta1_met <- tmp
            }
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[2,],
                
                true = beta2_true
              )
            )
            
            if(!is.null(tmp)){
              beta2_met <- tmp
            }
          }
          
          
          ####################################################
          ## zeta
          ####################################################
          
          tmp <- safe_eval(
            
            evaluate_parameter(
              
              est = fit$zeta_hat,
              
              true = zeta_true
            )
          )
          
          if(!is.null(tmp)){
            zeta_met <- tmp
          }
        }
        
        
        ######################################################
        ## result row
        ######################################################
        
        method_row <- data.frame(
          
          rep = r,
          
          method = method,
          
          success = TRUE,
          
          convergence = obj$convergence,
          
          runtime = obj$runtime,
          
          ci_runtime = ci_runtime,
          
          runtime_total = runtime_total,
          
          ci_success = ci_sum$ci_success,
          
          
          Frobenius =
            comp$Frobenius,
          
          KL =
            comp$KL,
          
          Shannon_MSE =
            comp$Shannon_MSE,
          
          Simpson_MSE =
            comp$Simpson_MSE,
          
          
          Procrustes =
            proc,
          
          
          beta0_Bias =
            beta0_met$Bias,
          
          beta0_AbsBias =
            beta0_met$AbsBias,
          
          beta0_RMSE =
            beta0_met$RMSE,
          
          
          beta1_Bias =
            beta1_met$Bias,
          
          beta1_AbsBias =
            beta1_met$AbsBias,
          
          beta1_RMSE =
            beta1_met$RMSE,
          
          
          beta2_Bias =
            beta2_met$Bias,
          
          beta2_AbsBias =
            beta2_met$AbsBias,
          
          beta2_RMSE =
            beta2_met$RMSE,
          
          
          zeta_Bias =
            zeta_met$Bias,
          
          zeta_AbsBias =
            zeta_met$AbsBias,
          
          zeta_RMSE =
            zeta_met$RMSE,
          
          
          beta0_Coverage = ci_sum$beta0_Coverage,
          
          
          beta0_meanSE = ci_sum$beta0_meanSE,
          
          
          beta0_meanCIWidth = ci_sum$beta0_meanWidth,
          
          
          beta0_nValid = ci_sum$beta0_nValid,
          
          
          beta0_validProp = ci_sum$beta0_validProp,
          
          
          beta1_Coverage = ci_sum$beta1_Coverage,
          
          
          beta1_meanSE = ci_sum$beta1_meanSE,
          
          
          beta1_meanCIWidth = ci_sum$beta1_meanWidth,
          
          
          beta1_nValid = ci_sum$beta1_nValid,
          
          
          beta1_validProp = ci_sum$beta1_validProp,
          
          
          beta2_Coverage = ci_sum$beta2_Coverage,
          
          
          beta2_meanSE = ci_sum$beta2_meanSE,
          
          
          beta2_meanCIWidth = ci_sum$beta2_meanWidth,
          
          
          beta2_nValid = ci_sum$beta2_nValid,
          
          
          beta2_validProp = ci_sum$beta2_validProp,

          
          
          observed_zero_prop =
            dat$observed_zero_prop,
          
          structural_zero_prop =
            dat$structural_zero_prop,
          
          sampling_zero_prop =
            dat$sampling_zero_prop,
          
          generation_runtime =
            generation_runtime,
          
          
          stringsAsFactors = FALSE
        )

        method_row

        }, silent = TRUE)

        if(inherits(method_eval, "try-error")){
          cat("  method", method, "evaluation FAILED -> saving NA for this method only\n")
          cat("    ", as.character(method_eval), "\n")
          out_rep[[method]] <- make_na_row(
            r = r, method = method,
            runtime = if(!is.null(obj$runtime)) obj$runtime else NA_real_,
            convergence = if(!is.null(obj$convergence)) obj$convergence else NA,
            observed_zero_prop = dat$observed_zero_prop,
            structural_zero_prop = dat$structural_zero_prop,
            sampling_zero_prop = dat$sampling_zero_prop,
            generation_runtime = generation_runtime
          )
        } else {
          out_rep[[method]] <- method_eval
        }
      }
      
      
      rbind_fill(out_rep)
      
    }, silent = TRUE)
    
    
    ##########################################################
    ## whole replicate failed
    ##########################################################
    
    if(inherits(rep_result, "try-error")){
      
      cat(
        "  replicate",
        r,
        "FAILED completely\n"
      )
      
      
      rep_result <- rbind_fill(
        lapply(
                  
                  methods_all,
                  
                  function(mm){
                    
                    make_na_row(
                      
                      r = r,
                      
                      method = mm
                    )
                  }
                )
      )
    }
    
    
    results[[r]] <- rep_result
  }
  
  
  ############################################################
  ## combine
  ############################################################
  
  results_df <- rbind_fill(results)
  
  
  rownames(results_df) <- NULL
  
  
  ############################################################
  ## scenario information
  ############################################################
  
  results_df$DGP <- "ZIMN"
  
  results_df$n <- n
  results_df$m <- m
  results_df$k <- k
  
  results_df$zi_level <- zi_level
  
  results_df$depth_level <- depth_level
  
  results_df$ZIDM_paper_iterations <- if(include_ZIDM_paper) ZIDM_paper_iterations else NA_integer_
  results_df$ZIDM_paper_thin <- if(include_ZIDM_paper) ZIDM_paper_thin else NA_integer_
  results_df$ZIDM_paper_burnin <- if(include_ZIDM_paper) ZIDM_paper_burnin else NA_real_
  
  
  results_df$compute_CI <- compute_CI
  results_df$CI_level <- if(compute_CI) CI_level else NA_real_

  id_cols <- c(
    "DGP",
    "n",
    "m",
    "k",
    "zi_level",
    "depth_level",
    "rep",
    "method"
  )
  
  
  results_df <- results_df[
    ,
    c(
      id_cols,
      setdiff(
        names(results_df),
        id_cols
      )
    )
  ]
  
  
  return(results_df)
}
})

# Original ZIMN scenario table and SLURM output handling.
run_scenario_ZIMN <- function() {




############################################################
## Puhti scenario control
## DGP = ZIMN
## FINAL calibrated design
## NO dataset filtering
############################################################

scenario_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(scenario_id)) scenario_id <- 1
cat("SLURM_ARRAY_TASK_ID =", scenario_id, "\n")

scenario_table <- data.frame(
  scenario = 1:6,
  n = c(50,50,50,100,100,100),
  m = c(100,100,100,50,50,50),
  zi_level = c("low","medium","high","low","medium","high"),
  depth_level = c("low","medium","high","low","medium","high"),
  zeta_low = c(0.05,0.20,0.35,0.05,0.20,0.35),
  zeta_high = c(0.15,0.30,0.45,0.15,0.30,0.45),
  depth_low = c(200,320,800,100,140,400),
  depth_high = c(300,480,1200,150,210,600),
  stringsAsFactors = FALSE
)

scenario_table$rep_k <- 100

scenario_table$label <- paste0(
  "ZIMN_n", scenario_table$n,
  "_m", scenario_table$m,
  "_ZI", scenario_table$zi_level,
  "_N", scenario_table$depth_low, "-", scenario_table$depth_high,
  "_R", scenario_table$rep_k
)

print(scenario_table)

if(scenario_id < 1 || scenario_id > nrow(scenario_table)){
  stop("scenario_id out of range.")
}

row <- scenario_table[scenario_id, , drop = FALSE]

cat("\n========================================\n")
cat("Running FINAL ZIMN scenario\n")
cat("========================================\n")
cat("scenario =", row$scenario, "\n")
cat("n =", row$n, "\n")
cat("m =", row$m, "\n")
cat("ZI =", row$zi_level, "\n")
cat("zeta_j ~ U(", row$zeta_low, ", ", row$zeta_high, ")\n", sep = "")
cat("N_i ~ discrete U{", row$depth_low, ", ", row$depth_high, "}\n", sep = "")
cat("replications =", row$rep_k, "\n")
cat("NO dataset filtering\n")
cat("label =", row$label, "\n")

output_dir <- "results_ZIMN"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outfile_rds <- file.path(output_dir, paste0(row$label, ".rds"))
outfile_csv <- file.path(output_dir, paste0(row$label, "_full.csv"))

if(file.exists(outfile_rds)){
  cat("Skipping existing file:", outfile_rds, "\n")
  quit(save = "no")
}

res <- try(
  run_simulation_ZIMN(
    rep_k = row$rep_k,
    n = row$n,
    m = row$m,
    k = 2,
    zi_level = row$zi_level,
    depth_level = row$depth_level,
    init_method = "DS",
    seed_base = 3036000 + scenario_id * 10000
  ),
  silent = TRUE
)

if(inherits(res, "try-error")){
  cat("ZIMN scenario failed.\n")
  failed_file <- file.path(output_dir, paste0("FAILED_", row$label, ".rds"))
  saveRDS(list(scenario = as.list(row), error = as.character(res)), failed_file)
  stop(res)
}

## Metadata only; no replicate is screened or rejected.
res$design_zeta_low <- row$zeta_low
res$design_zeta_high <- row$zeta_high
res$design_depth_low <- row$depth_low
res$design_depth_high <- row$depth_high
res$design_depth_distribution <- "discrete_uniform"

saveRDS(res, outfile_rds)
write.csv(res, outfile_csv, row.names = FALSE)

cat("\nSaved:\n")
cat(outfile_rds, "\n")
cat(outfile_csv, "\n")
cat("Finished ZIMN scenario:", row$label, "\n")

}


# ==================== ZINB generator ====================
run_simulation_ZINB <- local({
  # Generator-specific implementation: generate_ZIDM_data
generate_ZIDM_data <- function(
    
  # ============================================================
  # Basic simulation dimensions
  # ============================================================
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.40)
  # high   : U(0.50, 0.70)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional user-specified true values
  # ============================================================
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter distributions from simulation table
  # ============================================================
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenario settings
  # If supplied, override zi_level/depth_level
  # ============================================================
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Structural-zero handling
  # ============================================================
  zero_inflation = TRUE,
  
  # "zero_row":
  # if all taxa structurally absent, Y_i = 0
  #
  # "resample":
  # redraw Z until at least one taxon survives
  all_zero_action = c("zero_row", "resample"),
  
  max_resample = 1000,
  
  # ============================================================
  # Output control
  # ============================================================
  return_Z = TRUE,
  return_lambda = TRUE,
  return_eta = TRUE,
  return_alpha_ij = TRUE,
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  all_zero_action <- match.arg(all_zero_action)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Scenario definitions
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.40),
      
      high = c(0.50, 0.70)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be a positive number.")
  }
  
  
  # ============================================================
  # 2. Generate observed covariates X
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 3. Generate latent variables
  #
  # u_i ~ N(0, I_k)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 4. Generate beta0
  #
  # beta_0j ~ U(-1,1)
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Generate covariate effects
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Generate factor loadings
  #
  # gamma_j ~ U(-2,2)
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Zero-inflation probabilities
  #
  # zeta_j ~ Uniform(...)
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 8. Linear predictor
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 9. Dirichlet/Gamma shape
  #
  # alpha_ij = exp(eta_ij)
  # ============================================================
  
  alpha_ij <- exp(eta)
  
  
  if (
    any(!is.finite(alpha_ij)) ||
    any(alpha_ij <= 0)
  ) {
    stop(
      "Non-finite or non-positive exp(eta) encountered."
    )
  }
  
  
  # ============================================================
  # 10. Gamma latent variables
  #
  # lambda0_ij ~ Gamma(exp(eta_ij), rate = 1)
  # ============================================================
  
  lambda0 <- matrix(
    
    rgamma(
      n * m,
      shape = as.vector(alpha_ij),
      rate = 1
    ),
    
    nrow = n,
    ncol = m
  )
  
  
  # ============================================================
  # 11. Composition BEFORE zero inflation
  #
  # theta0_ij =
  # lambda0_ij / sum_l lambda0_il
  # ============================================================
  
  lambda0_rowsum <- rowSums(lambda0)
  
  rho0_true <- lambda0 /
    lambda0_rowsum
  
  
  # ============================================================
  # 12. Structural-zero process
  #
  # zeta_j = structural-zero probability
  #
  # Therefore:
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  #
  # Z = 1 : taxon survives
  # Z = 0 : structural zero
  # ============================================================
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      rep(
        1 - zeta,
        each = n
      ),
      nrow = n,
      ncol = m
    )
    
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
    
  } else {
    
    Z <- matrix(
      1L,
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 13. Optional handling of rows with all structural zeros
  # ============================================================
  
  if (
    zero_inflation &&
    all_zero_action == "resample"
  ) {
    
    zero_rows <- which(
      rowSums(Z) == 0
    )
    
    iter <- 0
    
    while (length(zero_rows) > 0) {
      
      iter <- iter + 1
      
      if (iter > max_resample) {
        stop(
          "Exceeded max_resample for all-zero structural rows."
        )
      }
      
      for (i in zero_rows) {
        
        Z[i, ] <- rbinom(
          m,
          size = 1,
          prob = 1 - zeta
        )
      }
      
      zero_rows <- which(
        rowSums(Z) == 0
      )
    }
  }
  
  
  # ============================================================
  # 14. Apply structural zeros
  # ============================================================
  
  lambdaZI <- lambda0 * Z
  
  row_mass <- rowSums(lambdaZI)
  
  
  # ============================================================
  # 15. ZI composition
  # ============================================================
  
  rhoZI_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  
  positive_rows <- row_mass > 0
  
  
  rhoZI_true[positive_rows, ] <-
    lambdaZI[positive_rows, , drop = FALSE] /
    row_mass[positive_rows]
  
  
  # ============================================================
  # 16. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.integer(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  
  N <- as.integer(N)
  
  
  # ============================================================
  # 17. Multinomial observation
  #
  # If all taxa are structurally zero:
  #
  # Y_i = 0
  #
  # otherwise:
  #
  # Y_i ~ Multinomial(N_i, rhoZI_i)
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  
  for (i in seq_len(n)) {
    
    if (row_mass[i] > 0) {
      
      Y[i, ] <- as.integer(
        
        rmultinom(
          n = 1,
          size = N[i],
          prob = rhoZI_true[i, ]
        )
      )
    }
  }
  
  
  # ============================================================
  # 18. Useful simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  all_structural_zero_prop <-
    mean(rowSums(Z) == 0)
  
  
  # ============================================================
  # 19. Return EVERYTHING needed for parameter recovery
  # ============================================================
  
  out <- list(
    
    # --------------------------------------------------------
    # observed data
    # --------------------------------------------------------
    
    Y = Y,
    X = X,
    N = N,
    
    
    # --------------------------------------------------------
    # TRUE composition
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for model recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # TRUE latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # TRUE zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # true linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    alpha_ij = alpha_ij,
    
    
    # --------------------------------------------------------
    # scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    all_structural_zero_prop =
      all_structural_zero_prop
  )
  
  
  # ============================================================
  # Optional latent outputs
  # ============================================================
  
  if (return_Z) {
    
    out$Z <- Z
  }
  
  
  if (return_lambda) {
    
    out$lambda0 <- lambda0
    
    out$lambdaZI <- lambdaZI
  }
  
  
  if (!return_eta) {
    
    out$eta <- NULL
  }
  
  
  if (!return_alpha_ij) {
    
    out$alpha_ij <- NULL
  }
  
  
  return(out)
}

  # Generator-specific implementation: generate_ZIMN_data
generate_ZIMN_data <- function(
    
  # ============================================================
  # Basic dimensions
  # ============================================================
  
  n,
  m,
  k = 2,
  
  # ============================================================
  # Simulation scenarios
  #
  # Zero inflation:
  # low    : U(0.05, 0.15)
  # medium : U(0.20, 0.40)
  # high   : U(0.50, 0.70)
  #
  # Sequencing depth:
  # low    : Poisson(500)
  # medium : Poisson(2000)
  # high   : Poisson(5000)
  # ============================================================
  
  zi_level = c("low", "medium", "high"),
  depth_level = c("low", "medium", "high"),
  
  # ============================================================
  # Optional true parameter values
  #
  # IMPORTANT:
  # parameters are now defined directly in TAXON SPACE
  # ============================================================
  
  beta0 = NULL,             # length m
  beta = NULL,              # m x 2
  Gamma = NULL,             # m x k
  zeta = NULL,              # length m
  
  X = NULL,                 # n x 2
  U = NULL,                 # n x k
  N = NULL,                 # length n
  
  # ============================================================
  # Parameter ranges
  # ============================================================
  
  beta0_range = c(-1, 1),
  beta_range = c(-1, 1),
  gamma_range = c(-2, 2),
  
  # ============================================================
  # Optional custom scenarios
  # ============================================================
  
  zeta_range = NULL,
  depth_lambda = NULL,
  
  # ============================================================
  # Zero inflation
  # ============================================================
  
  zero_inflation = TRUE,
  
  # ============================================================
  # Random seed
  # ============================================================
  
  seed = NULL
){
  
  # ============================================================
  # 0. Basic setup
  # ============================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  zi_level <- match.arg(zi_level)
  depth_level <- match.arg(depth_level)
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (m <= 1 || m != as.integer(m)) {
    stop("m must be an integer greater than 1.")
  }
  
  if (k < 0 || k != as.integer(k)) {
    stop("k must be a nonnegative integer.")
  }
  
  
  # ============================================================
  # 1. Simulation scenario settings
  # ============================================================
  
  # ------------------------------------------------------------
  # Zero-inflation level
  # ------------------------------------------------------------
  
  if (is.null(zeta_range)) {
    
    zeta_range <- switch(
      
      zi_level,
      
      low = c(0.05, 0.15),
      
      medium = c(0.20, 0.40),
      
      high = c(0.50, 0.70)
    )
  }
  
  if (
    length(zeta_range) != 2 ||
    zeta_range[1] < 0 ||
    zeta_range[2] > 1 ||
    zeta_range[1] > zeta_range[2]
  ) {
    stop("zeta_range must be a valid probability interval.")
  }
  
  
  # ------------------------------------------------------------
  # Sequencing depth
  # ------------------------------------------------------------
  
  if (is.null(depth_lambda)) {
    
    depth_lambda <- switch(
      
      depth_level,
      
      low = 500,
      
      medium = 2000,
      
      high = 5000
    )
  }
  
  if (
    length(depth_lambda) != 1 ||
    !is.finite(depth_lambda) ||
    depth_lambda <= 0
  ) {
    stop("depth_lambda must be positive.")
  }
  
  
  # ============================================================
  # 2. Latent variables
  #
  # u_i ~ N(0, I_k)
  #
  # U: n x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(U)) {
      
      U <- matrix(
        rnorm(n * k),
        nrow = n,
        ncol = k
      )
      
    } else {
      
      U <- as.matrix(U)
      
      if (!all(dim(U) == c(n, k))) {
        stop("U must have dimension n x k.")
      }
    }
    
  } else {
    
    U <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 3. Observed covariates
  #
  # x_.1 ~ N(0,1)
  # x_.2 ~ Exp(1)
  #
  # X: n x 2
  # ============================================================
  
  if (is.null(X)) {
    
    X <- cbind(
      
      x1 = rnorm(
        n,
        mean = 0,
        sd = 1
      ),
      
      x2 = rexp(
        n,
        rate = 1
      )
    )
    
  } else {
    
    X <- as.matrix(X)
    
    if (!all(dim(X) == c(n, 2))) {
      stop("X must have dimension n x 2.")
    }
  }
  
  
  # ============================================================
  # 4. Species-specific intercepts
  #
  # beta_0j ~ U(-1,1)
  #
  # beta0: length m
  # ============================================================
  
  if (is.null(beta0)) {
    
    beta0 <- runif(
      m,
      min = beta0_range[1],
      max = beta0_range[2]
    )
    
  } else {
    
    beta0 <- as.numeric(beta0)
    
    if (length(beta0) != m) {
      stop("beta0 must have length m.")
    }
  }
  
  
  # ============================================================
  # 5. Covariate coefficients
  #
  # beta_1j, beta_2j ~ U(-1,1)
  #
  # beta: m x 2
  #
  # rows    = taxa
  # columns = covariates
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- matrix(
      
      runif(
        m * 2,
        min = beta_range[1],
        max = beta_range[2]
      ),
      
      nrow = m,
      ncol = 2
    )
    
    colnames(beta) <- c(
      "beta1",
      "beta2"
    )
    
  } else {
    
    beta <- as.matrix(beta)
    
    if (!all(dim(beta) == c(m, 2))) {
      stop("beta must have dimension m x 2.")
    }
  }
  
  
  # ============================================================
  # 6. Latent loadings
  #
  # gamma_j ~ U(-2,2)
  #
  # Gamma: m x k
  # ============================================================
  
  if (k > 0) {
    
    if (is.null(Gamma)) {
      
      Gamma <- matrix(
        
        runif(
          m * k,
          min = gamma_range[1],
          max = gamma_range[2]
        ),
        
        nrow = m,
        ncol = k
      )
      
    } else {
      
      Gamma <- as.matrix(Gamma)
      
      if (!all(dim(Gamma) == c(m, k))) {
        stop("Gamma must have dimension m x k.")
      }
    }
    
  } else {
    
    Gamma <- matrix(
      numeric(0),
      nrow = m,
      ncol = 0
    )
  }
  
  
  # ============================================================
  # 7. Linear predictor DIRECTLY in taxon space
  #
  # eta_ij =
  # beta_0j
  # + x_i^T beta_j
  # + u_i^T gamma_j
  #
  # eta: n x m
  # ============================================================
  
  eta <- matrix(
    beta0,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  # ------------------------------------------------------------
  # X beta
  #
  # X:        n x 2
  # t(beta):  2 x m
  #
  # result:   n x m
  # ------------------------------------------------------------
  
  eta <- eta +
    X %*% t(beta)
  
  
  # ------------------------------------------------------------
  # U Gamma
  #
  # U:         n x k
  # t(Gamma):  k x m
  #
  # result:    n x m
  # ------------------------------------------------------------
  
  if (k > 0) {
    
    eta <- eta +
      U %*% t(Gamma)
  }
  
  
  # ============================================================
  # 8. TRUE baseline composition
  #
  # rho_ij =
  #
  # exp(eta_ij)
  # ------------------
  # sum_l exp(eta_il)
  #
  # stable softmax
  # ============================================================
  
  rho0_true <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    row_max <- max(
      eta[i, ]
    )
    
    tmp <- exp(
      eta[i, ] - row_max
    )
    
    rho0_true[i, ] <-
      tmp / sum(tmp)
  }
  
  
  # ============================================================
  # 9. Zero-inflation probabilities
  #
  # zeta_j = P(structural zero)
  #
  # zeta_j ~ Uniform(...)
  #
  # zeta: length m
  # ============================================================
  
  if (is.null(zeta)) {
    
    zeta <- runif(
      m,
      min = zeta_range[1],
      max = zeta_range[2]
    )
    
  } else {
    
    zeta <- as.numeric(zeta)
    
    if (length(zeta) != m) {
      stop("zeta must have length m.")
    }
  }
  
  if (any(zeta < 0 | zeta > 1)) {
    stop("All zeta values must lie between 0 and 1.")
  }
  
  
  # ============================================================
  # 10. Structural-zero indicators
  #
  # zeta_j = P(structural zero)
  #
  # Define:
  #
  # Z_ij = 1 : taxon is PRESENT
  # Z_ij = 0 : structural zero
  #
  # therefore
  #
  # Z_ij ~ Bernoulli(1 - zeta_j)
  # ============================================================
  
  Z <- matrix(
    1L,
    nrow = n,
    ncol = m
  )
  
  if (zero_inflation) {
    
    prob_present <- matrix(
      
      rep(
        1 - zeta,
        each = n
      ),
      
      nrow = n,
      ncol = m
    )
    
    Z <- matrix(
      
      rbinom(
        n * m,
        size = 1,
        prob = as.vector(prob_present)
      ),
      
      nrow = n,
      ncol = m
    )
  }
  
  
  # ============================================================
  # 11. Apply structural zeros
  #
  # rhoZI_ij proportional to
  #
  # Z_ij * rho0_ij
  # ============================================================
  
  rhoZI_raw <-
    rho0_true * Z
  
  rs <- rowSums(
    rhoZI_raw
  )
  
  
  # ------------------------------------------------------------
  # Extremely rare case:
  # all taxa removed in one sample
  #
  # Activate one taxon so multinomial probabilities are valid
  # ------------------------------------------------------------
  
  bad <- which(
    rs == 0
  )
  
  if (length(bad) > 0) {
    
    for (i in bad) {
      
      j <- sample(
        seq_len(m),
        size = 1
      )
      
      Z[i, j] <- 1L
      
      rhoZI_raw[i, j] <-
        rho0_true[i, j]
    }
    
    rs <- rowSums(
      rhoZI_raw
    )
  }
  
  
  # ============================================================
  # 12. TRUE zero-inflated composition
  #
  # rhoZI_ij =
  #
  # Z_ij rho0_ij
  # --------------------
  # sum_l Z_il rho0_il
  # ============================================================
  
  rhoZI_true <-
    rhoZI_raw / rs
  
  
  # ============================================================
  # 13. Sequencing depth
  #
  # N_i ~ Poisson(depth_lambda)
  # ============================================================
  
  if (is.null(N)) {
    
    N <- rpois(
      n,
      lambda = depth_lambda
    )
    
  } else {
    
    N <- as.numeric(N)
    
    if (length(N) != n) {
      stop("N must have length n.")
    }
  }
  
  N <- as.integer(
    round(N)
  )
  
  if (any(N < 0)) {
    stop("N must be nonnegative.")
  }
  
  
  # ============================================================
  # 14. Multinomial sampling
  #
  # Y_i ~ Multinomial(
  #     N_i,
  #     rhoZI_true_i
  # )
  # ============================================================
  
  Y <- matrix(
    0L,
    nrow = n,
    ncol = m
  )
  
  for (i in seq_len(n)) {
    
    Y[i, ] <- as.integer(
      
      rmultinom(
        n = 1,
        size = N[i],
        prob = rhoZI_true[i, ]
      )
    )
  }
  
  
  # ============================================================
  # 15. Simulation summaries
  # ============================================================
  
  observed_zero_prop <-
    mean(Y == 0)
  
  structural_zero_prop <-
    mean(Z == 0)
  
  sampling_zero_prop <-
    mean(
      (Y == 0) &
        (Z == 1)
    )
  
  
  # ============================================================
  # 16. Return EVERYTHING needed for simulation evaluation
  # ============================================================
  
  return(list(
    
    # --------------------------------------------------------
    # Observed data
    # --------------------------------------------------------
    
    Y = Y,
    
    X = X,
    
    N = N,
    
    
    # --------------------------------------------------------
    # True compositions
    # --------------------------------------------------------
    
    rho0_true = rho0_true,
    
    rhoZI_true = rhoZI_true,
    
    # Primary composition for composition recovery
    rho_true = rho0_true,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    P_true = rho0_true,
    
    P_ZI = rhoZI_true,
    
    
    # --------------------------------------------------------
    # TRUE regression parameters
    #
    # NOW ALL PARAMETERS ARE IN TAXON SPACE
    # --------------------------------------------------------
    
    beta0 = beta0,
    
    beta = beta,
    
    beta1 = beta[, 1],
    
    beta2 = beta[, 2],
    
    
    # --------------------------------------------------------
    # True latent variables / loadings
    # --------------------------------------------------------
    
    U = U,
    
    Gamma = Gamma,
    
    
    # --------------------------------------------------------
    # Compatibility aliases
    # --------------------------------------------------------
    
    Beta = t(beta),
    
    Lambda = Gamma,
    
    
    # --------------------------------------------------------
    # True zero-inflation parameters
    # --------------------------------------------------------
    
    zeta = zeta,
    
    lgt_zeta = qlogis(zeta),
    
    
    # --------------------------------------------------------
    # Linear predictor
    # --------------------------------------------------------
    
    eta = eta,
    
    eta_taxon = eta,
    
    
    # --------------------------------------------------------
    # Scenario information
    # --------------------------------------------------------
    
    zi_level = zi_level,
    
    zeta_range = zeta_range,
    
    depth_level = depth_level,
    
    depth_lambda = depth_lambda,
    
    zero_inflation = zero_inflation,
    
    
    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------
    
    observed_zero_prop =
      observed_zero_prop,
    
    structural_zero_prop =
      structural_zero_prop,
    
    sampling_zero_prop =
      sampling_zero_prop,
    
    
    # --------------------------------------------------------
    # Structural-zero indicator
    # --------------------------------------------------------
    
    Z = Z
  ))
}

  # Generator-specific implementation: fit_gllvm_ZINB
fit_gllvm_ZINB <- function(
    Y,
    X = NULL,
    k = 2,
    seed = 123
){
  
  Y <- as.matrix(Y)
  
  n <- nrow(Y)
  m <- ncol(Y)
  
  if(!is.null(X)){
    
    X <- as.matrix(X)
    
    if(nrow(X) != n){
      stop("X must have the same number of rows as Y.")
    }
  }
  
  
  ##########################################################
  ## 1. Fit ZINB-GLLVM
  ##
  ## Improved convergence:
  ## - more iterations
  ## - multiple initial values
  ## - residual-based starting values
  ##########################################################
  
  fit <- gllvm::gllvm(
    y = Y,
    X = X,
    family = "ZINB",
    num.lv = k,
    row.eff = "fixed",
    method = "VA",
    sd.errors = FALSE,
    seed = seed,
    
    control = list(
      maxit = 20000,
      max.iter = 20000,
      reltol = 1e-10
    ),
    
    control.start = list(
      starting.val = "res",
      n.init = 3
    )
  )
  
  
  ##########################################################
  ## 2. Extract fitted model parameters
  ##########################################################
  
  # species-specific intercepts beta_0j
  beta0_hat <- fit$params$beta0
  
  # environmental coefficients beta_j
  beta_hat <- fit$params$Xcoef
  
  # latent scores u_i
  u_hat <- as.matrix(
    gllvm::getLV(fit)
  )
  
  # properly scaled loading vectors gamma_j
  gamma_hat <- as.matrix(
    gllvm::getLoadings(fit)
  )
  
  # fixed sample / row effects alpha_i
  # sample 1 is the reference, alpha_1 = 0
  alpha_hat <- c(
    0,
    fit$params$row.params.fixed
  )
  
  
  ##########################################################
  ## 3. Zero-inflation and dispersion parameters
  ##########################################################
  
  # ZI probability p_j
  # This corresponds to zeta_j in our notation
  zeta_hat <- fit$params$phi
  
  # NB dispersion parameter
  dispersion_hat <- fit$params$ZINB.phi
  
  # inverse NB dispersion
  inv_dispersion_hat <- fit$params$ZINB.inv.phi
  
  
  ##########################################################
  ## 4. Variational covariance matrices A_i
  ##
  ## A_hat[i,,] = A_i
  ## dimension n x k x k
  ##########################################################
  
  A_hat <- fit$A
  
  
  ##########################################################
  ## 5. Construct fitted linear predictor eta
  ##
  ## eta_ij =
  ## alpha_i
  ## + beta_0j
  ## + x_i^T beta_j
  ## + u_i^T gamma_j
  ##########################################################
  
  eta_hat <- matrix(
    beta0_hat,
    nrow = n,
    ncol = m,
    byrow = TRUE
  )
  
  
  ##########################################################
  ## add row effect alpha_i
  ##########################################################
  
  eta_hat <- eta_hat +
    matrix(
      alpha_hat,
      nrow = n,
      ncol = m,
      byrow = FALSE
    )
  
  
  ##########################################################
  ## add environmental contribution X beta
  ##########################################################
  
  if(!is.null(X) && !is.null(beta_hat)){
    
    beta_mat <- as.matrix(beta_hat)
    
    # usual gllvm orientation: m x q
    if(
      nrow(beta_mat) == m &&
      ncol(beta_mat) == ncol(X)
    ){
      
      eta_hat <- eta_hat +
        X %*% t(beta_mat)
      
    } else if(
      nrow(beta_mat) == ncol(X) &&
      ncol(beta_mat) == m
    ){
      
      eta_hat <- eta_hat +
        X %*% beta_mat
      
    } else {
      
      stop(
        "Dimension of fit$params$Xcoef is not compatible with X."
      )
    }
  }
  
  
  ##########################################################
  ## add latent contribution U Gamma
  ##########################################################
  
  # convert loadings to k x m
  if(
    nrow(gamma_hat) == m &&
    ncol(gamma_hat) == k
  ){
    
    Gamma_mat <- t(gamma_hat)
    
  } else if(
    nrow(gamma_hat) == k &&
    ncol(gamma_hat) == m
  ){
    
    Gamma_mat <- gamma_hat
    
  } else {
    
    stop(
      "Dimension of getLoadings(fit) is not compatible with k and m."
    )
  }
  
  eta_hat <- eta_hat +
    u_hat %*% Gamma_mat
  
  
  ##########################################################
  ## 6. Variational quadratic correction
  ##
  ## cQ_ij =
  ## 0.5 * gamma_j^T A_i gamma_j
  ##########################################################
  
  cQ_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    Ai <- A_hat[i, , ]
    
    for(j in seq_len(m)){
      
      gamma_j <- Gamma_mat[, j]
      
      cQ_hat[i, j] <-
        0.5 *
        as.numeric(
          t(gamma_j) %*%
            Ai %*%
            gamma_j
        )
    }
  }
  
  
  ##########################################################
  ## 7. Variationally corrected predictor
  ##########################################################
  
  eta_corrected_hat <-
    eta_hat + cQ_hat
  
  
  ##########################################################
  ## 8. Corrected fitted composition rho
  ##########################################################
  
  rho_hat <- matrix(
    0,
    nrow = n,
    ncol = m
  )
  
  for(i in seq_len(n)){
    
    row_max <- max(
      eta_corrected_hat[i, ]
    )
    
    tmp <- exp(
      eta_corrected_hat[i, ] -
        row_max
    )
    
    rho_hat[i, ] <-
      tmp / sum(tmp)
  }
  
  
  ##########################################################
  ## 9. Return everything needed for simulation
  ##########################################################
  
  return(list(
    
    alpha_hat = alpha_hat,
    
    beta0_hat = beta0_hat,
    
    beta_hat = beta_hat,
    
    u_hat = u_hat,
    
    gamma_hat = gamma_hat,
    
    zeta_hat = zeta_hat,
    
    dispersion_hat = dispersion_hat,
    
    inv_dispersion_hat = inv_dispersion_hat,
    
    A_hat = A_hat,
    
    cQ_hat = cQ_hat,
    
    eta_hat = eta_hat,
    
    eta_corrected_hat = eta_corrected_hat,
    
    rho_hat = rho_hat,
    
    theta_raw = fit$params$theta,
    
    sigma_lv = fit$params$sigma.lv,
    
    ordination = u_hat,
    
    loadings = gamma_hat,
    
    model = fit
  ))
}

run_simulation_ZINB <- function(
    rep_k = 500,
    n = 50,
    m = 100,
    k = 2,
    zi_level = c("low", "medium", "high"),
    depth_level = c("low", "medium", "high"),
    init_method = "DS",
    phi = 1,
    seed_base = 2026,
    compute_CI = TRUE,
    CI_level = 0.95,
    include_ZIDM_paper = TRUE,
    ZIDM_paper_iterations = 20000,
    ZIDM_paper_thin = 10,
    ZIDM_paper_burnin = 0.5
){
  
  zi_level <- match.arg(zi_level)
  
  depth_level <- match.arg(depth_level)
  
  
  ############################################################
  ## ZI range
  ############################################################
  
  zeta_range <- switch(
    
    zi_level,
    
    low = c(0.05, 0.15),
    
    medium = c(0.20, 0.30),
    
    high = c(0.35, 0.45)
  )
  
  
  ############################################################
  ## ZINB row effect / depth setting
  ############################################################
  
  if(n == 50 && m == 100){
    alpha_i_value <- switch(
      zi_level,
      low = 0.25,
      medium = 1.00,
      high = 2.25
    )
  } else if(n == 100 && m == 50){
    alpha_i_value <- switch(
      zi_level,
      low = 0.25,
      medium = 1.00,
      high = 2.00
    )
  } else {
    stop("No calibrated ZINB alpha_i for this (n,m) setting.")
  }
  
  
  ############################################################
  ## TMB models: NO STR
  ############################################################
  
  tmb_model_list <- list(
    
    list(
      name = "ZIMN_soft",
      model = "ZIMN",
      extra = c(0,0)
    ),
    
    list(
      name = "ZIMN_hard",
      model = "ZIMN",
      extra = c(0,1)
    ),
    
    list(
      name = "ZIDM_soft",
      model = "ZIDM",
      extra = c(0,0)
    ),
    
    list(
      name = "ZIDM_hard",
      model = "ZIDM",
      extra = c(0,1)
    )
  )
  
  
  methods_all <- c(
    "ZINB_gllvm",
    "ZIMN_soft",
    "ZIMN_hard",
    "ZIDM_soft",
    "ZIDM_hard"
  )
  
  if(include_ZIDM_paper){
    methods_all <- c(
      methods_all,
      "ZIDM_paper"
    )
  }
  
  methods_all <- c(
    methods_all,
    "ZR"
  )
  
  
  ############################################################
  ## SAFE
  ############################################################
  
  safe_fit <- function(expr){
    
    out <- try(
      expr,
      silent = TRUE
    )
    
    if(inherits(out, "try-error")){
      return(NULL)
    }
    
    out
  }
  
  
  safe_eval <- function(expr){
    
    out <- try(
      expr,
      silent = TRUE
    )
    
    if(inherits(out, "try-error")){
      return(NULL)
    }
    
    out
  }
  
  
  ############################################################
  ## convergence checks
  ############################################################
  
  valid_TMB <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    
    if(is.null(fit$convergence)) return(FALSE)
    
    if(fit$convergence != 0) return(FALSE)
    
    if(is.null(fit$rho_hat)) return(FALSE)
    
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    TRUE
  }
  
  
  valid_gllvm <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    
    if(is.null(fit$model)) return(FALSE)
    
    if(is.null(fit$model$convergence)) return(FALSE)
    
    if(!isTRUE(fit$model$convergence)) return(FALSE)
    
    if(is.null(fit$rho_hat)) return(FALSE)
    
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    TRUE
  }
  
  
  ############################################################
  ## paper ZIDM validity
  ############################################################
  
  valid_ZIDM_paper <- function(fit){
    
    if(is.null(fit)) return(FALSE)
    if(is.null(fit$rho_hat)) return(FALSE)
    if(!is.matrix(fit$rho_hat)) return(FALSE)
    if(!all(dim(fit$rho_hat) == c(n, m))) return(FALSE)
    if(any(!is.finite(fit$rho_hat))) return(FALSE)
    
    if(is.null(fit$beta0_hat)) return(FALSE)
    if(length(fit$beta0_hat) != m) return(FALSE)
    if(any(!is.finite(fit$beta0_hat))) return(FALSE)
    
    if(is.null(fit$beta_hat)) return(FALSE)
    if(!all(dim(fit$beta_hat) == c(m, 2))) return(FALSE)
    if(any(!is.finite(fit$beta_hat))) return(FALSE)
    
    if(is.null(fit$zeta_hat)) return(FALSE)
    if(length(fit$zeta_hat) != m) return(FALSE)
    if(any(!is.finite(fit$zeta_hat))) return(FALSE)
    
    TRUE
  }
  
  
  ############################################################
  ## NA row
  ############################################################
  
    ############################################################
  ## CI summary helper
  ############################################################

  empty_ci_summary <- function(requested = compute_CI){
    list(
      ci_success = if(requested) FALSE else NA,
      beta0_Coverage = NA_real_,
      beta1_Coverage = NA_real_,
      beta2_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta1_meanSE = NA_real_,
      beta2_meanSE = NA_real_,
      beta0_meanWidth = NA_real_,
      beta1_meanWidth = NA_real_,
      beta2_meanWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta1_nValid = NA_integer_,
      beta2_nValid = NA_integer_,
      beta0_validProp = NA_real_,
      beta1_validProp = NA_real_,
      beta2_validProp = NA_real_
    )
  }


  ############################################################
  ## NA row
  ############################################################
  
  make_na_row <- function(
    r,
    method,
    runtime = NA_real_,
    convergence = NA,
    ci_runtime = NA_real_,
    observed_zero_prop = NA_real_,
    structural_zero_prop = NA_real_,
    sampling_zero_prop = NA_real_,
    generation_runtime = NA_real_
  ){

    runtime_total <- if(
      is.finite(runtime) && is.finite(ci_runtime)
    ){
      runtime + ci_runtime
    } else if(is.finite(runtime)){
      runtime
    } else {
      NA_real_
    }

    data.frame(
      rep = r,
      method = method,
      success = FALSE,
      convergence = convergence,

      runtime = runtime,
      ci_runtime = ci_runtime,
      runtime_total = runtime_total,
      ci_success = if(compute_CI && method != "ZR") FALSE else NA,

      Frobenius = NA_real_,
      KL = NA_real_,
      Shannon_MSE = NA_real_,
      Simpson_MSE = NA_real_,
      Procrustes = NA_real_,

      beta0_Bias = NA_real_,
      beta0_AbsBias = NA_real_,
      beta0_RMSE = NA_real_,

      beta1_Bias = NA_real_,
      beta1_AbsBias = NA_real_,
      beta1_RMSE = NA_real_,

      beta2_Bias = NA_real_,
      beta2_AbsBias = NA_real_,
      beta2_RMSE = NA_real_,

      zeta_Bias = NA_real_,
      zeta_AbsBias = NA_real_,
      zeta_RMSE = NA_real_,

      beta0_Coverage = NA_real_,
      beta0_meanSE = NA_real_,
      beta0_meanCIWidth = NA_real_,
      beta0_nValid = NA_integer_,
      beta0_validProp = NA_real_,

      beta1_Coverage = NA_real_,
      beta1_meanSE = NA_real_,
      beta1_meanCIWidth = NA_real_,
      beta1_nValid = NA_integer_,
      beta1_validProp = NA_real_,

      beta2_Coverage = NA_real_,
      beta2_meanSE = NA_real_,
      beta2_meanCIWidth = NA_real_,
      beta2_nValid = NA_integer_,
      beta2_validProp = NA_real_,

      observed_zero_prop = observed_zero_prop,
      structural_zero_prop = structural_zero_prop,
      sampling_zero_prop = sampling_zero_prop,
      generation_runtime = generation_runtime,

      stringsAsFactors = FALSE
    )
  }

  ############################################################
  ## safe bind by column name
  ## Different method rows are padded with NA before rbind.
  ############################################################

  rbind_fill <- function(x){

    x <- x[
      !vapply(
        x,
        is.null,
        logical(1)
      )
    ]

    if(length(x) == 0L){
      return(NULL)
    }

    all_cols <- unique(
      unlist(
        lapply(
          x,
          names
        ),
        use.names = FALSE
      )
    )

    x <- lapply(
      x,
      function(z){

        miss <- setdiff(
          all_cols,
          names(z)
        )

        if(length(miss) > 0L){
          for(nm in miss){
            z[[nm]] <- NA
          }
        }

        z[
          ,
          all_cols,
          drop = FALSE
        ]
      }
    )

    out <- do.call(
      rbind,
      x
    )

    rownames(out) <- NULL

    out
  }
  
  
  ############################################################
  ## storage
  ############################################################
  
  results <- vector(
    "list",
    rep_k
  )
  
  
  ############################################################
  ## repetitions
  ############################################################
  
  for(r in seq_len(rep_k)){
    
    cat(
      "\n========================================\n"
    )
    
    cat(
      "ZINB replicate:",
      r,
      "/",
      rep_k,
      "\n"
    )
    
    cat(
      "n =", n,
      "| m =", m,
      "| zi =", zi_level,
      "| depth =", depth_level,
      "| alpha_i =", alpha_i_value,
      "\n"
    )
    
    cat(
      "========================================\n"
    )
    
    
    rep_seed <- seed_base + r
    
    
    ##########################################################
    ## FULL REP PROTECTION
    ##########################################################
    
    rep_result <- try({
      
      ########################################################
      ## 1. Generate ZINB data
      ########################################################
      
      t0_generate <- Sys.time()
      
      
      dat <- generate_ZINB_data(
        
        n = n,
        
        m = m,
        
        k = k,
        
        p = 2,
        
        zeta_range = zeta_range,
        
        alpha_i = alpha_i_value,
        
        phi = phi,
        
        zero_inflation = TRUE,
        
        seed = rep_seed
      )
      
      
      generation_runtime <- as.numeric(
        Sys.time() - t0_generate,
        units = "secs"
      )
      
      
      ########################################################
      ## observed data
      ########################################################
      
      Y <- dat$Y
      
      X <- dat$X
      
      
      ########################################################
      ## true quantities
      ########################################################
      
      rho_true <- dat$truth$rho
      
      U_true <- dat$truth$U
      
      beta0_true <- dat$truth$beta0
      
      
      beta_true <- as.matrix(
        dat$truth$beta
      )
      
      
      beta1_true <-
        beta_true[,1]
      
      
      beta2_true <-
        beta_true[,2]
      
      
      zeta_true <-
        dat$truth$zeta
      
      
      ########################################################
      ## fitted models
      ########################################################
      
      fit_results <- list()
      
      
      ########################################################
      ## 2. gllvm ZINB
      ########################################################
      
      cat("  fitting ZINB_gllvm ...\n")
      
      
      t0 <- Sys.time()
      
      
      fit_g <- safe_fit(
        
        fit_gllvm_ZINB(
          
          Y = Y,
          
          X = X,
          
          k = k,
          
          seed = rep_seed
        )
      )
      
      
      runtime_g <- as.numeric(
        Sys.time() - t0,
        units = "secs"
      )
      
      
      if(valid_gllvm(fit_g)){
        
        fit_results$ZINB_gllvm <- list(
          
          fit = fit_g,
          
          runtime = runtime_g,
          
          success = TRUE,
          
          convergence = TRUE,
          
          type = "gllvm"
        )
        
      } else {
        
        fit_results$ZINB_gllvm <- list(
          
          fit = NULL,
          
          runtime = runtime_g,
          
          success = FALSE,
          
          convergence = FALSE,
          
          type = "gllvm"
        )
      }
      
      
      ########################################################
      ## 3. Four TMB models
      ########################################################
      
      for(mod in tmb_model_list){
        
        cat(
          "  fitting",
          mod$name,
          "...\n"
        )
        
        
        t0 <- Sys.time()
        
        
        fit_tmp <- safe_fit(
          
          fit_TMB(
            
            Y = Y,
            
            X = X,
            
            k = k,
            
            model = mod$model,
            
            extra = mod$extra,
            
            init_method = init_method
          )
        )
        
        
        runtime_tmp <- as.numeric(
          Sys.time() - t0,
          units = "secs"
        )
        
        
        if(valid_TMB(fit_tmp)){
          
          fit_results[[mod$name]] <- list(
            
            fit = fit_tmp,
            
            runtime = runtime_tmp,
            
            success = TRUE,
            
            convergence =
              fit_tmp$convergence,
            
            type = "TMB"
          )
          
        } else {
          
          fit_results[[mod$name]] <- list(
            
            fit = NULL,
            
            runtime = runtime_tmp,
            
            success = FALSE,
            
            convergence =
              if(!is.null(fit_tmp$convergence)){
                fit_tmp$convergence
              } else {
                NA
              },
            
            type = "TMB"
          )
        }
      }
      
      
      ########################################################
      ## 4. Koslovsky paper ZIDM
      ## standard regression, no BVS, no CI
      ########################################################
      
      if(include_ZIDM_paper){
        
        cat("  fitting ZIDM_paper ...\n")
        
        t0 <- Sys.time()
        
        fit_paper <- safe_fit(
          fit_ZIDM_paper(
            Y = Y,
            X = X,
            X_theta = NULL,
            iterations = ZIDM_paper_iterations,
            thin = ZIDM_paper_thin,
            burnin_prop = ZIDM_paper_burnin,
            seed = rep_seed
          )
        )
        
        runtime_paper <- as.numeric(
          Sys.time() - t0,
          units = "secs"
        )
        
        if(valid_ZIDM_paper(fit_paper)){
          
          fit_results$ZIDM_paper <- list(
            fit = fit_paper,
            runtime = runtime_paper,
            success = TRUE,
            convergence = TRUE,
            type = "ZIDM_paper"
          )
          
        } else {
          
          fit_results$ZIDM_paper <- list(
            fit = NULL,
            runtime = runtime_paper,
            success = FALSE,
            convergence = FALSE,
            type = "ZIDM_paper"
          )
        }
      }
      
      
      ########################################################
      ## 5. ZR
      ########################################################
      
      cat("  fitting ZR ...\n")
      
      
      t0 <- Sys.time()
      
      
      fit_zr <- safe_fit(
        fit_ZR(Y)
      )
      
      
      runtime_zr <- as.numeric(
        Sys.time() - t0,
        units = "secs"
      )
      
      
      if(
        !is.null(fit_zr) &&
        !is.null(fit_zr$rho_hat) &&
        all(is.finite(fit_zr$rho_hat))
      ){
        
        fit_results$ZR <- list(
          
          fit = fit_zr,
          
          runtime = runtime_zr,
          
          success = TRUE,
          
          convergence = TRUE,
          
          type = "ZR"
        )
        
      } else {
        
        fit_results$ZR <- list(
          
          fit = NULL,
          
          runtime = runtime_zr,
          
          success = FALSE,
          
          convergence = FALSE,
          
          type = "ZR"
        )
      }
      
      
      ########################################################
      ## 6. Evaluation
      ########################################################
      
      ########################################################
      ## Confidence intervals / credible intervals
      ########################################################

      ci_results <- list()

      for(method in methods_all){

        ci_results[[method]] <- list(
          ci = NULL,
          summary = empty_ci_summary(compute_CI),
          runtime = NA_real_
        )

        if(!compute_CI){
          next
        }

        if(method == "ZR"){
          ci_results[[method]] <- list(
            ci = NULL,
            summary = empty_ci_summary(FALSE),
            runtime = NA_real_
          )
          next
        }

        fit_obj <- fit_results[[method]]

        if(is.null(fit_obj) || !isTRUE(fit_obj$success) || is.null(fit_obj$fit)){
          next
        }

        t0_ci <- Sys.time()
        ci_tmp <- NULL

        if(fit_obj$type == "gllvm"){

          ci_tmp <- safe_eval(
            compute_gllvm_ZINB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )

        } else if(fit_obj$type == "TMB"){

          ci_tmp <- safe_eval(
            compute_TMB_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )

        } else if(fit_obj$type == "ZIDM_paper"){

          ci_tmp <- safe_eval(
            compute_ZIDM_paper_beta_CI(
              fit = fit_obj$fit,
              beta0_true = beta0_true,
              beta1_true = beta1_true,
              beta2_true = beta2_true,
              level = CI_level
            )
          )
        }

        runtime_ci <- as.numeric(Sys.time() - t0_ci, units = "secs")

        if(is.null(ci_tmp)){
          ci_sum <- empty_ci_summary(TRUE)
        } else {
          ci_sum <- safe_eval(summarize_beta_CI(ci_tmp))
          if(is.null(ci_sum)){
            ci_sum <- empty_ci_summary(TRUE)
          }
        }

        ci_results[[method]] <- list(
          ci = ci_tmp,
          summary = ci_sum,
          runtime = runtime_ci
        )
      }


      out_rep <- list()
      
      
      for(method in methods_all){
        
        obj <- fit_results[[method]]
        
        
        ######################################################
        ## failure
        ######################################################
        
        if(
          is.null(obj) ||
          !isTRUE(obj$success) ||
          is.null(obj$fit)
        ){
          
          out_rep[[method]] <-
            make_na_row(
              
              r = r,
              
              method = method,
              
              runtime =
                if(!is.null(obj$runtime)){
                  obj$runtime
                } else {
                  NA_real_
                }
            )
          
          next
        }
        
        
        method_eval <- try({

        fit <- obj$fit

        ci_obj <- ci_results[[method]]
        ci_sum <- ci_obj$summary
        ci_runtime <- ci_obj$runtime

        runtime_total <- if(is.finite(ci_runtime)){
          obj$runtime + ci_runtime
        } else {
          obj$runtime
        }
        
        
        ######################################################
        ## composition
        ######################################################
        
        comp <- safe_eval(
          
          evaluate_composition(
            
            rho_hat = fit$rho_hat,
            
            rho_true = rho_true
          )
        )
        
        
        if(is.null(comp)){
          
          comp <- list(
            Frobenius = NA,
            KL = NA,
            Shannon_MSE = NA,
            Simpson_MSE = NA
          )
        }
        
        
        ######################################################
        ## ordination
        ######################################################
        
        proc <- NA_real_
        
        
        if(method != "ZR" && method != "ZIDM_paper"){
          
          proc_tmp <- safe_eval(
            
            compute_procrustes(
              
              true_ord = U_true,
              
              est_ord = fit$ordination
            )
          )
          
          
          if(
            !is.null(proc_tmp) &&
            length(proc_tmp) == 1 &&
            is.finite(proc_tmp)
          ){
            
            proc <- proc_tmp
          }
        }
        
        
        ######################################################
        ## parameter defaults
        ######################################################
        
        empty_met <- list(
          Bias = NA,
          AbsBias = NA,
          RMSE = NA
        )
        
        
        beta0_met <- empty_met
        beta1_met <- empty_met
        beta2_met <- empty_met
        zeta_met <- empty_met
        
        
        ######################################################
        ## parameter recovery
        ######################################################
        
        if(method != "ZR"){
          
          
          ####################################################
          ## beta0
          ####################################################
          
          tmp <- safe_eval(
            
            evaluate_parameter(
              
              est = fit$beta0_hat,
              
              true = beta0_true
            )
          )
          
          if(!is.null(tmp)){
            beta0_met <- tmp
          }
          
          
          ####################################################
          ## beta
          ####################################################
          
          if(obj$type %in% c("gllvm", "ZIDM_paper")){
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[,1],
                
                true = beta1_true
              )
            )
            
            if(!is.null(tmp)){
              beta1_met <- tmp
            }
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[,2],
                
                true = beta2_true
              )
            )
            
            if(!is.null(tmp)){
              beta2_met <- tmp
            }
            
            
          } else if(obj$type == "TMB"){
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[1,],
                
                true = beta1_true
              )
            )
            
            if(!is.null(tmp)){
              beta1_met <- tmp
            }
            
            
            tmp <- safe_eval(
              
              evaluate_parameter(
                
                est = fit$beta_hat[2,],
                
                true = beta2_true
              )
            )
            
            if(!is.null(tmp)){
              beta2_met <- tmp
            }
          }
          
          
          ####################################################
          ## zero inflation
          ####################################################
          
          tmp <- safe_eval(
            
            evaluate_parameter(
              
              est = fit$zeta_hat,
              
              true = zeta_true
            )
          )
          
          
          if(!is.null(tmp)){
            zeta_met <- tmp
          }
        }
        
        
        ######################################################
        ## result row
        ######################################################
        
        method_row <- data.frame(
          
          rep = r,
          
          method = method,
          
          success = TRUE,
          
          convergence =
            obj$convergence,
          
          runtime =
            obj$runtime,
          
          ci_runtime = ci_runtime,
          
          runtime_total = runtime_total,
          
          ci_success = ci_sum$ci_success,
          
          
          Frobenius =
            comp$Frobenius,
          
          KL =
            comp$KL,
          
          Shannon_MSE =
            comp$Shannon_MSE,
          
          Simpson_MSE =
            comp$Simpson_MSE,
          
          
          Procrustes =
            proc,
          
          
          beta0_Bias =
            beta0_met$Bias,
          
          beta0_AbsBias =
            beta0_met$AbsBias,
          
          beta0_RMSE =
            beta0_met$RMSE,
          
          
          beta1_Bias =
            beta1_met$Bias,
          
          beta1_AbsBias =
            beta1_met$AbsBias,
          
          beta1_RMSE =
            beta1_met$RMSE,
          
          
          beta2_Bias =
            beta2_met$Bias,
          
          beta2_AbsBias =
            beta2_met$AbsBias,
          
          beta2_RMSE =
            beta2_met$RMSE,
          
          
          zeta_Bias =
            zeta_met$Bias,
          
          zeta_AbsBias =
            zeta_met$AbsBias,
          
          zeta_RMSE =
            zeta_met$RMSE,
          
          
          beta0_Coverage = ci_sum$beta0_Coverage,
          
          
          beta0_meanSE = ci_sum$beta0_meanSE,
          
          
          beta0_meanCIWidth = ci_sum$beta0_meanWidth,
          
          
          beta0_nValid = ci_sum$beta0_nValid,
          
          
          beta0_validProp = ci_sum$beta0_validProp,
          
          
          beta1_Coverage = ci_sum$beta1_Coverage,
          
          
          beta1_meanSE = ci_sum$beta1_meanSE,
          
          
          beta1_meanCIWidth = ci_sum$beta1_meanWidth,
          
          
          beta1_nValid = ci_sum$beta1_nValid,
          
          
          beta1_validProp = ci_sum$beta1_validProp,
          
          
          beta2_Coverage = ci_sum$beta2_Coverage,
          
          
          beta2_meanSE = ci_sum$beta2_meanSE,
          
          
          beta2_meanCIWidth = ci_sum$beta2_meanWidth,
          
          
          beta2_nValid = ci_sum$beta2_nValid,
          
          
          beta2_validProp = ci_sum$beta2_validProp,

          
          
          observed_zero_prop =
            dat$observed_zero_prop,
          
          structural_zero_prop =
            dat$structural_zero_prop,
          
          sampling_zero_prop =
            dat$sampling_zero_prop,
          
          generation_runtime =
            generation_runtime,
          
          
          stringsAsFactors = FALSE
        )

        method_row

        }, silent = TRUE)

        if(inherits(method_eval, "try-error")){
          cat("  method", method, "evaluation FAILED -> saving NA for this method only\n")
          cat("    ", as.character(method_eval), "\n")
          out_rep[[method]] <- make_na_row(
            r = r, method = method,
            runtime = if(!is.null(obj$runtime)) obj$runtime else NA_real_,
            convergence = if(!is.null(obj$convergence)) obj$convergence else NA,
            observed_zero_prop = dat$observed_zero_prop,
            structural_zero_prop = dat$structural_zero_prop,
            sampling_zero_prop = dat$sampling_zero_prop,
            generation_runtime = generation_runtime
          )
        } else {
          out_rep[[method]] <- method_eval
        }
      }
      
      
      rbind_fill(out_rep)
      
    }, silent = TRUE)
    
    
    ##########################################################
    ## whole replicate failure
    ##########################################################
    
    if(inherits(rep_result, "try-error")){
      
      cat(
        "  replicate",
        r,
        "FAILED completely\n"
      )
      
      
      rep_result <- rbind_fill(
        lapply(
                  
                  methods_all,
                  
                  function(mm){
                    
                    make_na_row(
                      r = r,
                      method = mm
                    )
                  }
                )
      )
    }
    
    
    results[[r]] <- rep_result
  }
  
  
  ############################################################
  ## combine
  ############################################################
  
  results_df <- rbind_fill(results)
  
  
  rownames(results_df) <- NULL
  
  
  ############################################################
  ## scenario information
  ############################################################
  
  results_df$DGP <- "ZINB"
  
  results_df$n <- n
  results_df$m <- m
  results_df$k <- k
  
  results_df$zi_level <- zi_level
  
  results_df$depth_level <- depth_level
  
  results_df$alpha_i <- alpha_i_value
  
  results_df$phi_true <- phi
  
  results_df$ZIDM_paper_iterations <- if(include_ZIDM_paper) ZIDM_paper_iterations else NA_integer_
  results_df$ZIDM_paper_thin <- if(include_ZIDM_paper) ZIDM_paper_thin else NA_integer_
  results_df$ZIDM_paper_burnin <- if(include_ZIDM_paper) ZIDM_paper_burnin else NA_real_
  
  
  results_df$compute_CI <- compute_CI
  results_df$CI_level <- if(compute_CI) CI_level else NA_real_

  id_cols <- c(
    "DGP",
    "n",
    "m",
    "k",
    "zi_level",
    "depth_level",
    "alpha_i",
    "phi_true",
    "rep",
    "method"
  )
  
  
  results_df <- results_df[
    ,
    c(
      id_cols,
      setdiff(
        names(results_df),
        id_cols
      )
    )
  ]
  
  
  return(results_df)
}
})

# Original ZINB scenario table and SLURM output handling.
run_scenario_ZINB <- function() {



############################################################
## Puhti scenario control
## DGP = ZINB
## FINAL calibrated design
## NO dataset filtering
############################################################

scenario_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(scenario_id)) scenario_id <- 1
cat("SLURM_ARRAY_TASK_ID =", scenario_id, "\n")

scenario_table <- data.frame(
  scenario = 1:6,
  n = c(50,50,50,100,100,100),
  m = c(100,100,100,50,50,50),
  zi_level = c("low","medium","high","low","medium","high"),
  depth_level = c("low","medium","high","low","medium","high"),
  zeta_low = c(0.05,0.20,0.35,0.05,0.20,0.35),
  zeta_high = c(0.15,0.30,0.45,0.15,0.30,0.45),
  alpha_i = c(0.25,1.00,2.25,0.25,1.00,2.00),
  stringsAsFactors = FALSE
)

scenario_table$rep_k <- 100

scenario_table$label <- paste0(
  "ZINB_n", scenario_table$n,
  "_m", scenario_table$m,
  "_ZI", scenario_table$zi_level,
  "_ALPHA", scenario_table$alpha_i,
  "_R", scenario_table$rep_k
)

print(scenario_table)

if(scenario_id < 1 || scenario_id > nrow(scenario_table)){
  stop("scenario_id out of range.")
}

row <- scenario_table[scenario_id, , drop = FALSE]

cat("\n========================================\n")
cat("Running FINAL ZINB scenario\n")
cat("========================================\n")
cat("scenario =", row$scenario, "\n")
cat("n =", row$n, "\n")
cat("m =", row$m, "\n")
cat("ZI =", row$zi_level, "\n")
cat("zeta_j ~ U(", row$zeta_low, ", ", row$zeta_high, ")\n", sep = "")
cat("alpha_i =", row$alpha_i, "\n")
cat("phi = 1\n")
cat("replications =", row$rep_k, "\n")
cat("NO dataset filtering\n")
cat("label =", row$label, "\n")

output_dir <- "results_ZINB"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outfile_rds <- file.path(output_dir, paste0(row$label, ".rds"))
outfile_csv <- file.path(output_dir, paste0(row$label, "_full.csv"))

if(file.exists(outfile_rds)){
  cat("Skipping existing file:", outfile_rds, "\n")
  quit(save = "no")
}

res <- try(
  run_simulation_ZINB(
    rep_k = row$rep_k,
    n = row$n,
    m = row$m,
    k = 2,
    zi_level = row$zi_level,
    depth_level = row$depth_level,
    phi = 1,
    init_method = "DS",
    seed_base = 4046000 + scenario_id * 10000
  ),
  silent = TRUE
)

if(inherits(res, "try-error")){
  cat("ZINB scenario failed.\n")
  failed_file <- file.path(output_dir, paste0("FAILED_", row$label, ".rds"))
  saveRDS(list(scenario = as.list(row), error = as.character(res)), failed_file)
  stop(res)
}

## Metadata only; no replicate is screened or rejected.
res$design_zeta_low <- row$zeta_low
res$design_zeta_high <- row$zeta_high
res$design_alpha_i <- row$alpha_i
res$design_phi <- 1

saveRDS(res, outfile_rds)
write.csv(res, outfile_csv, row.names = FALSE)

cat("\nSaved:\n")
cat(outfile_rds, "\n")
cat(outfile_csv, "\n")
cat("Finished ZINB scenario:", row$label, "\n")

}

# ==================== Command-line entry point ====================
# Running with no arguments is safe when the file is sourced in RStudio.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L || !args[[1L]] %in% c("ZIDM", "ZIMN", "ZINB")) {
    stop("Specify a data generator: Rscript simulation_main.R ZIDM|ZIMN|ZINB")
  }
  switch(args[[1L]],
         ZIDM = run_scenario_ZIDM(),
         ZIMN = run_scenario_ZIMN(),
         ZINB = run_scenario_ZINB())
}
