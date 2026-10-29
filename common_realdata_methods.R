# ==============================================================================
# common_realdata_methods.R
#
# COMMON MODEL-METHOD LAYER FOR REAL-DATA ANALYSES
#
# Shared model fitting and result-saving utilities.
# ==============================================================================

options(stringsAsFactors = FALSE, timeout = 600)

user_lib <- path.expand("~/R/library_452")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(TMB)
  library(gllvm)
  library(MASS)
  library(compositions)
  library(robCompositions)
  library(vegan)
  #library(ZIM)
  library(ZIPPCAlnm)
})

# The TMB model source zimnVA.cpp must be in the working directory.
TMB::compile(
  "zimnVA.cpp",
  flags = "-O3 -march=native -mtune=native",
  openmp = TRUE,
  framework = "TMBad"
)

dyn.load(dynlib("zimnVA"))
openmp(max = TRUE, autopar = TRUE, DLL = "zimnVA")

# ==============================================================================
# Model fitting functions
# ==============================================================================


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
  # X handling: always include an intercept
  # =========================================================
  if (is.null(X)) {
    
    # Include an intercept even when no covariates are supplied
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
    
    # Important
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
      
      # Include an intercept
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
  ## Dunn-Smyth residuals
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
      ## Dunn-Smyth randomization
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

# ==============================================================================
# Result extraction and saving helpers
# ==============================================================================


safe_fit <- function(expr) {
  tryCatch(
    list(ok = TRUE, fit = eval.parent(substitute(expr)), error = NULL),
    error = function(e) list(ok = FALSE, fit = NULL, error = conditionMessage(e))
  )
}

safe_save_rds <- function(object, path, compress = "gzip") {
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

make_portable <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  
  if (method %in% c("gZIDM", "gZIMN")) {
    keep <- c(
      "beta0_hat", "beta_hat", "b_hat", "u_hat", "gamma_hat",
      "theta_hat", "Au_hat", "lgt_zeta_hat", "zeta_hat",
      "lg_Si_hat", "Si_hat", "rho_hat", "rhoZI_hat", "eta_hat",
      "cQ_hat", "ordination", "loadings", "convergence", "objective",
      "logl", "nll", "model", "extra", "init_method", "all_parameters"
    )
    return(fit[intersect(keep, names(fit))])
  }
  
  if (method == "gZINB") {
    keep <- c(
      "alpha_hat", "beta0_hat", "beta_hat", "u_hat", "gamma_hat",
      "zeta_hat", "dispersion_hat", "inv_dispersion_hat", "A_hat",
      "cQ_hat", "eta_hat", "eta_corrected_hat", "rho_hat",
      "theta_raw", "sigma_lv", "ordination", "loadings"
    )
    return(fit[intersect(keep, names(fit))])
  }
  
  if (method == "ZIPPCA") {
    return(list(
      Q2 = fit$Q2,
      factor_scores = if (!is.null(fit$lvs$factor_scores)) fit$lvs$factor_scores else NULL,
      fit_names = names(fit)
    ))
  }
  
  if (method == "ZIDM-K") {
    # The full MCMC output is saved separately.
    # Keep the portable checkpoint light by excluding raw_fit.
    return(fit[setdiff(names(fit), "raw_fit")])
  }
  
  fit
}

normalize_rho <- function(rho) {
  rho <- as.matrix(rho)
  rs <- rowSums(rho)
  if (any(!is.finite(rs)) || any(rs <= 0)) stop("Invalid row sums in composition matrix.")
  rho / rs
}

extract_rho <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  rho <- switch(
    method,
    gZIDM = fit$rho_hat,
    gZIMN = fit$rho_hat,
    gZINB = fit$rho_hat,
    ZIPPCA = fit$Q2,
    `ZIDM-K` = fit$rho_hat,
    NULL
  )
  if (is.null(rho)) return(NULL)
  normalize_rho(rho)
}

extract_ord <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  U <- switch(
    method,
    gZIDM = fit$ordination,
    gZIMN = fit$ordination,
    gZINB = fit$ordination,
    ZIPPCA = if (!is.null(fit$lvs$factor_scores)) fit$lvs$factor_scores else NULL,
    `ZIDM-K` = NULL,
    NULL
  )
  if (is.null(U)) return(NULL)
  U <- as.matrix(U)
  if (ncol(U) < 2L || anyNA(U) || any(!is.finite(U))) return(NULL)
  scale(U)
}

extract_beta <- function(fit, method) {
  if (is.null(fit)) return(NULL)
  B <- switch(
    method,
    gZIDM = fit$beta_hat,
    gZIMN = fit$beta_hat,
    gZINB = fit$beta_hat,
    `ZIDM-K` = fit$beta_hat,
    NULL
  )
  if (is.null(B)) return(NULL)
  as.matrix(B)
}

# ==============================================================================
# Shared model fitting and saving interface
# ==============================================================================

save_model_bundle <- function(fit, method, tag, out_dir, runtime) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  full_file <- file.path(out_dir, paste0(tag, "_fit_full.rds"))
  portable_file <- file.path(out_dir, paste0(tag, "_fit_portable.rds"))
  
  saveRDS(fit, full_file, compress = "gzip")
  portable <- make_portable(fit, method)
  
  if (method == "ZIPPCA") {
    if (is.null(portable$rho_hat) && !is.null(portable$Q2)) {
      portable$rho_hat <- normalize_rho(portable$Q2)
    }
    if (is.null(portable$ordination) && !is.null(portable$factor_scores)) {
      portable$ordination <- as.matrix(portable$factor_scores)
    }
  }
  
  portable$runtime <- runtime
  portable$model <- method
  saveRDS(portable, portable_file, compress = "gzip")
  
  rho <- try(extract_rho(fit, method), silent = TRUE)
  if (!inherits(rho, "try-error") && !is.null(rho)) {
    saveRDS(rho, file.path(out_dir, paste0(tag, "_rho.rds")))
  }
  
  ord <- try(extract_ord(fit, method), silent = TRUE)
  if (!inherits(ord, "try-error") && !is.null(ord)) {
    saveRDS(ord, file.path(out_dir, paste0(tag, "_ordination.rds")))
  }
  
  beta <- try(extract_beta(fit, method), silent = TRUE)
  if (!inherits(beta, "try-error") && !is.null(beta)) {
    saveRDS(beta, file.path(out_dir, paste0(tag, "_beta.rds")))
  }
  
  list(
    fit = portable,
    runtime = runtime,
    full_file = full_file,
    portable_file = portable_file
  )
}

fit_and_save <- function(method, tag, Y, X = NULL, V = NULL, k = 2,
                         out_dir, extra = NULL,
                         zidmk_iterations = 10000,
                         zidmk_thin = 10,
                         zidmk_burnin = 0.5,
                         seed = 2026) {
  cat("\n============================================================\n")
  cat("Fitting:", tag, "\n")
  cat("============================================================\n")
  
  set.seed(seed)
  t0 <- Sys.time()
  
  fit <- switch(
    method,
    gZIDM = fit_TMB(
      Y = Y, X = X, k = k, model = "ZIDM",
      extra = extra, init_method = "DS"
    ),
    gZIMN = fit_TMB(
      Y = Y, X = X, k = k, model = "ZIMN",
      extra = extra, init_method = "DS"
    ),
    gZINB = fit_gllvm_ZINB(
      Y = Y, X = X, k = k, seed = seed
    ),
    ZIPPCA = ZIPPCAlnm::ZIPPCAlnm(
      X = Y,
      V = V,
      n.factors = k,
      rank = FALSE,
      trace = FALSE,
      maxit = 500,
      parallel = TRUE,
      sd.errors = FALSE,
      level = 0.95
    ),
    `ZIDM-K` = fit_ZIDM_paper(
      Y = Y,
      X = X,
      X_theta = NULL,
      iterations = zidmk_iterations,
      thin = zidmk_thin,
      burnin_prop = zidmk_burnin,
      seed = seed
    ),
    stop("Unknown method: ", method)
  )
  
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  save_model_bundle(fit, method, tag, out_dir, runtime)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
