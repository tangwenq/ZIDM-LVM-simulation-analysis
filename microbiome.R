# Microbiome data: preprocessing, model fitting, and result saving
#
# All models are fitted without covariates.
# Only the soft settings (extra = c(0, 0)) are used for gZIDM and gZIMN.

source("common_realdata_methods.R")

# Configuration
k <- 2L
seed <- 2026L
zidmk_iterations <- 10000L
zidmk_thin <- 10L
zidmk_burnin <- 0.5

out_dir <- file.path("results", "microbiome")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the microbiome dataset and remove taxa with zero counts in all samples.
data("microbialdata", package = "gllvm")
Y_raw <- as.matrix(microbialdata$Y)
keep_taxa <- colSums(Y_raw != 0) > 0
Y <- Y_raw[, keep_taxa, drop = FALSE]

if (ncol(Y) <= k || any(rowSums(Y) == 0)) {
  stop("Insufficient taxa or samples with zero total counts.")
}

saveRDS(
  list(Y = Y, keep_taxa = keep_taxa, metadata = microbialdata$Xenv),
  file.path(out_dir, "model_input.rds")
)

cat("Microbiome: n =", nrow(Y), "p =", ncol(Y),
    "zero proportion =", mean(Y == 0), "\n")

# Fit the five models without covariates.
methods <- c("gZIDM", "gZIMN", "gZINB", "ZIDM-K", "ZIPPCA")
fits <- list()
status_rows <- list()

for (method in methods) {
  tag <- paste0("A_XNULL_", gsub("-", "", method, fixed = TRUE))
  started <- Sys.time()

  result <- tryCatch(
    fit_and_save(
      method = method,
      tag = tag,
      Y = Y,
      X = NULL,
      V = NULL,
      k = k,
      out_dir = out_dir,
      extra = if (method %in% c("gZIDM", "gZIMN")) c(0, 0) else NULL,
      zidmk_iterations = zidmk_iterations,
      zidmk_thin = zidmk_thin,
      zidmk_burnin = zidmk_burnin,
      seed = seed
    ),
    error = function(e) e
  )

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  if (inherits(result, "error")) {
    message(method, " failed: ", conditionMessage(result))
    status_rows[[method]] <- data.frame(
      Method = method,
      Returned_fit = FALSE,
      Runtime_sec = elapsed,
      Convergence = NA_character_,
      Error = conditionMessage(result)
    )
  } else {
    fits[[method]] <- result$fit
    convergence <- result$fit$convergence
    status_rows[[method]] <- data.frame(
      Method = method,
      Returned_fit = TRUE,
      Runtime_sec = result$runtime,
      Convergence = if (is.null(convergence)) NA_character_ else
        paste(convergence, collapse = ";"),
      Error = ""
    )
  }

  saveRDS(fits, file.path(out_dir, "all_fits_portable.rds"))
  status <- do.call(rbind, status_rows)
  rownames(status) <- NULL
  utils::write.csv(status, file.path(out_dir, "fit_status.csv"), row.names = FALSE)
  gc(verbose = FALSE)
}

saveRDS(
  list(
    k = k,
    seed = seed,
    covariates = NULL,
    zidmk_iterations = zidmk_iterations,
    zidmk_thin = zidmk_thin,
    zidmk_burnin = zidmk_burnin
  ),
  file.path(out_dir, "analysis_settings.rds")
)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(status)
