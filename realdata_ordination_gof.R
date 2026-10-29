# Real-data ordination figures and goodness-of-fit evaluation
#
# Datasets: Nugent (Ravel) and microbiome (gllvm::microbialdata).
# Uses previously fitted models; no model fitting is performed.
# Models: ZIDM-LVM, ZIMN-LVM, ZINB-LVM, and ZIPPCA.
# CLR-PCA is included in the ordination figures only.
# Goodness-of-fit metrics: MARNE, gCOR, and mCOR.
# Information criteria (AIC, AICc, BIC) are not calculated.
#
# Run from the project root after nugent.R and microbiome.R.
# Update the paths below if your fitted results are elsewhere.

# Configuration ---------------------------------------------------------------

# Set this to models_filtered_prev0.05_XNULL for the filtered analysis.
nugent_dir <- file.path("result", "Ravel", "models_full_XNULL")
microbiome_dir <- file.path("results", "microbiome")
output_dir <- file.path("results", "realdata_ordination_gof")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (package in c("ggplot2", "patchwork", "gllvm")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package)
  }
}

method_names <- c(gZIDM = "ZIDM", gZIMN = "ZIMN",
                  gZINB = "ZINB", ZIPPCA = "ZIPPCA")
method_order <- c(unname(method_names), "CLR-PCA (+1)")

# Shared data and fit loading --------------------------------------------------

read_required <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  readRDS(path)
}

find_fit_file <- function(folder, tags, suffix = "_fit_full.rds") {
  candidates <- file.path(folder, paste0(tags, suffix))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("No saved fit found. Checked:\n", paste(candidates, collapse = "\n"))
  }
  found[[1L]]
}

load_dataset <- function(dataset) {
  if (dataset == "Nugent") {
    input <- read_required(file.path(nugent_dir, "model_input.rds"))
    Y <- as.matrix(input$Y)
    metadata <- input$metadata
    tag <- sub("^models_", "", basename(nugent_dir))
    tags <- list(
      gZIDM = c(paste0("Ravel_", tag, "_gZIDM"), "gZIDM_continued"),
      gZIMN = paste0("Ravel_", tag, "_gZIMN"),
      gZINB = paste0("Ravel_", tag, "_gZINB"),
      ZIPPCA = paste0("Ravel_", tag, "_ZIPPCA")
    )
    # A continued ZIDM fit is supported for previously completed analyses.
    zidm_continued <- file.path(nugent_dir, "gZIDM_continued",
                                "gZIDM_continued_fit_full.rds")
    fit_files <- vapply(names(tags), function(method) {
      if (method == "gZIDM" && file.exists(zidm_continued)) {
        return(zidm_continued)
      }
      find_fit_file(nugent_dir, tags[[method]])
    }, character(1))
  } else {
    input_path <- file.path(microbiome_dir, "model_input.rds")
    if (file.exists(input_path)) {
      input <- readRDS(input_path)
    } else {
      input <- read_required(file.path(microbiome_dir, "data_processed.rds"))
    }
    Y <- as.matrix(input$Y)
    metadata <- input$metadata
    tags <- list(
      gZIDM = c("A_XNULL_gZIDM", "A_XNULL_gZIDM_retain"),
      gZIMN = c("A_XNULL_gZIMN", "A_XNULL_gZIMN_retain"),
      gZINB = "A_XNULL_gZINB",
      ZIPPCA = "A_XNULL_ZIPPCA"
    )
    fit_files <- vapply(names(tags), function(method) {
      find_fit_file(microbiome_dir, tags[[method]])
    }, character(1))
  }

  if (!is.numeric(Y) || any(!is.finite(Y)) || any(Y < 0) ||
      any(Y != floor(Y)) || any(rowSums(Y) <= 0)) {
    stop(dataset, ": invalid count matrix.")
  }
  if (is.null(metadata) || nrow(metadata) != nrow(Y)) {
    stop(dataset, ": metadata is missing or has the wrong number of rows.")
  }
  list(Y = Y, metadata = metadata, fit_files = fit_files)
}

get_fit <- function(dataset, method, fit_file) {
  read_required(fit_file)
}

# Goodness-of-fit --------------------------------------------------------------
#
# MARNE: average across samples of MAE divided by each sample's count range.
# gCOR: Pearson correlation across all observed and predicted count cells.
# mCOR: mean taxon-wise Pearson correlation; NA if any taxon is undefined.
# A second mean across valid taxa is reported separately.

safe_cor <- function(x, y) {
  if (length(x) < 2L || any(!is.finite(x)) || any(!is.finite(y)) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = "pearson")
}

align_matrix <- function(M, Y, label) {
  M <- as.matrix(M)
  if (!identical(dim(M), dim(Y))) stop(label, ": dimensions differ from Y.")
  for (margin in 1:2) {
    from <- dimnames(M)[[margin]]
    to <- dimnames(Y)[[margin]]
    if (!is.null(from) && !is.null(to)) {
      if (anyDuplicated(from) || anyDuplicated(to) || !setequal(from, to)) {
        stop(label, ": sample or taxon names do not match.")
      }
      if (margin == 1L) M <- M[to, , drop = FALSE]
      if (margin == 2L) M <- M[, to, drop = FALSE]
    }
  }
  if (any(!is.finite(M))) stop(label, ": non-finite values.")
  dimnames(M) <- dimnames(Y)
  M
}

prepare_rho <- function(M, Y, label, normalize = FALSE) {
  rho <- align_matrix(M, Y, label)
  if (any(rho < 0)) stop(label, ": negative proportions.")
  totals <- rowSums(rho)
  if (any(totals <= 0)) stop(label, ": non-positive row sums.")
  if (!normalize && max(abs(totals - 1)) > 1e-6) {
    stop(label, ": rows do not sum to one.")
  }
  sweep(rho, 1, totals, "/")
}

predict_counts <- function(method, fit, Y, dataset) {
  if (method %in% c("gZIDM", "gZIMN")) {
    if (is.null(fit$rhoZI_hat)) {
      stop(method, ": missing rhoZI_hat.")
    }
    rho <- prepare_rho(fit$rhoZI_hat, Y, paste(method, "rhoZI_hat"))
    return(sweep(rho, 1, rowSums(Y), "*"))
  }
  if (method == "ZIPPCA") {
    rho_file <- if (dataset == "Nugent") NULL else
      file.path(microbiome_dir, "A_XNULL_ZIPPCA_rho.rds")
    if (!is.null(rho_file) && file.exists(rho_file)) {
      rho <- prepare_rho(readRDS(rho_file), Y, "ZIPPCA saved rho")
    } else {
      if (is.null(fit$Q2)) stop("ZIPPCA: missing Q2 and saved rho.")
      rho <- prepare_rho(fit$Q2, Y, "ZIPPCA Q2", normalize = TRUE)
    }
    return(sweep(rho, 1, rowSums(Y), "*"))
  }
  if (method == "gZINB") {
    model <- if (inherits(fit, "gllvm")) fit else fit$model
    if (!inherits(model, "gllvm")) stop("gZINB: missing gllvm model.")
    family <- unique(as.character(model$family))
    if (length(family) != 1L) stop("gZINB: inconsistent model families.")
    model$family <- family
    model_Y <- as.matrix(model$y)
    if (!isTRUE(all.equal(unname(align_matrix(model_Y, Y, "gZINB data")),
                          unname(Y), check.attributes = FALSE))) {
      stop("gZINB: fitted model data do not match Y.")
    }
    prediction <- as.matrix(stats::predict(model, type = "response", level = 1))
    if (!identical(dim(prediction), dim(model_Y))) {
      stop("gZINB: prediction dimensions do not match model data.")
    }
    # Predictions are generated in the fitted model's original sample order.
    dimnames(prediction) <- dimnames(model_Y)
    return(align_matrix(prediction, Y, "gZINB predictions"))
  }
  stop("Unsupported method: ", method)
}

calculate_gof <- function(Y, Y_hat, dataset, method) {
  if (!identical(dim(Y), dim(Y_hat)) ||
      any(!is.finite(Y_hat)) || any(Y_hat < 0)) {
    stop(dataset, " / ", method, ": invalid predicted counts.")
  }
  observed_range <- apply(Y, 1, function(x) diff(range(x)))
  valid_samples <- observed_range > 0
  mae <- rowMeans(abs(Y - Y_hat))
  marne <- rep(NA_real_, nrow(Y))
  marne[valid_samples] <- mae[valid_samples] / observed_range[valid_samples]
  correlations <- vapply(seq_len(ncol(Y)), function(j) {
    safe_cor(Y[, j], Y_hat[, j])
  }, numeric(1))
  valid_taxa <- is.finite(correlations)
  list(
    summary = data.frame(
      Dataset = dataset, Method = unname(method_names[[method]]),
      MARNE = if (any(valid_samples)) mean(marne[valid_samples]) else NA_real_,
      gCOR = safe_cor(as.vector(Y), as.vector(Y_hat)),
      mCOR = if (all(valid_taxa)) mean(correlations) else NA_real_,
      mCOR_valid_species = if (any(valid_taxa)) mean(correlations[valid_taxa]) else NA_real_,
      N_samples = nrow(Y), N_species = ncol(Y),
      MARNE_samples_used = sum(valid_samples),
      Correlation_species_used = sum(valid_taxa)
    ),
    species = data.frame(
      Dataset = dataset, Method = unname(method_names[[method]]),
      Taxon = if (is.null(colnames(Y))) paste0("Taxon_", seq_len(ncol(Y))) else colnames(Y),
      Pearson_correlation = correlations,
      Observed_mean = colMeans(Y), Fitted_mean = colMeans(Y_hat)
    ),
    samples = data.frame(
      Dataset = dataset, Method = unname(method_names[[method]]),
      SampleID = if (is.null(rownames(Y))) paste0("Sample_", seq_len(nrow(Y))) else rownames(Y),
      Observed_total = rowSums(Y), Fitted_total = rowSums(Y_hat),
      Observed_range = observed_range, MAE = mae, MARNE = marne
    )
  )
}

# Ordination ------------------------------------------------------------------
#
# The ordination figures use the same transformations as the source analysis:
# microbiome: SVD rotation without centering; CLR replaces zeros by one and
#             uses centered, scaled PCA.
# Nugent:     centered SVD rotation; CLR uses log(count + 1) and unscaled PCA.
# ZIDM-K has no ordination coordinates and is excluded from the figures.

ordination_scores <- function(fit, method, dataset, n) {
  U <- fit$ordination
  if (is.null(U)) U <- fit$factor_scores
  if (is.null(U) && method == "ZIPPCA" && !is.null(fit$lvs)) {
    U <- fit$lvs$factor_scores
  }
  if (is.null(U) && method == "gZINB") {
    model <- if (inherits(fit, "gllvm")) fit else fit$model
    if (!is.null(model)) {
      U <- model$lvs
      if (is.null(U)) U <- model$lv
    }
  }
  if (is.null(U)) stop(dataset, " / ", method, ": no saved latent scores.")
  U <- as.matrix(U)
  if (nrow(U) != n || ncol(U) < 2L || any(!is.finite(U))) {
    stop(dataset, " / ", method, ": invalid latent scores.")
  }
  U <- U[, 1:2, drop = FALSE]
  if (dataset == "Nugent") {
    U <- scale(U, center = TRUE, scale = FALSE)
  }
  U %*% svd(U)$v
}

clr_pca_scores <- function(Y, dataset) {
  if (dataset == "Nugent") {
    log_counts <- log(Y + 1)
    clr <- sweep(log_counts, 1, rowMeans(log_counts), "-")
    pca <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
    U <- pca$x[, 1:2, drop = FALSE]
    U <- scale(U, center = TRUE, scale = FALSE)
    return(U %*% svd(U)$v)
  }
  adjusted <- Y
  adjusted[adjusted == 0] <- 1
  log_counts <- log(adjusted)
  clr <- sweep(log_counts, 1, rowMeans(log_counts), "-")
  if (any(apply(clr, 2, stats::sd) == 0)) {
    stop("Microbiome CLR has constant taxa; scaled PCA is undefined.")
  }
  stats::prcomp(clr, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE]
}

plot_theme <- function() {
  ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13),
      aspect.ratio = 1
    )
}

make_ordination_figures <- function(dataset, input, fits, folder) {
  Y <- input$Y
  metadata <- input$metadata
  sample_ids <- rownames(Y)
  if (is.null(sample_ids)) sample_ids <- paste0("Sample_", seq_len(nrow(Y)))

  scores <- lapply(names(method_names), function(method) {
    U <- ordination_scores(fits[[method]], method, dataset, nrow(Y))
    data.frame(SampleID = sample_ids, LV1 = U[, 1], LV2 = U[, 2],
               Method = unname(method_names[[method]]))
  })
  U <- clr_pca_scores(Y, dataset)
  scores[[length(scores) + 1L]] <- data.frame(
    SampleID = sample_ids, LV1 = U[, 1], LV2 = U[, 2],
    Method = "CLR-PCA (+1)"
  )
  score_table <- do.call(rbind, scores)
  rownames(score_table) <- NULL

  if (dataset == "Nugent") {
    if (!"SampleID" %in% names(metadata) ||
        !setequal(metadata$SampleID, sample_ids)) {
      stop("Nugent metadata sample IDs do not match the count matrix.")
    }
    metadata <- metadata[match(sample_ids, metadata$SampleID), , drop = FALSE]
    groups <- list(
      Community_group = factor(metadata$Community_group,
                               levels = c("I", "II", "III", "IV", "V")),
      Nugent_category = factor(tolower(as.character(metadata$Nugent_score_category)),
                               levels = c("low", "intermediate", "high"))
    )
    palettes <- list(
      Community_group = c(I = "#D55E00", II = "#0072B2", III = "#009E73",
                          IV = "#CC79A7", V = "#E69F00"),
      Nugent_category = c(low = "#2166AC", intermediate = "#E69F00",
                          high = "#D73027")
    )
    for (variable in names(groups)) {
      if (anyNA(groups[[variable]])) stop("Missing Nugent group: ", variable)
      score_table[[variable]] <- rep(groups[[variable]], times = length(method_order))
      panels <- lapply(method_order, function(method) {
        dat <- score_table[score_table$Method == method, , drop = FALSE]
        ggplot2::ggplot(dat, ggplot2::aes(LV1, LV2, colour = .data[[variable]])) +
          ggplot2::geom_point(size = 1.6, alpha = 0.82) +
          ggplot2::scale_colour_manual(values = palettes[[variable]], drop = FALSE) +
          ggplot2::labs(title = method, x = "Ordination score 1",
                        y = "Ordination score 2", colour = variable) +
          plot_theme()
      })
      combined <- patchwork::wrap_plots(
        c(panels[1:3], panels[4:5], list(patchwork::plot_spacer())),
        ncol = 3, nrow = 2, guides = "collect"
      ) & ggplot2::theme(legend.position = "bottom")
      for (extension in c("pdf", "png")) {
        ggplot2::ggsave(file.path(folder, paste0("Nugent_", variable, ".", extension)),
                        combined, width = 13, height = 9, dpi = 300)
      }
    }
  } else {
    region <- factor(metadata$Region, levels = c("Aus", "Kil", "NyA"),
                     labels = c("Mayrhofen", "Kilpisjarvi", "Ny-Alesund"))
    if (anyNA(region)) stop("Missing or unknown microbiome regions.")
    score_table$Region <- rep(region, times = length(method_order))
    panels <- lapply(method_order, function(method) {
      dat <- score_table[score_table$Method == method, , drop = FALSE]
      ggplot2::ggplot(dat, ggplot2::aes(LV1, LV2, colour = Region, shape = Region)) +
        ggplot2::geom_point(size = 2.2, alpha = 0.95) +
        ggplot2::scale_colour_manual(values = c(
          Mayrhofen = "black", Kilpisjarvi = "blue", `Ny-Alesund` = "grey75"
        ), drop = FALSE) +
        ggplot2::scale_shape_manual(values = c(
          Mayrhofen = 16, Kilpisjarvi = 17, `Ny-Alesund` = 18
        ), drop = FALSE) +
        ggplot2::labs(title = method, x = "Ordination score 1",
                      y = "Ordination score 2", colour = NULL, shape = NULL) +
        plot_theme()
    })
    combined <- patchwork::wrap_plots(
      c(panels[1:3], panels[4:5], list(patchwork::plot_spacer())),
      ncol = 3, nrow = 2, guides = "collect"
    ) & ggplot2::theme(legend.position = "bottom")
    for (extension in c("pdf", "png")) {
      ggplot2::ggsave(file.path(folder, paste0("Microbiome_Region.", extension)),
                      combined, width = 13, height = 9, dpi = 300)
    }
  }
  utils::write.csv(score_table, file.path(folder, "ordination_scores_used.csv"),
                   row.names = FALSE)
  invisible(score_table)
}

# Run analyses and save outputs -----------------------------------------------
#
# A dataset failure does not prevent the other dataset from being processed.
# Each model's GOF is saved as soon as it has been calculated.

run_dataset <- function(dataset) {
  input <- load_dataset(dataset)
  folder <- file.path(output_dir, dataset)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  fits <- list()
  summaries <- list()
  species <- list()
  samples <- list()

  for (method in names(method_names)) {
    fit <- get_fit(dataset, method, input$fit_files[[method]])
    Y_hat <- predict_counts(method, fit, input$Y, dataset)
    gof <- calculate_gof(input$Y, Y_hat, dataset, method)
    summaries[[method]] <- gof$summary
    species[[method]] <- gof$species
    samples[[method]] <- gof$samples
    fits[[method]] <- fit
    utils::write.csv(do.call(rbind, summaries),
                     file.path(folder, "goodness_of_fit.csv"), row.names = FALSE)
  }

  utils::write.csv(do.call(rbind, species),
                   file.path(folder, "species_correlations.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, samples),
                   file.path(folder, "sample_fit_diagnostics.csv"), row.names = FALSE)

  make_ordination_figures(dataset, input, fits, folder)
  do.call(rbind, summaries)
}

results <- list()
errors <- list()
for (dataset in c("Nugent", "microbiome")) {
  message("Processing ", dataset, " ...")
  result <- tryCatch(run_dataset(dataset), error = function(e) e)
  if (inherits(result, "error")) {
    message("Failed: ", dataset, " - ", conditionMessage(result))
    errors[[dataset]] <- data.frame(Dataset = dataset,
                                    Error = conditionMessage(result))
  } else {
    results[[dataset]] <- result
  }
  gc(verbose = FALSE)
}

if (length(results)) {
  utils::write.csv(do.call(rbind, results),
                   file.path(output_dir, "all_goodness_of_fit.csv"),
                   row.names = FALSE)
}
if (length(errors)) {
  utils::write.csv(do.call(rbind, errors),
                   file.path(output_dir, "failed_datasets.csv"),
                   row.names = FALSE)
}
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Finished. Results: ", normalizePath(output_dir, mustWork = FALSE))
