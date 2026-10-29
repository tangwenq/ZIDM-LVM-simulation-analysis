# ==============================================================================
# Nugent (Ravel) data: preprocessing, model fitting, and result saving
#
# Models: gZIDM, gZIMN, gZINB, and ZIPPCA
# No covariates: X = NULL and V = NULL
# Fit each model once with k = 2
# ==============================================================================


# ==============================================================================
# 0. Configuration
# ==============================================================================

USE_FILTER <- FALSE          # TRUE: apply filtering; FALSE: use the full dataset
MIN_PREVALENCE <- 0.05       # Present in at least 5% of samples
K <- 2L
SEED <- 2026L

COMMON_FILE <- "common_realdata_methods.R"

options(timeout = 600)

if (!file.exists(COMMON_FILE)) {
  stop("Common script not found: ", COMMON_FILE,
       "\nPlease set the correct working directory.")
}

if (!file.exists("zimnVA.cpp")) {
  stop("zimnVA.cpp is missing from the working directory.")
}

# Load shared model functions and compile/load TMB as configured in the common script
source(COMMON_FILE)

stopifnot(
  exists("fit_and_save", mode = "function")
)


# ==============================================================================
# 1. Output directories
#
# Store filtered and full-data results separately.
# ==============================================================================

data_dir <- file.path("result", "Ravel", "source_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

analysis_tag <- if (USE_FILTER) {
  paste0(
    "filtered_prev",
    format(MIN_PREVALENCE, scientific = FALSE, trim = TRUE)
  )
} else {
  "full"
}

out_dir <- file.path(
  "result", "Ravel",
  paste0("models_", analysis_tag, "_XNULL")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Download data from GitHub
#
# Reuse existing local files.
# ==============================================================================

github_root <- paste0(
  "https://raw.githubusercontent.com/",
  "knights-lab/MLRepo/master/datasets/ravel/"
)

source_urls <- c(
  counts = paste0(github_root, "refseq/otutable.txt"),
  metadata = paste0(github_root, "mapping-orig.txt")
)

download_if_needed <- function(url, filename) {
  
  destination <- file.path(data_dir, filename)
  
  if (!file.exists(destination) ||
      file.info(destination)$size == 0) {
    
    tmp <- paste0(destination, ".download")
    
    status <- utils::download.file(
      url = url,
      destfile = tmp,
      mode = "wb",
      method = "libcurl"
    )
    
    if (status != 0L ||
        !file.exists(tmp) ||
        file.info(tmp)$size == 0) {
      stop("Download failed: ", url)
    }
    
    if (!file.copy(tmp, destination, overwrite = TRUE)) {
      stop("Unable to save downloaded file: ", destination)
    }
    
    unlink(tmp)
  }
  
  destination
}

otu_file <- download_if_needed(
  source_urls[["counts"]],
  "otutable.txt"
)

metadata_file <- download_if_needed(
  source_urls[["metadata"]],
  "mapping-orig.txt"
)


# ==============================================================================
# 3. Read data
#
# comment.char = "" preserves headers beginning with #
# check.names = FALSE preserves original sample IDs
# ==============================================================================

otu_table <- utils::read.delim(
  otu_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

metadata <- utils::read.delim(
  metadata_file,
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

required_fields <- c(
  "#SampleID", "Subject.ID", "Ethnic_Group", "pH",
  "Nugent_score", "Nugent_score_category",
  "Community_group", "Total_number_reads"
)

missing_fields <- setdiff(required_fields, names(metadata))

if (length(missing_fields) > 0L) {
  stop(
    "Missing metadata fields: ",
    paste(missing_fields, collapse = ", ")
  )
}

names(metadata)[names(metadata) == "#SampleID"] <- "SampleID"

# Transpose the OTU-by-sample table to sample-by-OTU format
Y_raw <- t(as.matrix(otu_table))
storage.mode(Y_raw) <- "double"

if (anyDuplicated(rownames(Y_raw)) ||
    anyDuplicated(colnames(Y_raw)) ||
    anyDuplicated(metadata$SampleID)) {
  stop("Duplicate sample or OTU IDs detected.")
}

if (!all(is.finite(Y_raw)) ||
    any(Y_raw < 0) ||
    any(Y_raw != floor(Y_raw))) {
  stop("Counts must be finite nonnegative integers.")
}

if (any(rowSums(Y_raw) == 0)) {
  stop("The original count matrix contains all-zero samples.")
}

if (!setequal(rownames(Y_raw), metadata$SampleID)) {
  
  print(list(
    Counts_without_metadata = setdiff(
      rownames(Y_raw), metadata$SampleID
    ),
    Metadata_without_counts = setdiff(
      metadata$SampleID, rownames(Y_raw)
    )
  ))
  
  stop("The count matrix and metadata have different sample sets.")
}


# ==============================================================================
# 4. Align metadata and prepare variables
# ==============================================================================

metadata <- metadata[
  match(rownames(Y_raw), metadata$SampleID),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$SampleID

stopifnot(
  identical(rownames(Y_raw), metadata$SampleID)
)

for (field in c("pH", "Nugent_score", "Total_number_reads")) {
  metadata[[field]] <- as.numeric(metadata[[field]])
}

metadata$Nugent_score_category <- factor(
  metadata$Nugent_score_category,
  levels = c("low", "intermediate", "high")
)

metadata$Community_group <- factor(
  metadata$Community_group,
  levels = c("I", "II", "III", "IV", "V")
)

metadata$Ethnic_Group <- factor(
  metadata$Ethnic_Group,
  levels = c("Asian", "Black", "Hispanic", "White")
)

if (anyNA(metadata)) {
  stop("Metadata contain missing values or unrecognized categories; please check.")
}


# ==============================================================================
# 5. Optional feature filtering
#
# Filter by prevalence only, without using Nugent scores, pH, or groups.
# ==============================================================================

stopifnot(
  length(USE_FILTER) == 1L,
  is.logical(USE_FILTER),
  !is.na(USE_FILTER),
  MIN_PREVALENCE > 0,
  MIN_PREVALENCE <= 1
)

otu_prevalence <- colSums(Y_raw > 0)
minimum_samples <- ceiling(nrow(Y_raw) * MIN_PREVALENCE)

keep_otu <- if (USE_FILTER) {
  otu_prevalence >= minimum_samples
} else {
  rep(TRUE, ncol(Y_raw))
}

Y <- Y_raw[, keep_otu, drop = FALSE]

if (ncol(Y) <= K || any(rowSums(Y) == 0)) {
  stop("Filtering left too few features or produced all-zero samples; adjust the filtering settings.")
}

# Calculate library sizes from the count matrix actually used for fitting
metadata$Count_total_raw <- rowSums(Y_raw)
metadata$Count_total_used <- rowSums(Y)

feature_summary <- data.frame(
  OTU = colnames(Y_raw),
  Total_count = colSums(Y_raw),
  Samples_present = otu_prevalence,
  Prevalence = otu_prevalence / nrow(Y_raw),
  Retained = keep_otu
)

data_summary <- data.frame(
  Version = c("Raw", "Used_for_models"),
  N_samples = c(nrow(Y_raw), nrow(Y)),
  N_features = c(ncol(Y_raw), ncol(Y)),
  Zero_proportion = c(mean(Y_raw == 0), mean(Y == 0)),
  Minimum_depth = c(min(rowSums(Y_raw)), min(rowSums(Y_raw[, keep_otu, drop = FALSE]))),
  Maximum_depth = c(max(rowSums(Y_raw)), max(rowSums(Y)))
)

print(data_summary)
print(table(metadata$Nugent_score_category))

cat(
  "\nFiltering enabled: ", USE_FILTER, "\n",
  "Fitting matrix: ", nrow(Y), " x ", ncol(Y), "\n",
  "Zero proportion: ", round(mean(Y == 0), 4), "\n",
  "Output directory: ", out_dir, "\n",
  sep = ""
)

utils::write.csv(
  data_summary,
  file.path(out_dir, "data_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  feature_summary,
  file.path(out_dir, "feature_filter_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    Y_raw = Y_raw,
    Y = Y,
    metadata = metadata,
    settings = list(
      use_filter = USE_FILTER,
      min_prevalence = MIN_PREVALENCE,
      minimum_samples = minimum_samples,
      k = K,
      seed = SEED,
      X = NULL,
      V = NULL
    ),
    source_urls = source_urls
  ),
  file.path(out_dir, "model_input.rds")
)


# ==============================================================================
# 6. Fit four models
#
# Full and portable fits are saved after each model.
# Continue with the remaining models if one model fails.
# Use extra = c(0, 0) for gZIDM and gZIMN.
# ==============================================================================

methods <- c("gZIDM", "gZIMN", "gZINB", "ZIPPCA")

fits <- list()
status_rows <- list()

for (method in methods) {
  
  tag <- paste0("Ravel_", analysis_tag, "_XNULL_", method)
  started <- Sys.time()
  
  result <- tryCatch(
    fit_and_save(
      method = method,
      tag = tag,
      Y = Y,
      X = NULL,
      V = NULL,
      k = K,
      out_dir = out_dir,
      extra = c(0, 0),
      seed = SEED
    ),
    error = function(e) e
  )
  
  elapsed <- as.numeric(
    difftime(Sys.time(), started, units = "secs")
  )
  
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
      Convergence = if (is.null(convergence)) {
        NA_character_
      } else {
        paste(convergence, collapse = ";")
      },
      Error = ""
    )
  }
  
  runtime_table <- do.call(rbind, status_rows)
  rownames(runtime_table) <- NULL
  
  saveRDS(
    fits,
    file.path(out_dir, "all_fits_portable.rds")
  )
  
  utils::write.csv(
    runtime_table,
    file.path(out_dir, "fit_status.csv"),
    row.names = FALSE
  )
  
  gc()
}

print(runtime_table)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))


