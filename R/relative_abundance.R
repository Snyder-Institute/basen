#' Calculate species-level relative abundance from Kraken2 outputs
#'
#' @param report_dir Directory containing *.report files. If NULL, kraken_dir
#'   will be used.
#' @param kraken_dir Directory containing *.kraken files. If NULL, report_dir
#'   will be used.
#' @param k_mer_length Positive integer specifying the k-mer length used by
#'   Kraken2. This value is subtracted from genome_size_bp when computing
#'   coverage proxies. Default is 35.
#' @param genome_stats_file Full path to genome stats file. The first column
#'   must be taxid and the second column genome_size_bp.
#' @param sample_ids Optional character vector of sample IDs (file basenames)
#'   to process. If NULL, all samples found in the directory are processed.
#' @param output_folder Optional directory to write per-sample result files.
#'   If provided, one tab-delimited file named <sample_id>.basen is written
#'   for each processed sample.
#' @param verbose Logical; print progress and sanity checks.
#'
#' @return A data.table with columns:
#'   sample_id, name, taxid, genome_size_bp, bases_assigned,
#'   reads_assigned, coverage_proxy, relative_abundance
#'
#' @export
#'
#' @importFrom dplyr filter mutate rowwise ungroup group_by summarise left_join arrange select
#' @importFrom data.table as.data.table fwrite
#'
#' @examples
#' ## Example 1: Process all samples and return results only
#' result_in_list <- kraken_relative_abundance(
#'         report_dir = report_dir,
#'         genome_stats_file = genome_stats,
#'         verbose = TRUE
#' )
#'
#' ## Example 2: Process a subset of samples and export results
#' result_in_list <- kraken_relative_abundance(
#'         report_dir = report_dir,
#'         kraken_dir = kraken_dir,
#'         genome_stats_file = genome_stats,
#'         sample_ids = c("SampleA", "SampleB"),
#'         output_folder = "relative_abundance_results",
#'         verbose = TRUE
#' )
kraken_relative_abundance <- function(report_dir = NULL,
                                      kraken_dir = NULL,
                                      k_mer_length = 35,
                                      genome_stats_file,
                                      sample_ids = NULL,
                                      output_folder = NULL,
                                      verbose = TRUE) {
        ## Directory handling
        if (is.null(report_dir) && is.null(kraken_dir)) {
                stop("At least one of report_dir or kraken_dir must be provided.")
        }

        if (is.null(report_dir)) report_dir <- kraken_dir
        if (is.null(kraken_dir)) kraken_dir <- report_dir

        if (!is.numeric(k_mer_length) ||
                length(k_mer_length) != 1 ||
                is.na(k_mer_length) ||
                k_mer_length <= 0 ||
                k_mer_length %% 1 != 0) {
                stop("k_mer_length must be a single positive integer.")
        }

        if (!dir.exists(report_dir)) {
                stop("Report directory does not exist: ", report_dir)
        }
        if (!dir.exists(kraken_dir)) {
                stop("Kraken directory does not exist: ", kraken_dir)
        }
        if (!file.exists(genome_stats_file)) {
                stop("Genome stats file does not exist: ", genome_stats_file)
        }
        if (!is.null(output_folder)) {
                if (!dir.exists(output_folder)) {
                        dir.create(output_folder, recursive = TRUE)
                        if (verbose) {
                                message("Created output folder: ", output_folder)
                        }
                }
        }

        ## Discover available files
        report_files <- list.files(
                report_dir,
                pattern = "\\.report$",
                full.names = TRUE
        )

        if (length(report_files) == 0) {
                stop("No .report files found in report directory.")
        }

        available_ids <- sub("\\.report$", "", basename(report_files))

        ## Subset samples if requested
        if (!is.null(sample_ids)) {
                missing <- setdiff(sample_ids, available_ids)
                if (length(missing) > 0) {
                        stop(
                                "The following sample_ids were not found as .report files: ",
                                paste(missing, collapse = ", ")
                        )
                }
                sample_ids <- intersect(sample_ids, available_ids)
        } else {
                sample_ids <- available_ids
        }

        if (verbose) {
                message("Processing ", length(sample_ids), " sample(s).")
        }

        ## Begin
        genome_stats <- read.table(
                genome_stats_file,
                header = FALSE,
                sep = "\t",
                stringsAsFactors = FALSE,
                comment.char = "",
                quote = ""
        )

        ## Sanity check: minimum required columns
        if (ncol(genome_stats) < 2) {
                stop(
                        "The genome stats file must contain at least two columns: ",
                        "column 1 = taxid, column 2 = genome_size_bp."
                )
        }

        ## Keep only required columns and standardize names
        genome_stats <- genome_stats[, 1:2]
        colnames(genome_stats) <- c("taxid", "genome_size_bp")

        ## Type coercion and checks
        genome_stats$taxid <- as.character(genome_stats$taxid)
        genome_stats$genome_size_bp <- as.numeric(genome_stats$genome_size_bp)

        if (anyNA(genome_stats$genome_size_bp)) {
                stop("Non-numeric genome_size_bp values detected in genome stats file.")
        }

        ## Adjust genome size by k-mer length
        genome_stats$genome_size_bp <- genome_stats$genome_size_bp - k_mer_length - 1

        if (any(genome_stats$genome_size_bp <= 1, na.rm = TRUE)) {
                stop(
                        "Adjusted genome_size_bp (genome_size_bp - k_mer_length - 1) ",
                        "must be positive for all taxa."
                )
        }

        if (verbose) {
                message(
                        "Loaded genome stats for ",
                        nrow(genome_stats),
                        " taxa (adjusted genome_size_bp = genome_size_bp - ",
                        k_mer_length,
                        ")."
                )
        }

        resultL <- lapply(sample_ids, function(sample_id) {
                if (verbose) message("Processing sample: ", sample_id)

                report_path <- file.path(report_dir, paste0(sample_id, ".report"))
                kraken_path <- file.path(kraken_dir, paste0(sample_id, ".kraken"))

                if (!file.exists(kraken_path)) {
                        warning("Missing .kraken file for ", sample_id, "; skipping...")
                        next
                }

                report <- read.table(
                        report_path,
                        sep = "\t",
                        header = FALSE,
                        quote = "",
                        comment.char = "",
                        stringsAsFactors = FALSE
                )
                colnames(report) <- c("percent", "reads_clade", "reads_direct", "rank_code", "taxid", "name")
                report$taxid <- as.character(report$taxid)

                kraken <- read.table(
                        kraken_path,
                        sep = "\t",
                        header = FALSE,
                        quote = "",
                        comment.char = "",
                        stringsAsFactors = FALSE
                )
                colnames(kraken) <- c("status", "read_id", "taxid", "read_length", "kmer_string")
                kraken$taxid <- as.character(kraken$taxid)

                species_report <- subset(report, rank_code == "S")

                if (nrow(species_report) == 0) {
                        stop("No species-level taxa found for sample ", sample_id)
                }

                if (verbose) {
                        message(" - Species-level TaxIDs: ", nrow(species_report))
                }

                if (verbose) {
                        message(" - Total reads: ", nrow(kraken))
                }

                kraken_c <- subset(kraken, status == "C")

                if (verbose) {
                        message(" - Classified reads: ", nrow(kraken_c), " (", round(nrow(kraken_c) / nrow(kraken) * 100, digits = 2), "%)")
                }

                species_kmers <- kraken_c |>
                        dplyr::filter(taxid %in% species_report$taxid) |>
                        dplyr::rowwise() |>
                        dplyr::mutate(assigned_kmers = extract_kmer_counts(kmer_string, taxid)) |>
                        dplyr::ungroup()

                if (verbose) {
                        message(" - Directly assigned reads: ", nrow(species_kmers), " (", round(nrow(species_kmers) / nrow(kraken) * 100, digits = 2), "%)")
                }

                base_counts <- species_kmers |>
                        dplyr::group_by(taxid) |>
                        dplyr::summarise(bases_assigned = sum(assigned_kmers), reads_assigned = dplyr::n(), .groups = "drop") |>
                        dplyr::left_join(species_report[, c("taxid", "name")], by = "taxid") |>
                        dplyr::mutate(name = trimws(name)) |>
                        dplyr::left_join(genome_stats[, c("taxid", "genome_size_bp")], by = "taxid")

                if (any(is.na(base_counts$genome_size_bp))) {
                        warning(
                                "Missing genome sizes for TaxIDs: ",
                                paste(base_counts$taxid[is.na(base_counts$genome_size_bp)],
                                        collapse = ", "
                                )
                        )
                }

                base_counts <- base_counts |>
                        dplyr::mutate(coverage_proxy = bases_assigned / genome_size_bp)

                total_coverage <- sum(base_counts$coverage_proxy, na.rm = TRUE)

                base_counts <- base_counts |>
                        dplyr::mutate(relative_abundance = coverage_proxy / total_coverage) |>
                        dplyr::arrange(desc(relative_abundance)) |>
                        dplyr::select(name, taxid, genome_size_bp, bases_assigned, reads_assigned, coverage_proxy, relative_abundance)

                result_dt <- data.table::as.data.table(base_counts)
                result_dt[, sample_id := sample_id]

                data.table::setcolorder(
                        result_dt,
                        c("sample_id",
                        "name",
                        "taxid",
                        "genome_size_bp",
                        "bases_assigned",
                        "reads_assigned",
                        "coverage_proxy",
                        "relative_abundance"
                        )
                )

                ## Optional export
                if (!is.null(output_folder)) {
                        out_file <- file.path(output_folder, paste0(sample_id, ".basen"))
                        data.table::fwrite(result_dt, file = out_file, sep = "\t")
                        if (verbose) {
                                message(" - Wrote output file: ", out_file)
                        }
                }

                return(result_dt)
        })

        combined <- data.table::rbindlist(
                resultL,
                use.names = TRUE,
                fill = TRUE
        )

        return(combined)
}

#' @keywords internal
extract_kmer_counts <- function(kmer_string, target_taxid = NULL) {

        pairs <- unlist(strsplit(kmer_string, " ", fixed = TRUE))

        taxids <- sub(":.*", "", pairs)
        counts <- as.numeric(sub(".*:", "", pairs))

        result <- NULL

        ## Case 1: target_taxid is provided
        if (!is.null(target_taxid)) {
                result <- sum(
                        counts[taxids == as.character(target_taxid)],
                        na.rm = TRUE
                )
        } else {
                ## Case 2: aggregate all taxids
                agg <- tapply(
                        counts,
                        taxids,
                        sum,
                        na.rm = TRUE
                )

                ## Ensure deterministic ordering
                ord <- order(as.numeric(names(agg)))
                result <- agg[ord]
        }

        return(result)
}

#' Collect BASEN per-sample metrics from *.basen outputs
#'
#' This function reads multiple per-sample BASEN result files and extracts the
#' requested metric(s) into a single long-format table.
#'
#' Users can choose to collect relative abundance, a coverage proxy, or both.
#' Genome metadata columns (e.g., name, genome_size) are intentionally not returned.
#'
#' @param output_folder Directory containing per-sample result files, or NULL
#'   if \code{files} is provided.
#' @param files Optional character vector of full paths to result files.
#'   If provided, \code{output_folder} is ignored.
#' @param pattern File name pattern used to identify result files.
#'   Defaults to \code{"\\.basen$"}.
#' @param metric Character; which metric(s) to collect. One of:
#'   \code{"relative_abundance"}, \code{"coverage_proxy"}, or \code{"both"}.
#' @param verbose Logical; print progress and sanity checks.
#'
#' @return A \code{data.table} with columns:
#'   \itemize{
#'     \item \code{sample_id}
#'     \item \code{taxid}
#'     \item \code{relative_abundance} (if selected)
#'     \item \code{coverage_proxy} (if selected)
#'   }
#'
#' @export
#'
#' @importFrom data.table fread rbindlist as.data.table
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' ## Example 1: Collect both metrics from a directory
#' output_folder <- "path/to/basen_outputs"
#' res_all <- collect_basen_metrics(
#'         output_folder = output_folder,
#'         metric = "both",
#'         verbose = TRUE
#' )
#'
#' ## Example 2: Collect coverage proxy from a specific set of files
#' files <- c(
#'         "path/to/SampleA.basen",
#'         "path/to/SampleB.basen"
#' )
#' res_subset <- collect_basen_metrics(
#'         files = files,
#'         metric = "coverage_proxy",
#'         verbose = FALSE
#' )
collect_basen_metrics <- function(output_folder = NULL,
                                  files = NULL,
                                  pattern = "\\.basen$",
                                  metric = c("relative_abundance", "coverage_proxy", "both"),
                                  verbose = TRUE) {
        ## Validate metric choice early (also gives a clean error message)
        metric <- match.arg(metric)

        ## Input handling
        if (is.null(output_folder) && is.null(files)) {
                stop("Either output_folder or files must be provided.")
        }

        if (!is.null(files)) {
                missing <- files[!file.exists(files)]
                if (length(missing) > 0) {
                        stop(
                                "The following files do not exist: ",
                                paste(missing, collapse = ", ")
                        )
                }
                result_files <- files
        } else {
                if (!dir.exists(output_folder)) {
                        stop("Results directory does not exist: ", output_folder)
                }

                result_files <- list.files(
                        output_folder,
                        pattern = pattern,
                        full.names = TRUE
                )

                if (length(result_files) == 0) {
                        stop("No result files found in directory.")
                }
        }

        if (verbose) {
                message("Reading ", length(result_files), " result file(s).")
                message("Selected metric: ", metric)
        }

        ## Define which columns to collect (taxid is the minimal stable identifier)
        metric_cols <- switch(
                metric,
                relative_abundance = "relative_abundance",
                coverage_proxy     = "coverage_proxy",
                both               = c("relative_abundance", "coverage_proxy")
        )
        required_cols <- c("taxid", metric_cols)

        ## Read and collect
        res <- lapply(result_files, function(f) {
                dt <- fread(f, sep = "\t")

                missing_cols <- setdiff(required_cols, colnames(dt))
                if (length(missing_cols) > 0) {
                        stop(
                                "File ", basename(f),
                                " is missing required columns: ",
                                paste(missing_cols, collapse = ", ")
                        )
                }

                sample_id <- tools::file_path_sans_ext(basename(f))

                ## Keep only the requested columns (no name/taxonomy strings/genome_size here)
                out <- dt[, ..required_cols]
                out[, sample_id := sample_id]

                ## Reorder columns in a consistent way
                out <- out[, c("sample_id", required_cols), with = FALSE]

                return(out)
        })

        combined <- rbindlist(res, use.names = TRUE)

        if (verbose) {
                message(
                        "Collected ",
                        nrow(combined),
                        " taxon-level entries across ",
                        length(unique(combined$sample_id)),
                        " samples."
                )
        }

        return(as.data.table(combined))
}

#' Calculate size factors from a coverage proxy matrix
#'
#' This function calculates per-sample size factors from a coverage proxy matrix
#' (rows = taxids, columns = samples, values = \code{coverage_proxy}).
#'
#' Important: size factors should be computed from \code{coverage_proxy} (a depth proxy),
#' not from \code{relative_abundance}. Relative abundance is compositional and does not
#' retain the total-depth information that size-factor scaling is designed to capture.
#'
#' The implementation follows the standard "median ratio" family of methods:
#' 1) build a per-taxon reference across samples (geometric mean of coverage_proxy),
#' 2) compute per-sample ratios (coverage_proxy / reference),
#' 3) summarize ratios into a single size factor per sample via median or geomean.
#'
#' @param x A numeric matrix (or coercible to matrix) with taxids as rownames and
#'   sample IDs as colnames. Cell values must be \code{coverage_proxy}.
#' @param method Character; one of \code{"median"} or \code{"geomean"}.
#' @param verbose Logical; print progress and sanity checks.
#'
#' @return A numeric vector of size factors named by sample (column) names.
#'
#' @export
#'
#' @import stats
#' @importFrom stats median
#'
#' @examples
#' ## x is a taxid-by-sample matrix of coverage_proxy values
#' ## rownames(x) = taxids, colnames(x) = sample IDs
#' ##
#' ## sf <- calc_size_factors(x, method = "median")
#' ## sf
#'
#' ## Optional: normalize coverage_proxy by size factors
#' ## x_norm <- sweep(x, 2, sf, FUN = "/")
calc_size_factors <- function(x,
                              method = c("median", "geomean"),
                              verbose = TRUE) {
        method <- match.arg(method)

        if (isTRUE(verbose)) {
                message("Note: calc_size_factors expects coverage_proxy input (not relative_abundance).")
        }

        ## Coerce to numeric matrix
        if (is.data.frame(x)) {
                x <- as.matrix(x)
        }
        if (!is.matrix(x)) {
                stop("Input x must be a matrix (or coercible to a matrix).")
        }

        ## Ensure numeric storage mode
        if (!is.numeric(x)) {
                suppressWarnings(storage.mode(x) <- "numeric")
        }
        if (!is.numeric(x)) {
                stop("Input x must be numeric (coverage_proxy values).")
        }

        ## Basic dimension checks
        if (nrow(x) < 2) {
                stop("x must have at least 2 rows (taxids) to compute a stable reference.")
        }
        if (ncol(x) < 2) {
                stop("x must have at least 2 columns (samples) to compute size factors.")
        }

        ## Require names (taxids and sample IDs)
        if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
                stop("x must have non-empty rownames representing taxids.")
        }
        if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
                stop("x must have non-empty colnames representing sample IDs.")
        }

        ## Sanity check for non-finite values
        if (any(!is.finite(x), na.rm = TRUE)) {
                stop("x contains non-finite values (NA/NaN/Inf). Please fix or remove them before computing size factors.")
        }

        ## Filter taxa that cannot contribute to geometric mean reference:
        ## geometric mean requires positive values across the contributing samples.
        ## A common convention: drop taxa with any zero across samples for the reference.
        taxa_ok <- apply(x, 1, function(v) all(v > 0))
        if (sum(taxa_ok) == 0) {
                stop(
                        "No taxa have strictly positive coverage_proxy across all samples.\n",
                        "Size-factor estimation via geometric-mean reference is undefined in this case.\n",
                        "Consider filtering to prevalent taxa or adding a small pseudocount upstream."
                )
        }
        x_use <- x[taxa_ok, , drop = FALSE]

        if (isTRUE(verbose)) {
                message("Using ", nrow(x_use), " / ", nrow(x), " taxa with strictly positive values across all samples.")
        }

        ## Per-taxon reference = geometric mean across samples (log-space)
        ## ref_i = exp(mean_j log(x_ij))
        ref <- exp(rowMeans(log(x_use)))

        ## Compute ratios per sample for each taxon: x_ij / ref_i
        ratios <- sweep(x_use, 1, ref, FUN = "/")

        ## Summarize ratios into a size factor per sample
        if (method == "median") {
                sf <- apply(ratios, 2, function(v) stats::median(v, na.rm = TRUE))
        } else {
                ## geomean of ratios per sample: exp(mean log(ratio))
                sf <- apply(ratios, 2, function(v) exp(mean(log(v), na.rm = TRUE)))
        }

        ## Final sanity checks
        if (any(!is.finite(sf)) || any(sf <= 0)) {
                stop(
                        "Computed invalid size factors (non-finite or non-positive). ",
                        "Check that coverage_proxy values are strictly positive for the taxa used."
                )
        }

        ## Name the vector by sample IDs
        sf <- setNames(as.numeric(sf), colnames(x_use))

        if (isTRUE(verbose)) {
                message("Computed size factors for ", length(sf), " sample(s) using method: ", method, ".")
        }

        return(sf)
}