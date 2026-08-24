#' Calculate per-site editing difference between conditions
#'
#' Computes the mean editing ratio per condition for each site and returns
#' the signed difference (case minus reference). Sites are ranked by absolute
#' editing difference. This is the R package equivalent of
#' \code{editing_difference.r} in the original pipeline.
#'
#' @param data_path Path to the filtered editing count file
#'   (\code{filtered_sites_clustered.txt} or similar). Required columns:
#'   \code{site}, \code{sample}, \code{edit_ratio}, \code{condition}.
#'   If \code{condition} is absent the metadata file is required via
#'   \code{meta_path}.
#' @param meta_path Path to the sample metadata file
#'   (\code{sample_metadata.txt}). Required columns: \code{sample},
#'   \code{condition}. Ignored if \code{condition} is already present in the
#'   data file.
#' @param out_path File path for the tab-delimited results table. Pass
#'   \code{NULL} (default) to skip writing.
#' @param reference_level Condition label for the reference group. Default
#'   \code{"control"}.
#' @param case_level Condition label for the case / disease group. Default
#'   \code{"diabetic"}.
#'
#' @return A \code{data.table} with columns \code{site},
#'   \code{<case_level>_mean}, \code{<reference_level>_mean},
#'   \code{editing_difference} (= case − reference), sorted by
#'   \code{abs(editing_difference)} descending.
#'
#' @examples
#' \dontrun{
#' eff <- editing_difference(
#'   data_path = "filtered_sites_clustered.txt",
#'   meta_path = "sample_metadata.txt"
#' )
#' head(eff)
#' }
#'
#' @export
editing_difference <- function(data_path,
                                meta_path      = NULL,
                                out_path       = NULL,
                                reference_level = "control",
                                case_level     = "diabetic") {

  # Load editing-site data
  data <- fread(data_path)

  # Calculate editing ratio if it has not already been provided.
  # This allows the function to work directly from edited and total read counts.
  if (!"edit_ratio" %in% names(data))
    data[, edit_ratio := edited / total]

  # Add condition information from the metadata if it is not already in the data.
  if (!"condition" %in% names(data)) {
    if (is.null(meta_path))
      stop("'condition' column not found in data file and no meta_path supplied.")
    meta <- fread(meta_path)
    data <- merge(data, meta, by = "sample")
  }

  # Use the condition names supplied by the user for the output columns.
  case_col <- paste0(case_level, "_mean")
  ref_col  <- paste0(reference_level, "_mean")

# Calculate the mean editing ratio for each condition at each site.
# Using the mean across samples gives a site-level estimate of the
# average editing level in each condition.
  editing_summary <- data[, .(
    case_mean = mean(edit_ratio[condition == case_level],      na.rm = TRUE),
    ref_mean  = mean(edit_ratio[condition == reference_level], na.rm = TRUE)
  ), by = site]

  setnames(editing_summary, c("case_mean", "ref_mean"), c(case_col, ref_col))

# Calculate the signed change in editing (case minus reference).
# Positive values indicate increased editing in the case condition;
# negative values indicate decreased editing.
  editing_summary$editing_difference <-
    editing_summary[[case_col]] - editing_summary[[ref_col]]

# Rank sites by the magnitude of their editing change.
# Absolute values ensure both increases and decreases are prioritised.
  editing_summary <- editing_summary[order(-abs(editing_difference))]

  if (!is.null(out_path))
    fwrite(editing_summary, out_path, sep = "\t")

  cat("Finished calculating editing differences\n")
  editing_summary
}
