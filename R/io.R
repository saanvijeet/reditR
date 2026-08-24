#' Read an editing-site read-count table
#'
#' Reads a tab-separated file of per-site, per-sample editing counts.
#' If the \code{edit_ratio} column is absent it is computed as
#' \code{edited / total}.
#'
#' @param path Path to a tab-separated read-count file. Required columns:
#'   \code{site}, \code{sample}, \code{edited}, \code{total}.
#'
#' @return A \code{data.table} with columns \code{site}, \code{sample},
#'   \code{edited}, \code{total}, \code{edit_ratio}.
#'
#' @examples
#' f <- system.file("extdata", "example_editing.txt", package = "reditR")
#' read_editing_table(f)
#'
#' @export
read_editing_table <- function(path) {
  # Check that the input file exists before attempting to read it.
  if (!file.exists(path))
    stop("File not found: ", path)

  dt <- fread(path)
  # Check that all columns required for downstream analysis are present.
  required <- c("site", "sample", "edited", "total")
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0)
    stop("Missing required columns: ", paste(missing, collapse = ", "))

  # Calculate editing ratio if it was not included in the input.
  if (!"edit_ratio" %in% names(dt))
    dt[, edit_ratio := edited / total]
  
  # Return only the columns required by the downstream workflow.
  dt[, .(site, sample, edited, total, edit_ratio)]
}


#' Read sample metadata
#'
#' Reads a tab-separated sample metadata file. Optionally converts the
#' \code{condition} column to a factor with a specified reference level first,
#' which makes the GLMM coefficient name in \code{\link{differential_editing}}
#' predictable.
#'
#' @param path Path to a tab-separated metadata file. Required columns:
#'   \code{sample}, \code{condition}.
#' @param reference_level Character string naming the reference condition level.
#'   If supplied, \code{condition} is converted to a factor with this level
#'   first. An error is raised if the value is not found in the data.
#'
#' @return A \code{data.table} with columns \code{sample}, \code{condition}.
#'
#' @examples
#' f <- system.file("extdata", "example_metadata.txt", package = "reditR")
#' read_metadata(f, reference_level = "control")
#'
#' @export
read_metadata <- function(path, reference_level = NULL) {
  if (!file.exists(path))
    stop("File not found: ", path)

  dt <- fread(path)

  # Check that the columns required to link samples to conditions are present.
  required <- c("sample", "condition")
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0)
    stop("Missing required columns: ", paste(missing, collapse = ", "))

  if (!is.null(reference_level)) {
    # Check the requested reference level is actually present, and list the
    # available levels if it is not.
    if (!reference_level %in% dt$condition)
      stop("reference_level '", reference_level,
           "' not found in condition column. Available levels: ",
           paste(unique(dt$condition), collapse = ", "))

    # Put the reference condition first in the factor levels. Coercing a
    # character column to a factor orders the levels alphabetically, and
    # under treatment contrasts the first level is the baseline -- so the
    # baseline would otherwise be whichever condition happens to sort first,
    # which flips the sign of the estimated effect and changes the name of
    # the fitted coefficient.
    #
    # This is for callers who fit their own models from this table.
    # differential_editing() and editing_difference() read the metadata file
    # directly rather than through read_metadata(), and differential_editing()
    # sets the contrast itself, so neither depends on the ordering set here.
    other_levels <- setdiff(unique(as.character(dt$condition)), reference_level)
    dt[, condition := factor(condition,
                             levels = c(reference_level, other_levels))]
  }
  # Return only the metadata required by the downstream workflow.
  dt[, .(sample, condition)]
}
