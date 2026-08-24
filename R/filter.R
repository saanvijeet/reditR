#' Filter editing sites by coverage, read count, and cluster proximity
#'
#' Reads a merged editing count table, applies per-row QC filters, and then
#' identifies genomically clustered sites — sites within \code{cluster_window}
#' bp of a neighbouring site on the same chromosome. This is the R package
#' equivalent of \code{filtered_sites.R} in the original pipeline.
#'
#' The HPC/SPRINT helper script that produces the input file is shipped at
#' \code{system.file("scripts", "extracting_read_counts.sh", package = "reditR")}.
#'
#' @param data_path Path to the merged editing count file. Required columns:
#'   \code{site}, \code{sample}, \code{edited}, and \code{total}.
#' @param out_dir Directory for saving the filtered results. Set to
#'   \code{NULL} to return the results without writing files.
#' @param min_coverage Minimum total read coverage required per observation.
#'   Default is 10.
#' @param min_edited Minimum number of edited reads required per observation.
#'   Default is 2.
#' @param min_groups Minimum number of samples in which a site must be
#'   observed. Default is 2.
#' @param min_edit_ratio Minimum editing ratio required per observation.
#'   Default is 0, meaning no ratio filter is applied.
#' @param cluster_window Maximum distance in bp between neighbouring sites
#'   on the same chromosome for them to be considered clustered. Default is 50.
#'
#' @details
#' Filters are based on data quality and site coverage and are independent
#' of the experimental comparison. Condition-dependent filtering is not
#' applied because the appropriate criteria depend on the study design.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{all}}{All QC-passing sites.}
#'     \item{\code{clustered}}{Subset of sites within a genomic cluster.}
#'   }
#'
#' @examples
#' \dontrun{
#' res   <- filter_editing_sites("all_samples_editing.txt", out_dir = "results/")
#' sites <- res$clustered
#' }
#'
#' @export
filter_editing_sites <- function(data_path,
                                  out_dir        = NULL,
                                  min_coverage   = 10,
                                  min_edited     = 2,
                                  min_groups     = 2,
                                  min_edit_ratio = 0,
                                  cluster_window = 50) {

  # Load the editing site counts
  data <- read.table(data_path, header = TRUE)

  # Calculate editing ratio
  data <- data %>%
    mutate(edit_ratio = edited / total)

  # Remove observations with insufficient read coverage or edited reads.
  # These thresholds reduce the influence of poorly supported observations.
  data_filtered <- data %>%
    filter(total >= min_coverage, edited >= min_edited)

  # Apply an editing-ratio filter only when a threshold is requested.
  # This preserves all observations when min_edit_ratio = 0.
  if (min_edit_ratio > 0) {
    data_filtered <- data_filtered %>%
      filter(edit_ratio >= min_edit_ratio)
  }
  # Keep sites observed in enough independent samples to support
  # downstream differential-editing analysis.
  data_filtered <- data_filtered %>%
    group_by(site) %>%
    filter(n_distinct(sample) >= min_groups) %>%
    ungroup()

  # Split the site identifier into chromosome and genomic position
  # so that neighbouring sites can be identified.
  data_filtered <- data_filtered %>%
    separate(site, into = c("chr", "pos"), sep = ":")

  data_filtered$pos <- as.numeric(data_filtered$pos)

  # Identify sites that have another site within the clustering window.
  # Nearby sites are retained because genuine RNA-editing sites often
  # occur in genomic clusters.
  unique_sites <- data_filtered %>%
    distinct(chr, pos)

  clustered_positions <- unique_sites %>%
    group_by(chr) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(
      dist_prev = pos - lag(pos),
      dist_next = lead(pos) - pos
    ) %>%
    filter(dist_prev <= cluster_window | dist_next <= cluster_window) %>%
    ungroup()

  # Keep all sample observations belonging to clustered sites.
  clustered_sites <- data_filtered %>%
    semi_join(clustered_positions, by = c("chr", "pos"))

  # Reconstruct the site identifier and keep the columns needed downstream.
  data_filtered <- data_filtered %>%
    mutate(site = paste(chr, pos, sep = ":")) %>%
    select(site, sample, edited, total, edit_ratio, chr, pos)

  clustered_sites <- clustered_sites %>%
    mutate(site = paste(chr, pos, sep = ":")) %>%
    select(site, sample, edited, total, edit_ratio, chr, pos)

  # Save both the QC-filtered sites and the clustered subset if requested.
  if (!is.null(out_dir)) {
    write.table(data_filtered,
                file.path(out_dir, "filtered_sites_all.txt"),
                sep       = "\t",
                row.names = FALSE,
                quote     = FALSE)

    write.table(clustered_sites,
                file.path(out_dir, "filtered_sites_clustered.txt"),
                sep       = "\t",
                row.names = FALSE,
                quote     = FALSE)
  }

  invisible(list(all = data_filtered, clustered = clustered_sites))
}
