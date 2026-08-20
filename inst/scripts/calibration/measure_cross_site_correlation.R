# measure_cross_site_correlation.R
#
# Is site-to-site correlation the structural difference between reditR's
# simulator and real editing data?
#
# WHY. The parametric sweep shows the GLMM holding 2-4% FPR up to sigma = 1,
# yet permutation on real data breaks it. Per-site overdispersion, singular
# fits and shared sample effects have each been tested and excluded. The
# remaining candidate is that simulate_editing_data() draws every site
# independently -- rnorm(), rnbinom() and rbinom() are all called inside the
# per-site lapply -- whereas in real data a library that runs hot does so
# across many sites at once. BH assumes independence (or positive regression
# dependence) across tested hypotheses, so this bears directly on whether
# FDR control is licensed.
#
# MEASURE. For each site, take the vector of per-sample editing ratios and
# subtract the site mean, leaving each sample's deviation. Correlate those
# deviation vectors between pairs of sites. Independent sites -> ~0.
#
# SEPARATING TWO CAUSES. filter_editing_sites() keeps sites clustered within
# 50 bp, and neighbouring sites literally share reads, so some positive
# correlation is expected mechanically rather than as a sample effect. This
# script therefore reports correlation twice: over all site pairs, and over
# DIFFERENT-CHROMOSOME pairs only. Read-sharing cannot act across
# chromosomes, so the second number isolates a genome-wide sample effect.
#
# Only sites observed in every sample are used, so the deviation vectors are
# the same length and comparable. Each correlation rests on n_samples points
# and is individually noisy; the estimate is the mean over many pairs.
#
# Output: cross_site_correlation.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })
setwd("/rds/general/user/sj1825/home/msc_prj")
set.seed(20260817)
NMAX <- as.integer(Sys.getenv("NMAX", "400"))   # sites per dataset, caps the O(n^2) cor()

xsite <- function(dt, label) {
  s <- unique(dt[, .(site, sample, r = edit_ratio)], by = c("site", "sample"))
  full <- s[, .N, by = site][N == max(N), site]          # observed in every sample
  keep <- head(full, NMAX)
  m <- dcast(s[site %in% keep], site ~ sample, value.var = "r")
  M <- as.matrix(m[, -1]); rownames(M) <- m$site
  M <- M - rowMeans(M)                                    # deviation from site mean
  ok <- apply(M, 1, function(z) sd(z) > 0)
  M <- M[ok, , drop = FALSE]
  if (nrow(M) < 10) return(data.table(dataset = label, n_sites = nrow(M),
                                      r_all = NA_real_, r_diff_chr = NA_real_, n_samples = ncol(M)))
  cc <- suppressWarnings(cor(t(M)))
  ut <- upper.tri(cc)
  chr <- sub(":.*$", "", rownames(M))
  diffchr <- outer(chr, chr, "!=") & ut
  data.table(dataset = label, n_sites = nrow(M), n_samples = ncol(M),
             r_all = mean(cc[ut], na.rm = TRUE),
             r_diff_chr = if (any(diffchr)) mean(cc[diffchr], na.rm = TRUE) else NA_real_)
}

res <- list()
for (sg in c(0, 0.25, 1.0)) {
  sim <- simulate_editing_data(n_null = NMAX, n_effects = c(`0.10` = 0),
                               n_per_condition = 6L, sample_re_sd = sg, seed = 7L)
  d <- as.data.table(sim$editing)
  # simulated site ids carry no chromosome, so only r_all is defined
  res[[length(res) + 1]] <- xsite(d, sprintf("simulated sigma=%.2f", sg))
}
H <- "/rds/general/user/sj1825/home"
res[[length(res) + 1]] <- xsite(fread(file.path(H, "diabetes_output/filtered_sites_clustered_t.txt")),
                                "REAL diabetes (6v6)")
res[[length(res) + 1]] <- xsite(fread(file.path(H, "endothelial2/output/filtered_sites_clustered.txt")),
                                "REAL endothelial2 (9 samples)")

out <- rbindlist(res)
fwrite(out, "cross_site_correlation.txt", sep = "\t")
print(out)
cat("\nSimulated sites are independent by construction; any departure from 0 is sampling noise.\n")
cat("r_diff_chr isolates a genome-wide sample effect from 50 bp read-sharing.\n")
