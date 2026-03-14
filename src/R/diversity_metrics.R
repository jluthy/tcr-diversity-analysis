# ============================================================================
# diversity_metrics.R
# Functions for computing TCR repertoire diversity indices
# Author: Joshua Luthy
# ============================================================================

library(tidyverse)

#' Compute Shannon Entropy
#' @param counts Integer vector of clone counts
#' @return Numeric Shannon entropy value
shannon_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log(p))
}

#' Compute Clonality Index
#' @param counts Integer vector of clone counts
#' @return Numeric clonality value (0 = perfectly even, 1 = monoclonal)
clonality_index <- function(counts) {
  H <- shannon_entropy(counts)
  N <- length(counts[counts > 0])
  if (N <= 1) return(0)
  1 - (H / log(N))
}

#' Compute D50 (number of clones comprising top 50% of reads)
#' @param counts Integer vector of clone counts
#' @return Integer D50 value
d50 <- function(counts) {
  sorted <- sort(counts, decreasing = TRUE)
  cum <- cumsum(sorted)
  total <- sum(sorted)
  which(cum >= total * 0.5)[1]
}

#' Compute Gini-Simpson Index
#' @param counts Integer vector of clone counts
#' @return Numeric Gini-Simpson value
gini_simpson <- function(counts) {
  p <- counts / sum(counts)
  1 - sum(p^2)
}

#' Compute all diversity metrics for a sample
#' @param df Data frame with clone_count column
#' @return Named list of diversity metrics
compute_diversity <- function(df) {
  counts <- df$clone_count
  list(
    n_clonotypes    = length(counts),
    total_reads     = sum(counts),
    shannon_entropy = shannon_entropy(counts),
    clonality       = clonality_index(counts),
    gini_simpson    = gini_simpson(counts),
    d50             = d50(counts)
  )
}

#' Compute diversity metrics for all samples in a dataset
#' @param tcr_df Full TCR data frame with patient_id, sample_type, clone_count
#' @return Tibble of per-sample diversity metrics
compute_all_diversity <- function(tcr_df) {
  tcr_df %>%
    group_by(patient_id, sample_type, clinical_response) %>%
    group_modify(~ {
      metrics <- compute_diversity(.x)
      as_tibble(metrics)
    }) %>%
    ungroup()
}
