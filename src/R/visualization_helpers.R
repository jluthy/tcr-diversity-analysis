# ============================================================================
# visualization_helpers.R
# Common plotting functions for TCR repertoire analysis
# Author: Joshua Luthy
# ============================================================================

library(ggplot2)

# Custom theme matching portfolio site aesthetics
theme_tcr <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.background = element_rect(fill = "#0a0f14", color = NA),
      panel.background = element_rect(fill = "#0a0f14", color = NA),
      panel.grid.major = element_line(color = "#1e3348", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      text = element_text(color = "#e0eef5"),
      axis.text = element_text(color = "#7b9db5"),
      axis.title = element_text(color = "#7b9db5"),
      plot.title = element_text(color = "#e0eef5", face = "bold", size = 16),
      plot.subtitle = element_text(color = "#7b9db5", size = 11),
      legend.background = element_rect(fill = "#0f1923", color = "#1e3348"),
      legend.text = element_text(color = "#7b9db5"),
      legend.title = element_text(color = "#e0eef5")
    )
}

# Color palettes
tcr_colors <- list(
  sample_type = c("Apheresis" = "#00d4ff", "Product" = "#7b2fff"),
  response    = c("CR" = "#00ff9d", "PR" = "#00d4ff", "PD" = "#ff6b6b"),
  accent      = c("#00d4ff", "#7b2fff", "#00ff9d", "#ffd93d", "#ff6b6b")
)

# Diversity boxplot
plot_diversity_boxplot <- function(df, metric, ylab = metric) {
  ggplot(df, aes(x = sample_type, y = .data[[metric]], fill = sample_type)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
    scale_fill_manual(values = tcr_colors$sample_type) +
    labs(y = ylab, x = NULL) +
    theme_tcr() +
    theme(legend.position = "none")
}

# Cumulative clonal dominance plot
plot_cumulative_dominance <- function(top_clones_df) {
  cumulative <- top_clones_df %>%
    dplyr::filter(sample_type == "Product") %>%
    dplyr::group_by(patient_id) %>%
    dplyr::arrange(rank) %>%
    dplyr::mutate(cum_fraction = cumsum(clone_fraction))

  ggplot(cumulative, aes(x = rank, y = cum_fraction,
                          color = clinical_response, group = patient_id)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = tcr_colors$response) +
    labs(x = "Clonotype Rank", y = "Cumulative Fraction", color = "Response") +
    theme_tcr()
}
