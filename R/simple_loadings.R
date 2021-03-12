#' @title HotLoadings.plot_loadings_simple
#' Plot Loadings the simple way
#' @importFrom mixOmics plotLoadings
#' @importFrom RColorBrewer brewer.pal
#' @inheritParams HotLoadings.names
#' @inheritParams HotLoadings.plot_loadings
#' @inheritParams mixOmics::plotLoadings
#' @param colors Vector of Hexadecimal codes for colors. Must be the same length
#' of Y levels.
#' @return The function return a ggplot2 loadings object.

HotLoadings.plot_loadings_simple <- function(PSOBJ, format = c("short", "last","long"),
                                             data.splsda,
                                             method = "mean",
                                             contrib = "max",
                                             Y_name, comp,
                                             ndisplay = 15,
                                             offset = 0.005,
                                             colors = NULL) {
  names_df <- cbind(taxa_names(PSOBJ), HotLoadings.names(PSOBJ, format = format))

  plot_data <- mixOmics::plotLoadings(data.splsda,
    comp = comp, method = "mean",
    contrib = "max", ndisplay = ndisplay,
    plot = FALSE
  )

  df_loadings <- plot_data[, c("GroupContrib", "color", "importance")]
  df_ordered <- df_loadings[order(df_loadings$importance, decreasing = FALSE), ]
  Y <- unlist(phyloseq::sample_data(PSOBJ)[, Y_name])
  if (is(Y, "factor")) {
    Y_levels <- levels(Y)
  } else {
    Y_levels <- unique(Y)
  }
  df_ordered$GroupContrib <- factor(df_ordered$GroupContrib, levels = Y_levels, ordered = TRUE)
  if (length(colors) != length(Y_levels)) {
    stop("Number of colors and number of levels are different.")
  }
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(max(length(Y_levels), 3), name = "Set1")[1:length(Y_levels)]
  }
  df_ordered$color <- factor(df_ordered$GroupContrib, levels = Y_levels, labels = colors)

  p_loadings <- ggplot(df_ordered, aes(x = rownames(df_ordered), y = importance, fill = GroupContrib)) +
    geom_col() +
    geom_text(aes(
      label = names_df[match(rownames(df_ordered), names_df[, 1]), 2],
      y = ifelse(importance > 0, -offset, offset)
    ),
    hjust = ifelse(df_ordered$importance > 0, "right", "left"),
    vjust = "mid"
    ) +
    scale_x_discrete(limits = rownames(df_ordered)) +
    scale_y_continuous(limits = c(-max(abs(df_ordered$importance)), max(abs(df_ordered$importance)))) +
    scale_fill_manual(values = levels(df_ordered$color), breaks = levels(df_ordered$GroupContrib)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank()
    ) +
    coord_flip() +
    labs(x = "Taxonomic features", y = "Loading", fill = paste0(Y_name, ": ")) +
    ggtitle(
      label = "sPLS-DA",
      subtitle = paste0("Loadings - ", comp, " component")
    )
  return(p_loadings)
}
