##' Draw Q-Q plot.
##' @description Draws a Q-Q plot of p-values from a data.frame.
##' @param daf data.frame containing p-values.
##' @param p_cols p-value column names
##' @param p_col_labels p-value column labels to apply to the plot legend
##' @param subset_col optional character string specifying a column to use for subsetting data
##' @param subset_values optional vector of values in subset_col to plot as separate lines
##' @param subset_labels optional vector of labels for the subsets (must match length of subset_values)
##' @param geom character string specifying the type of geom to use for the Q-Q plot (default is "line")
##' @return ggplot object
##' @examples
##' \dontrun{
##' library(data.table)
##' 
##' # Example 1: Multiple columns
##' gwas_data <- data.table(
##'   p_discovery = runif(10000),
##'   p_replication = runif(10000)
##' )
##' qqplot(gwas_data, 
##'        p_cols = c("p_discovery", "p_replication"),
##'        p_col_labels = c("Discovery", "Replication"))
##' 
##' # Example 2: Subsets of single column
##' gwas_data <- data.table(
##'   p_value = runif(10000),
##'   chromosome = sample(1:22, 10000, replace = TRUE)
##' )
##' qqplot(gwas_data,
##'        p_cols = "p_value",
##'        subset_col = "chromosome",
##'        subset_values = c(1, 2, 3),
##'        subset_labels = c("Chr 1", "Chr 2", "Chr 3"))
##' 
##' # Example 3: Multiple columns with subsetting
##' gwas_data <- data.table(
##'   p_males = runif(10000),
##'   p_females = runif(10000),
##'   significant = sample(c(TRUE, FALSE), 10000, replace = TRUE)
##' )
##' qqplot(gwas_data,
##'        p_cols = c("p_males", "p_females"),
##'        p_col_labels = c("Males", "Females"),
##'        subset_col = "significant",
##'        subset_values = c(TRUE, FALSE),
##'        subset_labels = c("Significant", "Not Significant"))
##' }
##' @importFrom ggplot2 ggplot aes geom_abline labs coord_fixed stat_qq
##' @importFrom ggtext element_markdown
##' @importFrom stats qunif
##' @importFrom tidyr pivot_longer
##' @importFrom tidyselect all_of
##' @importFrom dplyr mutate filter
##' @export
qqplot <- function(daf, p_cols, p_col_labels = p_cols, subset_col = NULL, 
                   subset_values = NULL, subset_labels = NULL, geom = "line") {
  
  # If subsetting is requested
  if (!is.null(subset_col) && !is.null(subset_values)) {
    if (is.null(subset_labels)) {
      subset_labels <- as.character(subset_values)
    }
    
    if (length(subset_values) != length(subset_labels)) {
      stop("subset_values and subset_labels must have the same length")
    }
    
    # Create a list to hold data for each subset
    subset_list <- list()
    
    for (i in seq_along(subset_values)) {
      subset_data <- daf[daf[[subset_col]] == subset_values[i], , drop = FALSE]
      
      subset_long <- tidyr::pivot_longer(
        subset_data,
        cols = tidyselect::all_of(p_cols),
        names_to = "Statistic",
        values_to = "value"
      ) %>%
        dplyr::mutate(
          Statistic = factor(Statistic, levels = p_cols, labels = p_col_labels),
          Subset = subset_labels[i]
        )
      
      subset_list[[i]] <- subset_long
    }
    
    # Combine all subsets
    plot_data <- do.call(rbind, subset_list)
    
    # Create combined label for color aesthetic
    plot_data <- dplyr::mutate(plot_data, 
                                Combined = interaction(Statistic, Subset, sep = " - "))
    
    pl <- ggplot2::ggplot(plot_data) +
      ggplot2::stat_qq(ggplot2::aes(sample = value, col = Combined), 
                       geom = geom, distribution = stats::qunif) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
      scale_x_neglog10() +
      scale_y_neglog10() +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(
        x = "Expected -log<sub>10</sub>(p)",
        y = "Observed -log<sub>10</sub>(p)",
        colour = "Group"
      ) +
      ggplot2::theme(
        axis.title = ggtext::element_markdown(),
      )
    
  } else {
    # Original behavior: just plot multiple columns
    pl <- tidyr::pivot_longer(daf, cols = tidyselect::all_of(p_cols), 
                              names_to = "Statistic", values_to = "value") %>%
      dplyr::mutate(Statistic = factor(Statistic, levels = p_cols, labels = p_col_labels)) %>%
      ggplot2::ggplot() +
      ggplot2::stat_qq(ggplot2::aes(sample = value, col = Statistic), 
                       geom = geom, distribution = stats::qunif) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
      scale_x_neglog10() +
      scale_y_neglog10() +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(
        x = "Expected -log<sub>10</sub>(p)",
        y = "Observed -log<sub>10</sub>(p)",
        colour = "Statistic"
      ) +
      ggplot2::theme(
        axis.title = ggtext::element_markdown(),
      )
  }
  
  return(pl)
}

#' Draw stratified Q-Q plot.
#'
#' @description Draws a stratified Q-Q plot of p-values from a data.frame, stratified by a secondary variable.
#' @param daf data.frame containing p-values for primary and secondary variables.
#' @param p_col p-value column label for primary variable
#' @param q_col p-value column label for secondary variable for stratification
#' @param breaks numeric vector of breaks for stratification
#' @param geom character string specifying the type of geom to use for the Q-Q plot (default is "line")
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_abline labs coord_fixed stat_qq
#' @importFrom stats qunif quantile
#' @importFrom ggtext element_markdown
#' @importFrom stats qunif
#' @importFrom dplyr mutate select
#' @export
stratified_qqplot <- function(daf, p_col, q_col, breaks = c(0, 10^(c(-5:-1, 0))), geom = "line") {
  tibble::as_tibble(daf) %>%
    dplyr::mutate(
      qbreaks = list(stats::quantile(!!rlang::sym(q_col), probs = breaks, na.rm = TRUE)),
      Q = cut(!!rlang::sym(q_col),
        breaks = qbreaks[[1]],
        labels = as.character(breaks[2:length(breaks)]),
        include.lowest = TRUE
      )
    ) %>%
    dplyr::select(-qbreaks) %>%
    ggplot2::ggplot() +
    ggplot2::stat_qq(ggplot2::aes(sample = !!rlang::sym(p_col), col = Q), geom = geom, distribution = stats::qunif) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
    scale_x_neglog10() +
    scale_y_neglog10() +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      x = "Expected -log<sub>10</sub>(p)",
      y = "Observed -log<sub>10</sub>(p)",
      colour = "Statistic"
    ) +
    ggplot2::theme(
      axis.title = ggtext::element_markdown(),
    )
}
