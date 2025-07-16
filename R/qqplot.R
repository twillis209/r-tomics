##' Draw Q-Q plot.
##' @description Draws a Q-Q plot of p-values from a data.frame.
##' @param daf data.frame containing p-values.
##' @param p_col p-value column label
##' @return ggplot object
##' @importFrom ggplot2 ggplot aes geom_abline labs coord_fixed stat_qq
##' @importFrom ggtext element_markdown
##' @importFrom stats qunif
##' @export
qqplot <- function(daf, p_col) {
  ggplot2::ggplot(data = daf) +
    ggplot2::stat_qq(ggplot2::aes(sample = !!rlang::sym(p_col)), geom = "line", distribution = stats::qunif) +
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

#' Draw stratified Q-Q plot.
#'
#' @description Draws a stratified Q-Q plot of p-values from a data.frame, stratified by a secondary variable.
#' @param daf data.frame containing p-values for primary and secondary variables.
#' @param p_col p-value column label for primary variable
#' @param q_col p-value column label for secondary variable for stratification
#' @param breaks numeric vector of breaks for stratification
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_abline labs coord_fixed stat_qq
#' @importFrom stats qunif quantile
#' @importFrom ggtext element_markdown
#' @importFrom stats qunif
#' @importFrom dplyr mutate select
#' @export
stratified_qqplot <- function(daf, p_col, q_col, breaks = c(0, 10^(c(-5:-1, 0)))) {
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
    ggplot2::stat_qq(ggplot2::aes(sample = !!rlang::sym(p_col), col = Q), geom = "line", distribution = stats::qunif) +
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
