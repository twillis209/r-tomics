##' Process summary statistics for Manhattan plot
##'
##' Code originally from Daniel Roelfs, modified by Tom Willis
##' https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##'
##' @param dat data.table containing summary statistics
##' @param chr_col character string containing the column name for chromosome
##' @param bp_col character string containing the column name for position
##' @param p_col character string containing the column name for p-value
##' @param x_chr_no index of X chromosome
##' @return list containing data.table with updated chromosome and position
##' columns, axis set, and y-axis limits
##' @author Daniel Roelfs
##' @author Tom Willis
##' @import data.table
##' @importFrom dplyr group_by summarise mutate select inner_join filter pull
##' @importFrom stats lag
##' @export
process_sumstats_for_manhattan <- function(dat, chr_col = "chromosome", bp_col = "base_pair_location", p_col = "p_value", x_chr_no = 23) {
  chr <- bp <- max_bp <- bp_add <- bp_cum <- p <- NULL

  dat <- data.table::copy(dat)

  data.table::setnames(dat, c(chr_col, bp_col, p_col), c("chr", "bp", "p"))

  dat[chr == "X", chr := x_chr_no]

  dat[, chr := as.integer(chr)]

  dat <- dat[chr %in% seq(1, 23)]

  data_cum <- dat %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_bp = max(bp)) %>%
    dplyr::mutate(bp_add = stats::lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    dplyr::select(chr, bp_add)

  data <- dat %>%
    dplyr::inner_join(data_cum, by = "chr") %>%
    dplyr::mutate(bp_cum = as.numeric(bp) + as.numeric(bp_add))

  axis_set <- data %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(center = mean(bp_cum))

  return(list(data = data, axis_set = axis_set))
}

##' Draw Manhattan plot using processed summary statistics
##'
##' Code originally from Daniel Roelfs, modified by Tom Willis
##' https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##'
##' @param processed_sumstats list containing data.table with updated chromosome and position values
##' @param palette character vector containing colors for chromosomes
##' @param title character string containing the title for the plot
##' @param y_axis_break numeric vector containing coordinates at which to break the y-axis
##' @param y_limits limits for y axis
##' @author Daniel Roelfs
##' @author Tom Willis
##' @importFrom ggplot2 ggplot geom_hline geom_point scale_x_continuous scale_color_manual scale_size_continuous scale_y_continuous labs theme ggtitle
##' @importFrom ggtext element_markdown
##' @importFrom ggbreak scale_y_break
##' @export
draw_manhattan <- function(processed_sumstats, palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), title = "", y_limits = c(1, 1e-10), y_axis_break = NULL) {
  chr <- bp_cum <- p <- NULL

  gwas_data <- processed_sumstats$data
  axis_set <- processed_sumstats$axis_set

  pl <- ggplot2::ggplot(ggplot2::aes(x = bp_cum, y = p, color = as.factor(chr)), data = gwas_data) +
    ggplot2::geom_hline(yintercept = 5e-8, color = "grey40", linetype = "dashed") +
    ggplot2::geom_point(size = 0.3) +
    ggplot2::scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    ggplot2::scale_color_manual(values = rep(palette, unique(length(axis_set$chr)))) +
    ggplot2::scale_size_continuous(range = c(0.5, 3)) +
    scale_y_neglog10(limits = y_limits) +
    ggplot2::labs(
      x = NULL,
      y = "-log<sub>10</sub>(p)"
    ) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggplot2::element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    ggplot2::ggtitle(title)

  if (!is.null(y_axis_break)) {
    pl <- pl + ggbreak::scale_y_break(y_axis_break)
  }

  pl
}

##' -log transformation
##'
##' @param base base of the logarithm
##' @importFrom scales trans_new
neglog_trans <- function(base = exp(1)) {
  scales::trans_new("neglog",
    transform = function(x) -log(x, base),
    inverse = function(x) base^(-x),
    domain = c(1e-100, Inf)
  )
}

##' -log10 transformation of the x-axis
##'
##' @param  ... additional arguments passed to [ggplot2::scale_x_continuous()]
##' @return ggplot2::scale_x_continuous() with -log10 transformation
##'
##' @author Tom Willis
##' @importFrom scales trans_breaks trans_format math_format
##' @importFrom ggplot2 scale_x_continuous
##' @export
scale_x_neglog10 <- function(...) {
  ggplot2::scale_x_continuous(...,
    trans = neglog_trans(base = 10),
    breaks = scales::trans_breaks(function(x) {
      log10(x) * -1
    }, function(x) {
      10^(-1 * x)
    }),
    labels = scales::trans_format(function(x) {
      log10(x) * -1
    }, scales::math_format(.x))
  )
}

##' -log10 transformation of the y-axis
##'
##' @param  ... additional arguments passed to [ggplot2::scale_y_continuous()]
##' @return ggplot2::scale_y_continuous() with -log10 transformation
##' @author Tom Willis
##' @importFrom scales trans_breaks trans_format math_format
##' @importFrom ggplot2 scale_y_continuous
##' @export
scale_y_neglog10 <- function(...) {
  .x <- NULL

  ggplot2::scale_y_continuous(...,
    trans = neglog_trans(base = 10),
    breaks = scales::trans_breaks(function(x) {
      log10(x) * -1
    }, function(x) {
      10^(-1 * x)
    }),
    labels = scales::trans_format(function(x) {
      log10(x) * -1
    }, scales::math_format(.x))
  )
}
