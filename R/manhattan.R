##' Process summary statistics for Manhattan plot
##'
##' Code originally from Daniel Roelfs, modified by Tom Willis
##' https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##'
##' @param dat data.table containing summary statistics
##' @param chr_col character string containing the column name for chromosome
##' @param bp_col character string containing the column name for position
##' @param stat_cols character vector containing the column name(s) for per-variant statistics such as p-values, chi-squared test statistics, etc.
##' @param chromosomes character vector containing chromosomes to include in the plot
##' @param inter_chr_offset size of gap between chromosomes
##' @return list containing data.table with updated chromosome and position
##' columns, axis set, and y-axis limits
##' @author Daniel Roelfs
##' @author Tom Willis
##' @import data.table
##' @export
process_sumstats_for_manhattan <- function(dat, chr_col = "chromosome", bp_col = "base_pair_location", stat_cols = "p_value", chromosomes = as.character(1:23), inter_chr_offset = 5e7) {
  chr <- bp <- max_bp <- bp_add <- i.bp_add <- bp_cum <- inter_chr <- NULL

  dat <- data.table::copy(dat)

  data.table::setnames(dat, c(chr_col, bp_col), c("chr", "bp"))

  dat <- dat[chr %in% chromosomes]

  dat[, chr := as.integer(chr)]

  data.table::setorder(dat, chr)

  bp_offsets <- dat[, .(max_bp = max(bp)), by = chr] %>%
    .[, .(chr, bp_add = cumsum(as.numeric(data.table::shift(max_bp, n = 1, fill = 0, type = "lag"))))]

  bp_offsets[, inter_chr := inter_chr_offset * (0:(.N - 1))]

  bp_offsets[, bp_add := bp_add + inter_chr]

  dat[bp_offsets, on = "chr", bp_cum := bp + i.bp_add]

  axis_set <- dat[, .(center = mean(bp_cum)), by = chr]

  cols_to_keep <- c("chr", "bp", "bp_cum", stat_cols)

  return(list(dat = dat[, ..cols_to_keep], axis_set = axis_set))
}

##' Draw Manhattan plot using processed summary statistics
##'
##' Code originally from Daniel Roelfs, modified by Tom Willis
##' https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##'
##' @param processed_sumstats list containing data.table with updated chromosome and position values
##' @param stat_col character string containing the column name for the test statistic to display
##' @param palette character vector containing colors for chromosomes
##' @param title character string containing the title for the plot
##' @param y_axis_break numeric vector containing coordinates at which to break the y-axis
##' @param y_limits limits for y axis
##' @param lead_snps data.table containing lead SNPs to annotate with gene names
##' @param repel_args list containing arguments for ggrepel::geom_text_repel()
##' @author Daniel Roelfs
##' @author Tom Willis
##' @importFrom ggplot2 ggplot aes geom_hline geom_point scale_x_continuous
##' scale_color_manual scale_size_continuous scale_y_continuous labs theme
##' expansion
##' @importFrom ggrepel geom_text_repel
##' @importFrom ggtext element_markdown
##' @importFrom ggbreak scale_y_break
##' @importFrom rlang exec sym
##' @export
draw_manhattan <- function(processed_sumstats,
                           stat_col = "p",
                           palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                           title = "", y_limits = c(1, 1e-10), y_axis_break = NULL,
                           lead_snps = NULL,
                           repel_args = list(
                             hjust = -0.2,
                             vjust = 0,
                             size = 6,
                             angle = 45,
                             colour = "black",
                             direction = "x",
                             nudge_y = 0
                           )) {
  chr <- bp <- bp_cum <- NULL

  if (!is.null(lead_snps)) {
    lead_snps <- data.table::copy(lead_snps)
    if (!all(c("chr", "bp", stat_col, "gene") %in% colnames(lead_snps))) {
      stop(sprintf("lead_snps must contain columns 'chr', 'bp', '%s', and 'gene'", stat_col))
    }
  }

  gwas_data <- processed_sumstats$dat
  axis_set <- processed_sumstats$axis_set

  pl <- ggplot2::ggplot(ggplot2::aes(x = bp_cum, y = !!rlang::sym(stat_col), color = as.factor(chr)), data = gwas_data) +
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

  if (!is.null(lead_snps)) {
    merged_lead_snps <- merge(lead_snps, gwas_data[, .(chr, bp, bp_cum)], by = c("chr", "bp"))

    pl <- pl +
      ## scale_y_neglog10(limits = y_limits, expand = ggplot2::expansion(mult = c(0, y_axis_space_mult))) +
      ggplot2::geom_point(ggplot2::aes(x = bp_cum, y = !!rlang::sym(stat_col), color = as.factor(chr)), size = 0.9, pch = 21, colour = "black", data = merged_lead_snps) +
      ggplot2::coord_cartesian(clip = "off") +
      rlang::exec(
        ggrepel::geom_text_repel,
        mapping = ggplot2::aes(label = .data$gene),
        data = merged_lead_snps,
        !!!repel_args,
      )
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
  .x <- NULL

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
