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
##' @param y_limits limits for y axis
##' @param y_axis_breaks numeric vector containing coordinates at which to place breaks on the y-axis
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
                           title = "", y_limits = c(1, 1e-10), y_axis_breaks = 10^(-seq(0, 10, by = 1)),
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
    scale_y_neglog10(limits = y_limits, breaks = y_axis_breaks) +
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



  if (!is.null(lead_snps)) {
    merged_lead_snps <- merge(lead_snps, gwas_data[, .(chr, bp, bp_cum)], by = c("chr", "bp"))

    pl <- pl +
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
scale_y_neglog10 <- function(..., breaks = scales::trans_breaks(function(x) {
  log10(x) * -1
}, function(x) {
  10^(-1 * x)
})) {
  .x <- NULL

  ggplot2::scale_y_continuous(...,
    trans = neglog_trans(base = 10),
    breaks = breaks,
    labels = scales::trans_format(function(x) {
      log10(x) * -1
    }, scales::math_format(.x))
  )
}

##' Draw interactive plotly Manhattan plot for GWAS data
##'
##' Creates an interactive Manhattan plot using plotly that allows zooming and
##' panning to explore GWAS association signals. Includes optional p-value
##' filtering to improve performance with large datasets.
##'
##' @param dat data.table containing summary statistics with columns for
##'   chromosome, position, p-value, and optionally SNP ID
##' @param chr_col character string containing the column name for chromosome
##' @param bp_col character string containing the column name for position
##' @param p_col character string containing the column name for p-value
##' @param snp_col character string containing the column name for SNP ID (optional)
##' @param chromosomes character vector containing chromosomes to include in the plot
##' @param region_chr character string specifying a single chromosome to plot (optional).
##'   If provided with region_start and region_end, only this region will be plotted.
##' @param region_start numeric start position in base pairs for region to plot (optional)
##' @param region_end numeric end position in base pairs for region to plot (optional)
##' @param p_threshold numeric p-value threshold for filtering SNPs (default 1,
##'   no filtering). SNPs with p > p_threshold are excluded to improve performance.
##' @param palette character vector containing colors for chromosomes
##' @param title character string containing the title for the plot
##' @param genome_wide_line numeric p-value for genome-wide significance line
##'   (default 5e-8)
##' @param suggestive_line numeric p-value for suggestive significance line
##'   (default 1e-5)
##' @param point_size numeric size of points (default 3)
##' @return plotly object with interactive Manhattan plot
##' @author Tom Willis
##' @examples
##' \dontrun{
##' library(data.table)
##' 
##' # Create example GWAS data
##' gwas_data <- data.table(
##'   chromosome = rep(1:22, each = 10000),
##'   base_pair_location = rep(1:10000, 22) * 1000,
##'   p_value = runif(220000, 0, 1),
##'   snp_id = paste0("rs", 1:220000)
##' )
##' 
##' # Basic interactive Manhattan plot
##' fig <- draw_plotly_manhattan(
##'   dat = gwas_data,
##'   chr_col = "chromosome",
##'   bp_col = "base_pair_location",
##'   p_col = "p_value",
##'   title = "My GWAS Study"
##' )
##' fig
##' 
##' # With p-value filtering for better performance (only show p < 0.01)
##' fig <- draw_plotly_manhattan(
##'   dat = gwas_data,
##'   chr_col = "chromosome",
##'   bp_col = "base_pair_location",
##'   p_col = "p_value",
##'   snp_col = "snp_id",
##'   p_threshold = 0.01,
##'   title = "Filtered Manhattan Plot"
##' )
##' fig
##' 
##' # Customize colors and significance lines
##' fig <- draw_plotly_manhattan(
##'   dat = gwas_data,
##'   p_threshold = 0.001,
##'   palette = c("#1f77b4", "#ff7f0e"),
##'   genome_wide_line = 5e-8,
##'   suggestive_line = 1e-5,
##'   point_size = 4
##' )
##' fig
##' 
##' # Plot a specific region (e.g., chromosome 6, 25-35 Mb)
##' fig <- draw_plotly_manhattan(
##'   dat = gwas_data,
##'   region_chr = "6",
##'   region_start = 25000000,
##'   region_end = 35000000,
##'   title = "Chromosome 6: 25-35 Mb"
##' )
##' fig
##' 
##' # Save to HTML file for viewing in browser
##' fig <- draw_plotly_manhattan(dat = gwas_data, p_threshold = 0.01)
##' save_plotly_html(fig, "manhattan_plot.html")
##' }
##' @export
draw_plotly_manhattan <- function(dat,
                                  chr_col = "chromosome",
                                  bp_col = "base_pair_location",
                                  p_col = "p_value",
                                  snp_col = NULL,
                                  chromosomes = as.character(1:22),
                                  region_chr = NULL,
                                  region_start = NULL,
                                  region_end = NULL,
                                  p_threshold = 1,
                                  palette = c("#E69F00", "#56B4E9"),
                                  title = "Manhattan Plot",
                                  genome_wide_line = 5e-8,
                                  suggestive_line = 1e-5,
                                  point_size = 3) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function. Please install it with: install.packages('plotly')")
  }
  
  chr <- bp <- p_value <- neglog_p <- bp_cum <- color_group <- NULL
  
  dat <- data.table::copy(dat)
  
  data.table::setnames(dat, c(chr_col, bp_col, p_col), c("chr", "bp", "p_value"))
  
  if (!is.null(region_chr)) {
    dat <- dat[chr == region_chr]
    if (!is.null(region_start)) {
      dat <- dat[bp >= region_start]
    }
    if (!is.null(region_end)) {
      dat <- dat[bp <= region_end]
    }
    chromosomes <- region_chr
  } else {
    dat <- dat[chr %in% chromosomes]
  }
  
  if (p_threshold < 1) {
    dat <- dat[p_value <= p_threshold]
  }
  
  dat[, chr := as.integer(chr)]
  data.table::setorder(dat, chr, bp)
  
  dat[, neglog_p := -log10(p_value)]
  
  if (!is.null(region_chr)) {
    dat[, bp_cum := bp]
    axis_set <- data.table::data.table(center = mean(dat$bp), chr = unique(dat$chr))
  } else {
    bp_offsets <- dat[, .(max_bp = max(bp)), by = chr]
    bp_offsets[, bp_add := cumsum(as.numeric(data.table::shift(max_bp, n = 1, fill = 0, type = "lag")))]
    
    dat[bp_offsets, on = "chr", bp_cum := bp + i.bp_add]
    
    axis_set <- dat[, .(center = mean(bp_cum), chr = unique(chr)), by = chr]
  }
  
  dat[, color_group := (chr %% 2) + 1]
  
  hover_text <- if (!is.null(snp_col) && snp_col %in% names(dat)) {
    sprintf("SNP: %s<br>Chr: %s<br>Pos: %s<br>P-value: %.2e",
            dat[[snp_col]], dat$chr, dat$bp, dat$p_value)
  } else {
    sprintf("Chr: %s<br>Pos: %s<br>P-value: %.2e",
            dat$chr, dat$bp, dat$p_value)
  }
  
  fig <- plotly::plot_ly()
  
  for (grp in unique(dat$color_group)) {
    dat_subset <- dat[color_group == grp]
    fig <- plotly::add_trace(
      fig,
      data = dat_subset,
      x = ~bp_cum,
      y = ~neglog_p,
      type = "scatter",
      mode = "markers",
      marker = list(
        size = point_size,
        color = palette[grp],
        opacity = 0.7
      ),
      text = hover_text[dat$color_group == grp],
      hovertemplate = "%{text}<extra></extra>",
      showlegend = FALSE
    )
  }
  
  if (!is.null(genome_wide_line) && genome_wide_line > 0) {
    fig <- plotly::add_trace(
      fig,
      x = c(min(dat$bp_cum), max(dat$bp_cum)),
      y = c(-log10(genome_wide_line), -log10(genome_wide_line)),
      type = "scatter",
      mode = "lines",
      line = list(color = "red", dash = "dash", width = 1),
      name = sprintf("Genome-wide (p=%.0e)", genome_wide_line),
      hoverinfo = "name",
      showlegend = TRUE
    )
  }
  
  if (!is.null(suggestive_line) && suggestive_line > 0) {
    fig <- plotly::add_trace(
      fig,
      x = c(min(dat$bp_cum), max(dat$bp_cum)),
      y = c(-log10(suggestive_line), -log10(suggestive_line)),
      type = "scatter",
      mode = "lines",
      line = list(color = "blue", dash = "dash", width = 1),
      name = sprintf("Suggestive (p=%.0e)", suggestive_line),
      hoverinfo = "name",
      showlegend = TRUE
    )
  }
  
  xaxis_config <- if (!is.null(region_chr)) {
    list(
      title = "Position (bp)",
      zeroline = FALSE
    )
  } else {
    list(
      title = "Chromosome",
      tickmode = "array",
      tickvals = axis_set$center,
      ticktext = axis_set$chr,
      zeroline = FALSE
    )
  }
  
  fig <- plotly::layout(
    fig,
    title = title,
    xaxis = xaxis_config,
    yaxis = list(
      title = "-log<sub>10</sub>(P-value)",
      zeroline = FALSE
    ),
    hovermode = "closest",
    dragmode = "zoom"
  )
  
  return(fig)
}

##' Save plotly figure as standalone HTML file
##'
##' Saves an interactive plotly figure as a self-contained HTML file that can
##' be opened in any web browser without needing R or any other dependencies.
##'
##' @param fig plotly figure object to save
##' @param file character string specifying the output file path
##' @param title character string for the HTML page title (optional)
##' @return NULL (invisibly). Called for side effect of saving HTML file.
##' @author Tom Willis
##' @examples
##' \dontrun{
##' library(data.table)
##' 
##' gwas_data <- data.table(
##'   chromosome = rep(1:22, each = 10000),
##'   base_pair_location = rep(1:10000, 22) * 1000,
##'   p_value = runif(220000, 0, 1),
##'   snp_id = paste0("rs", 1:220000)
##' )
##' 
##' fig <- draw_plotly_manhattan(dat = gwas_data, p_threshold = 0.01)
##' save_plotly_html(fig, "manhattan_plot.html")
##' 
##' # View in browser
##' browseURL("manhattan_plot.html")
##' }
##' @export
save_plotly_html <- function(fig, file, title = "Interactive Plot") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function. Please install it with: install.packages('plotly')")
  }
  
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("Package 'htmlwidgets' is required for this function. Please install it with: install.packages('htmlwidgets')")
  }
  
  htmlwidgets::saveWidget(
    widget = fig,
    file = file,
    selfcontained = TRUE,
    title = title
  )
  
  message(sprintf("Plot saved to: %s", file))
  invisible(NULL)
}

