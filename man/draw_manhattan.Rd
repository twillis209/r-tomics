% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manhattan.R
\name{draw_manhattan}
\alias{draw_manhattan}
\title{Draw Manhattan plot using processed summary statistics}
\usage{
draw_manhattan(
  processed_sumstats,
  stat_col = "p",
  palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7"),
  title = "",
  y_limits = c(1, 1e-10),
  y_axis_break = NULL,
  lead_snps = NULL,
  repel_args = list(hjust = -0.2, vjust = 0, size = 6, angle = 45, colour = "black",
    direction = "x", nudge_y = 0)
)
}
\arguments{
\item{processed_sumstats}{list containing data.table with updated chromosome and position values}

\item{stat_col}{character string containing the column name for the test statistic to display}

\item{palette}{character vector containing colors for chromosomes}

\item{title}{character string containing the title for the plot}

\item{y_limits}{limits for y axis}

\item{y_axis_break}{numeric vector containing coordinates at which to break the y-axis}

\item{lead_snps}{data.table containing lead SNPs to annotate with gene names}

\item{repel_args}{list containing arguments for ggrepel::geom_text_repel()}
}
\description{
Code originally from Daniel Roelfs, modified by Tom Willis
https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
}
\author{
Daniel Roelfs

Tom Willis
}
