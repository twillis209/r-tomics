#' @title Clumping SNPs to identify lead SNPs in association signals
#'
#' @description This function performs clumping of SNPs based on their physical distance and p-values.
#'
##' @param dat data.table containing summary statistics
##' @param chr_col character string containing the column name for chromosome
##' @param bp_col character string containing the column name for position
##' @param p_col character string containing the column name for p-value
##' @param distance_window numeric value specifying the distance window for clumping
##' @param index_threshold numeric value specifying the p-value threshold for candidate SNPs
##'
##' @return data.table with clumped SNPs
##' @author Tom Willis
##'
#' @import data.table
#' @export
distance_clump <- function(dat, chr_col, bp_col, p_col, distance_window = 1e6, index_threshold = 5e-8) {
  chr <- p <- bp <- NULL

  dat <- dat[p <= index_threshold, env = list(p = p_col)]

  data.table::setorderv(dat, p_col)

  # Loop through each chromosome
  for (i in unique(dat[[chr_col]])) {
    j <- 1

    while (j <= dat[chr == i, .N, env = list(chr = chr_col)]) {
      # Get the position of the current SNP
      pos <- dat[chr == i, env = list(chr = chr_col)][j, bp, env = list(bp = bp_col)]

      # Find and remove SNPs within the distance window
      dat <- dat[!(chr == i & bp != pos & bp %between% c(max(pos - distance_window / 2, 0), pos + distance_window / 2)),
        env = list(chr = chr_col, bp = bp_col)
      ]

      # Move to the next SNP
      j <- j + 1
    }
  }

  return(dat)
}
