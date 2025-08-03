##' Estimate continuous trait standard deviation given vectors of squared standard errors, MAF estimates, and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##'
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' @importFrom stats lm coef
##' @author Chris Wallace
##' @export
sdY.est <- function(vbeta, maf, n) {
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- stats::lm(nvx ~ oneover - 1)
  cf <- stats::coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}
