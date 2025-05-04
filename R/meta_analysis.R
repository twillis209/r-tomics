##' Update existing effect estimates and standard errors using inverse variance-weighted fixed-effects meta-analysis
##'
##' @param dat data.table containing effect estimates and standard errors
##' @param study_name character string containing the name of the study
##' @param effect_label character string containing the label for the effect estimate
##' @param se_label character string containing the label for the standard error
##' @param z2_label character string containing the label for the squared z-score
##' @param p_label character string containing the label for the p-value
##' @param wt_label character string containing the label for the weight
##' @return data.table with updated effect estimates, standard errors, z-scores, and p-values
##' @author Tom Willis
##' @import data.table
##' @importFrom stats pchisq
##' @export
update_meta_analysis_estimates <- function(dat, study_name, effect_label = "BETA", se_label = "SE", z2_label = "Z2", p_label = "P", wt_label = "wt") {
  beta.study <- se.study <- wt.study <- se <- wt <- p <- z2 <- NULL

  if (!data.table::is.data.table(dat)) {
    stop("dat must be a data.table")
  }

  beta_study_label <- paste(effect_label, study_name, sep = ".")
  se_study_label <- paste(se_label, study_name, sep = ".")
  wt_study_label <- paste(wt_label, study_name, sep = ".")

  # This is the first study
  if (!(effect_label %in% names(dat))) {
    dat[, `:=`(beta = beta.study, se = se.study),
      env = list(
        beta = effect_label,
        se = se_label,
        beta.study = beta_study_label,
        se.study = se_study_label
      )
    ]
  } else {
    # Already have estimate for these SNPs which we need to update
    dat[!is.na(beta), `:=`(wt = se^-2, wt.study = se.study^-2),
      env = list(
        beta = effect_label,
        se = se_label,
        wt = wt_label,
        se.study = se_study_label,
        wt.study = wt_study_label
      )
    ]
    dat[!is.na(beta) & !is.na(beta.study),
        `:=`(beta = (beta * wt + beta.study * wt.study) / (wt + wt.study), se = (wt + wt.study)^-0.5),
      env = list(
        beta = effect_label,
        se = se_label,
        beta.study = beta_study_label,
        wt = wt_label,
        wt.study = wt_study_label
      )
    ]

    # SNPs absent from existing data
    dat[is.na(beta) & !is.na(beta.study), `:=`(beta = beta.study, se = se.study),
      env = list(
        beta = effect_label,
        se = se_label,
        beta.study = beta_study_label,
        se.study = se_study_label
      )
    ]
  }

  dat[!is.na(beta), `:=`(z2 = (beta / se)^2),
    env = list(
      beta = effect_label,
      z2 = z2_label,
      se = se_label
    )
  ]

  dat[!is.na(beta), p := stats::pchisq(z2, df = 1, lower.tail = FALSE),
    env = list(
      beta = effect_label,
      p = p_label,
      z2 = z2_label
    )
  ]
}
