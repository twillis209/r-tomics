# `r-tomics`

This R package is intended to house disparate functions I have found useful in processing data from GWAS and related study types. The name is an unspired portmanteau of my name and "'omics", and mention of the latter is aspirational at present as the code only relates to GWAS summary statistics.

NB: As of October '25 or so, many of the new features of this package have been vibe-coded with Claude Code 4.5 (that's right, I'm afraid, it's AI slop now).

# Development workflow

See `release.sh` for the build and deployment script.

The YAML file defining the `tomics-dev` `conda` environment is located in the project's root directory. Note that at the moment it's still necessary to update the dependencies in `recipe.yml`, `release.sh` doesn't handle that yet.

# Adding `vdiffr` tests

Upon adding a new `vdiffr::expect_doppelganger` tests case:
* run `devtools::test()`
* if new, review new snapshots under `tests/testthat/_snaps`
* if a regression is detected, run `testthat::snapshot_review()` to review changes

# TODO Update `recipe.yml`

I use `attachments` to update `DESCRIPTION` with new packages, but I should do the same to update `recipe.yml`.


