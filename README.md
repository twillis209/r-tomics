# `r-tomics`

This R package is intended to house disparate functions I have found useful in processing data from GWAS and related study types. The name is an unspired portmanteau of my name and "'omics", and mention of the latter is aspirational at present as the code only relates to GWAS summary statistics.

# Development workflow

What I do after making changes to source:

```bash
  conda activate tomics-dev
  R -e "devtools::document()"
  R -e "attachment::att_amend_desc()"
  R CMD build .
  # No pdflatex available on my Mac atm
  R CMD check tomics_0.0.2.tar.gz --no-manual --no-build-vignettes
```

The YAML file defining the `tomics-dev` `conda` environment is located in the project's root directory.

The `recipe.yml` file needs to have the following fields updated for a new version:
- version field
- sha256 field

# Build process

Using `rattler-build`:

```bash
conda activate tomics-dev
rattler-build build --recipe recipe.yml --output-dir ../r-tomics
```
And then publishing on Anaconda:
```bash
rattler-build upload anaconda ../r-tomics/noarch/r-tomics-0.0.2-h60d57d3_0.conda --owner twillis209
```
