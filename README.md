# `r-tomics`

This R package is intended to house disparate functions I have found useful in processing data from GWAS and related study types. The name is an unspired portmanteau of my name and "'omics", and mention of the latter is aspirational at present as the code only relates to GWAS summary statistics.

# Development workflow

What I do after making changes to source:

- increment version number

```bash
  conda activate tomics-dev
  R -e "devtools::document()"
  R -e "attachment::att_amend_desc()"
  R CMD build .
  # No pdflatex available on my Mac atm
  R CMD check $(ls -v tomics_*.tar.gz | tail -n 1) --no-manual --no-build-vignettes
```

The YAML file defining the `tomics-dev` `conda` environment is located in the project's root directory.

The `recipe.yml` file needs to have the following fields updated for a new version:
- version field
- sha256 field

```bash
git tag v0.0.x
git push v0.0.x
# --notes can be used to specify notes string
gh release create v0.0.x tomics_0_0.x.tar.gz --notes-from-tag
```

# Build process

Using `rattler-build`:

```bash
conda activate tomics-dev
rattler-build build --recipe recipe.yml --output-dir ../r-tomics
```
And then publishing on Anaconda:
```bash
rattler-build upload anaconda ../r-tomics/noarch/r-tomics-0.0.1-h4616a5c_0.conda --owner twillis209
```
