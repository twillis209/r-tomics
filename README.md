# `r-tomics`

This R package is intended to house disparate functions I have found useful in processing data from GWAS and related study types. The name is an unspired portmanteau of my name and "'omics", and mention of the latter is aspirational at present as the code only relates to GWAS summary statistics.

# Development workflow

```bash
  conda activate tomics-dev
  # Need to commit changes prior to this as otherwise usethis::use_version() will complain
  R -e "usethis::use_version('minor')"
  R -e "devtools::document()"
  # NB: updates dependencies etc.
  R -e "attachment::att_amend_desc()"
  R CMD build .
  latest=$(ls -v tomics_*.tar.gz | tail -n 1)
  R CMD check $latest --no-manual --no-build-vignettes
  version=$(echo "$latest" | sed -E 's/tomics_(.*)\.tar\.gz/\1/')
  sha256=$(sha256sum "$latest" | awk '{print $1}')
  # Only want to change first version
  # NB: double backspaces for group because this is the osx version
  sed -i -E "s/^([[:space:]]*)version: [0-9]+\.[0-9]+\.[0-9]+/\\1version: $version/" recipe.yml
  sed -i -E "s/sha256: .*/sha256: $sha256/" recipe.yml
```

The YAML file defining the `tomics-dev` `conda` environment is located in the project's root directory. 

At the moment it's still necessary to update the dependencies in `recipe.yml`, should be able to fix it by parsing the YAML.

It's necessary to `commit` again after the earlier steps, usually just amend the initial commit.

```bash
git tag "v$version"
git push origin master
git push -f origin "v$version"
# --notes can be used to specify notes string
gh release create "v$version" $latest --notes-from-tag
```

# Build process

Using `rattler-build`:

```bash
rattler-build build --recipe recipe.yml --output-dir ../r-tomics
```
And then publishing on Anaconda:
```bash
rattler-build upload anaconda $(ls ../r-tomics/linux-64/r-tomics-$version-*.conda) --owner twillis209
```
Or
``` bash
rattler-build upload anaconda $(ls ../r-tomics/osx-arm64/r-tomics-$version-*.conda) --owner twillis209
```
