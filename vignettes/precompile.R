# Precompute the vignettes whose experiments take several minutes.
#
# Each `<name>.Rmd.orig` is knitted into a static `<name>.Rmd` (with figures in
# `vignettes/figures/`), which is what R CMD build and CRAN see. Run this
# script from the package root after changing the sampler or the vignettes:
#
#   Rscript vignettes/precompile.R
#
# Set the environment variable RTWALK_VIGNETTE_CACHE to keep the knitr cache
# outside the package directory.

pkgload::load_all(".")

old_wd <- setwd("vignettes")
on.exit(setwd(old_wd))

for (name in c("reproducing_christen_fox", "sampler_comparison")) {
  source_file <- paste0(name, ".Rmd.orig")
  if (file.exists(source_file)) {
    knitr::knit(source_file, paste0(name, ".Rmd"))
  }
}
