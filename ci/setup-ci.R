# ci/setup-ci.R

if (!"rnaturalearthhires" %in% installed.packages()) {
  try(devtools::install_github("ropensci/rnaturalearthhires"), silent = TRUE)
}
