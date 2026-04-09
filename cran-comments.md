## Test environments
* local Mac install, R 4.5.3
* local fedora install, R 4.5.2
* rhub: R-devel and R-* on windows, linux, mac-arm64

## R CMD check results

No changes to worse with rev. depends clam coffee rbacon rplum (I maintain these packages)

The package now suggests p3k14c, a GitHub‑only package which provides optional archaeological radiocarbon datasets. All examples, tests and vignettes run without p3k14c installed. Its use is optional and guarded by availability checks.
