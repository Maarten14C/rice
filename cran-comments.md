## Test environments
* local Mac install, R 4.5.2
* local fedora install, R 4.5.2
* rhub (R-devel, release and patch): windows, linux, mac-arm

## R CMD check results

1.6.2: The test files are now updated and pass all tests. Apologies for the oversight.

1.6.1: The error on fedora-gcc related to the github package rnaturalearthhires has been resolved by no longer suggesting its installation using R commands. Instead, users are now directed to the relevant github page. As per CRAN request.

No notes, warnings or errors. No revdep changes to worse (coffee, rbacon, rplum, clam).
