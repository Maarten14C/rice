## Test environments
* local Mac install, R 4.4.1
* local fedora 41 install, R 4.4.1
* win-builder (devel)
* rhub (R-devel): macos-arm64, macos, windows, linux, gcc14

## R CMD check results
There were no ERRORs or WARNINGs.

The NOTE about UTF-8 strings relates to a .csv file with URF-8 encoding for international names and places. 

The NOTE about non-standard things in the check directory (‘rnaturalearthhires’) has been resolved by adding ^rnaturalearthhires$ to the .Rbuildignore file.
