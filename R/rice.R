# to do:

# add functions to do with marine dR
# uncertainties in dR and diet?
# smoothing


# internal functions to speed up reading and writing files, using the data.table R package if present
fastread <- function(fl, ...)
  if("data.frame" %in% (.packages())) # some Macs have problems with this package
    as.data.frame(data.table::fread(fl), ...) else
      read.table(fl, ...)



fastwrite <- function(fl, ...)
  if("data.frame" %in% (.packages())) # some Macs have problems with this package
    data.table::fwrite(as.data.frame(fl), ...) else
      write.table(fl, ...)


