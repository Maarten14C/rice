#' shroud Data
#'
#' A dataset containing the radiocarbon dates on the Shroud of Turin, from three labs
#'
#' @docType data
#' @format A data frame with 1968 rows and 15 variables.
#' \describe{
#'   \item{ID}{Lab numbers. Replicates are indicates with .1, .2, etc.}
#'   \item{y}{Radiocarbon year}
#'   \item{er}{Lab error}
#' }
#' @source Data taken from Damon et al. 1989 [Nature] <doi:10.1038/337611a0>, see also Christen 1994 [Applied Statistics] <doi:10.2307/2986273>
#' @examples
#' data(shroud)
#' head(shroud)
"shroud"
