#' #' Fast testing if a vector contains elements of a set
#' #'
#' #' @param x The vector which you want to test
#' #' @param table The values with which you want to see if these are in the set
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @importFrom fastmatch fmatch
#' #'
#' #' @examples
#' `%fin%` <- function(x, table) {
#'     # stopifnot(require(fastmatch))
#'     fmatch(x, table, nomatch = 0L) > 0L
#' }

require(fastmatch)
