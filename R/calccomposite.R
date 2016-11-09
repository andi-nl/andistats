#' Title Calculate composite score using one of many functions
#'
#' @param inputfile JSON file that contains compVar, which is the name of the composite variable, and args, which contains the named arguments and values that go into the calculations
#'
#' @return JSON formatted compositescore
#' @importFrom jsonlite fromJSON toJSON
#' @export
#'
#' @examples calccomposite( "exampleJSON.json")

calccomposite <- function( inputfile){
  inputvals <- inputfile
  compositescore <- as.data.frame(do.call( compositefunctions[[inputvals$compVar]], as.list(as.numeric(inputvals$args)) ))
  names(compositescore) <- inputvals$compVar
  return(toJSON(compositescore))
  }
