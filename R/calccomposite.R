listcompositecalculatons <-
      list(
          "AVLT__total1_to_5" = function(AVLT__trial1, AVLT__trial2, AVLT__trial3, AVLT__trial4, AVLT__trial5, ...) {
            return(AVLT__trial1 + AVLT__trial2 + AVLT__trial3 + AVLT__trial4 + AVLT__trial5)
            },
          "TMT__interference" = function(TMTa, TMTb, ...) {
            return(TMT__b - TMT__a)
          }
          "BADS__ruleshift__1and2scorediff" = function( BADS__ruleshift__1score, BADS__ruleshift__2score){
            return(BADS__ruleshift__2score - BADS__ruleshift__1score)
          }
          "BADS__ruleshift__1and2timediff" = function( BADS__ruleshift__1time, BADS__ruleshift__2time){
            return(BADS__ruleshift__2time - BADS__ruleshift__1time)
          }
          )

#' Title Calculate composite score using lists of functions
#'
#' @param testname Name of the composite score
#' @param ... Arguments are substituent test scores that are needed in the calculation of the test score. May be named or unnamed (if named, names have to be correct)
#'
#' @return Value of the composite score
#' @export
#'
#' @examples calccomposite_func( "AVLT__total_1_to_5", AVLT_trial1 = 5, AVLT_trial2 = 5, AVLT_trial3 = 7, AVLT_trial4 = 8 , AVLT_trial5 = 12)
#' @examples calccomposite_func( "AVLT__total_1_to_5", 5, 5, 7, 8, 12)

calccomposite_func <- function( testname, ...){
  calccomposite[[testname]]( ...)
}



