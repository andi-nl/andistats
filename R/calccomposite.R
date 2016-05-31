##### These can all be preloaded in the sysdata.rda file I think #####

calccomposite <-
      list(
          "AVLTtotal1-5" = function(AVLT_trial1, AVLT_trial2, AVLT_trial3, AVLT_trial4, AVLT_trial5, ...) {
            return(AVLT_trial1 + AVLT_trial2 + AVLT_trial3 + AVLT_trial4 + AVLT_trial5)
            },
          "TMTinterference" = function(TMTa, TMTb, ...) {
            return(TMTb - TMTa)
            }
          )

# These are example calls to OpenCPU that return the appropriate values
calccomposite[["AVLTtotal1-5"]]( 5, 7, 8, 15, 15)
calccomposite[["TMTinterference"]]( 60, 180)



# This might already be enough. But, there may be advantages to calling a function.
calccomposite_func <- function( testname, ...){
  calccomposite[[testname]]( ...)
}

# This works correctly (returns 37)
calccomposite_func( "AVLTtotal1-5", AVLT_trial1 = 5, AVLT_trial2 = 5, AVLT_trial3 = 7, AVLT_trial4 = 8 , AVLT_trial5 = 12)
calccomposite_func( "AVLTtotal1-5", 5, 5, 7, 8, 12)

# This correctly doesn't work (returns argument "TMTb" is missing)
calccomposite_func( "TMTinterference", AVLT_trial1 = 5, AVLT_trial2 = 5, AVLT_trial3 = 7, AVLT_trial4 = 8 , AVLT_trial5 = 12)

# This works, but this might indicate  that it's too flexible
calccomposite_func( "TMTinterference", AVLT_trial1 = 5, 60, AVLT_trial2 = 5, AVLT_trial3 = 7, AVLT_trial4 = 8 , AVLT_trial5 = 12, 180)

# And in the CSV, we'll add a column named something like constituents, or args, or whatever
# which contains AVLT_trial1,AVLT_trial2,AVLT_trial3,AVLT_trial4,AVLT_trial5 for the AVLTtotal1-5 row
# and contains TMT_a,TMT_b for the TMTinterference row.




