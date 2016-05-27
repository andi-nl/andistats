##### These can all be preloaded in the sysdata.rda file I think #####
calccomposite <- function( compositename = "", v1, v2, v3, v4, v5){
  return( do.call( composite.list$func[composite.list$name == compositename], composite.list$args))
}

AVLTsum1to5.func <- function( v1, v2, v3, v4, v5, ...){ return( v1 + v2 + v3 + v4 + v5 )}
TMTint.func <- function( v1, v2, ...){ return( v2 - v1 )}

composite.list <- list( name = c("AVLTtotal1-5", "TMTinterference"),
                        func = c("AVLTsum1to5.func", "TMTint.func"),
                        args = list(as.name("v1"), as.name("v2"), as.name("v3"), as.name("v4"), as.name("v5")))
######################################################################

# These are example calls to OpenCPU that return the appropriate values
calccomposite( "AVLTtotal1-5", 5, 7, 8, 15, 15)
calccomposite( "TMTinterference", 60, 180)

# And in the CSV, we'll add a column named something like consistuents, or args, or whatever
# which contains AVLT_trial1,AVLT_trial2,AVLT_trial3,AVLT_trial4,AVLT_trial5 for the AVLTtotal1-5 row
# and contains TMT_a,TMT_b for the TMTinterference row.




