#' Perform multivariate and univariate normative comparisons
#'
#' @param myJSON hierarchically structured json file or json formatted text containing patient information with nested test information (including scores)
#' @importFrom jsonlite fromJSON toJSON
#' @return flat json file containing test information and test results (differences, statistics, degrees of freedom) with repeated results for the  multivariate test
#' @export
#'
#' @examples normcomp("jsonfile")
normcomp <- function( myJSON){

  json <- fromJSON(myJSON)
  no.patients <- nrow(json$patientScores)
  mypatdata <- NULL
  for( i in 1:(no.patients) ){
    demos <- c( json$patientScores$id[i], json$patientScores$age[i], json$patientScores$sex[i], json$patientScores$education[i])
    testinfo <- json$patientScores$test[i][[1]]
    no.tests <- nrow(testinfo)

    mypatdata <- rbind( mypatdata, cbind( matrix(rep(demos, no.tests), no.tests, byrow = T), testinfo))
  }

  colnames( mypatdata) <- c("patid", "age", "SEX", "EDU", "uniqueid", "label",
                            "Dataset", "SPSS_name", "highborder", "highweb", "lowborder", "lowweb",
                            "score")
  mypatdata[['patid']] <- as.character(mypatdata[['patid']])
  mypatdata[['SEX']] <- as.numeric(mypatdata[['SEX']])
  mypatdata[['age']] <- as.numeric(mypatdata[['age']]) - 65
  mypatdata[['EDU']] <- as.numeric(as.character(mypatdata[['EDU']]))
  mypatdata[['conf']] <- as.numeric(json$settings$conf)
  mypatdata[['sig']] <- json$settings$sig
  mypatdata[['normative']] <- json$settings$normative

  # defaultvalues
  uniqueID <- ANDImetadata[['uniqueid']]

  covariancemat <- betweencov + withincov
  rownames(covariancemat) <- uniqueID
  colnames(covariancemat) <- uniqueID

  totaloutputdataframe <- NULL
  whichtests <- unique(mypatdata[['uniqueid']])
  for( pat in unique(mypatdata[['patid']])){
    mydata <- mypatdata[mypatdata[['patid']] == pat,]
    mydata[['score']] <- ((mydata[['score']] ^ ANDImetadata[['mybestpowertransform']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]] *
                             sign(ANDImetadata[['mybestpowertransform']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]]) -
                             ANDImetadata[['mymean.transformedscores']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]]) /
                            ANDImetadata[['mysd.transformedscores']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]]) *
      ANDImetadata[['recode']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]]



    #C <- covariancemat[ rownames(covariancemat) %in% whichtests, colnames(covariancemat) %in% whichtests]
    C <- covariancemat[ rownames(covariancemat) %in% whichtests, colnames(covariancemat) %in% whichtests,drop=FALSE]

    P <- nrow(mydata)
    inv.C <- solve(C)

    #
    rownames(beta) <- rep(ANDImetadata[['uniqueid']],4)
    betaselection <- beta[rownames(beta) %in% whichtests]
    mydata[['pred']] <- (t( c(1, mydata[['SEX']][1], mydata[['age']][1], mydata[['EDU']][1])) %x% diag(1,P)) %*% betaselection

    tstatistics <- NULL
    pvalues <- NULL
    Tsquared <- NULL

    est.n <- ANDImetadata[['Nfinal']][ANDImetadata[['uniqueid']] %in% mydata[['uniqueid']]]
    min.est.n <- min(est.n)
    dfs <- est.n - 1
    g <- ( min.est.n  + 1 ) / min.est.n

    Tsquared <- ( 1 / g ) * ( ( min.est.n - P ) / ( ( min.est.n - 1 ) * P ) ) * t( mydata[['pred']] - mydata[['score']] ) %*% inv.C %*% ( mydata[['pred']] - mydata[['score']] )

    tstatistics <- ((mydata[['score']] - mydata[['pred']]) / ( sqrt(diag(C)) / sqrt(est.n))) * (1 / sqrt(est.n + 1))
    difference <- (mydata[['score']] - mydata[['pred']])

    tailed <- mydata[['sig']][1]
    if( tailed == "oneTailedLeft"){
      pvalues <- pt(tstatistics, dfs, lower = FALSE)
      MNCpvalue <- pf( Tsquared, P, min.est.n - P, lower = FALSE)
      if( sign(sum( mydata[['score']] - mydata[['pred']])) < 0){
        MNCpvalue <- pf( Tsquared, P, min.est.n - P, lower = FALSE) / 2
      } else if( sign(sum( difference)) >= 0){
        MNCpvalue <- 1
      }
    }

    if( tailed == "oneTailedRight"){
      pvalues <- pt(tstatistics, dfs, lower = TRUE)
      if( sign(sum( mydata[['score']] - mydata[['pred']])) > 0){
        MNCpvalue <- pf( Tsquared, P, min.est.n - P, lower = FALSE) / 2
      } else if( sign(sum( difference)) <= 0){
        MNCpvalue <- 1
      }
    }

    if( tailed == "twoTailed"){
      pvalues <- pt(abs(tstatistics), dfs, lower = FALSE) * 2
      MNCpvalue <- pf( Tsquared, P, min.est.n - P, lower = FALSE)
    }
    #
    # if( pvaluecorrection == "None"){
    #   pvalues[i,] <- pvalues[i,]
    # }
    # if( pvaluecorrection == "Bonferroni"){
    #   pvalues[i,] <- p.adjust(pvalues[i,], method = "bonferroni")
    # }
    # if( pvaluecorrection == "Holm"){
    #   pvalues[i,] <- p.adjust(pvalues[i,], method = "holm")
    # }



    if( tailed == "oneTailedLeft"){
      inneredge <- qt( (1 - ( mydata[['conf']][1]  / 100)), dfs, lower.tail = T)
      outeredge <- rep(99999,P)
    }

    if( tailed == "oneTailedRight"){
      inneredge <- rep(-99999,P)
      outeredge <- qt( (1 - ( mydata[['conf']][1]  / 100)), dfs, lower.tail = F)
    }

    if( tailed == "twoTailed"){
      inneredge <- qt( (1 - ( mydata[['conf']][1] / 100)) / 2, dfs)
      outeredge <- abs( qt( (1 - ( mydata[['conf']][1] / 100)) / 2, dfs))
    }

    tstatistics <- round(tstatistics, 2)


    pvalues <- format(round(pvalues, 3), nsmall = 3)
    MNCpvalue <- format(round(MNCpvalue, 3), nsmall = 3)


    longtestnames <- paste(ANDImetadata[['Long.name.1']], ANDImetadata[['Long.name.2']], ANDImetadata[['Long.name.3']][!is.na(ANDImetadata[['Long.name.3']])])
    plotnames <- trimws(paste(ANDImetadata[['ID1']], ANDImetadata[['Long.name.2']], ANDImetadata[['Long.name.3']][!is.na(ANDImetadata[['Long.name.3']])]))
    shortestpossiblenames <- trimws(paste(ANDImetadata[['ID1']], ANDImetadata[['ID2']], ANDImetadata[['ID3']][!is.na(ANDImetadata[['Long.name.3']])]))

    whichnames <- c()
    for( eachvar in mydata[['uniqueid']]){
      whichnames <- c( whichnames, which( ANDImetadata[['uniqueid']] == eachvar))
    }

    myoutputdataframe <- data.frame(
      id = mydata[['patid']],
      testname = mydata[['uniqueid']],
      longtestname = longtestnames[whichnames],
      plotname = plotnames[whichnames],
      shortestname = shortestpossiblenames[whichnames],
      tails = mydata[['sig']],
      inneredge = inneredge,
      outeredge = outeredge,
      univariatedifferences = difference,
      univariateT = tstatistics,
      univariatedf = dfs,
      univariatep = pvalues,
      multivariatedifference = sum(difference),
      multivariateT = rep(Tsquared,each = P),
      multivariatedf = rep(paste0(P, ", ", min.est.n - P),each = P),
      multivariatep = MNCpvalue
    )
    totaloutputdataframe <- rbind( totaloutputdataframe, myoutputdataframe)
  }
  myoutputdata <- toJSON( totaloutputdataframe,pretty = T)
  #cat(myoutputdata)
  return(myoutputdata)

}
