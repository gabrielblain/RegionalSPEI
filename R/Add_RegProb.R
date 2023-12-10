#' Add_RegProb
#'
#' @param dataset
#' A csv file with data from each sites.
#' The first column is the Years.
#' The second column is the Months (from 1 to 12).
#' The third column is the data coming from each site.
#' Use data(dataset) as example.
#' @param bestPar
#' Parameters of the regional distribution (selected by 'Add_Goodness').
#' Obtained from 'Add_RegPar'. See example.
#' @param distr
#' A string indicating the regional distribution
#' \dQuote{GEV} \dQuote{GLO} \dQuote{PE3}
#' \dQuote{GNO} \dQuote{GPA}.
#' @return
#' Cumulative probabilities for the sites forming
#' the homogeneous region or group.
#' @export
#' @importFrom lmom cdfgev cdfglo cdfpe3 cdfgno cdfgpa
#' @examples
#' data(dataset)
#' data(parGEV)
#' Add_RegProb(dataset,distr="GEV",bestPar=parGEV)

Add_RegProb=function(dataset,distr,bestPar){
  #if (!require(lmom)) install.packages('lmom')
  #library(lmom)
  n.sites=length(dataset[1,])-2
  Reg.prob=dataset
  if (n.sites<7){stop("The number of sites should be larger than 6.")}
  if (distr == "GEV"){
    for (month in 1:12){
      dataset.month=dataset[which(dataset[,2]==month),1:(n.sites+2)]
      probacum.month=dataset.month
      index.flood=colMeans(dataset.month,na.rm = TRUE)
      for (site in 3:(n.sites+2)){
        dataset.month[,site]=dataset.month[,site]-index.flood[site]
        probacum.month[,site]=as.matrix(cdfgev(dataset.month[,site],
                                               c(bestPar$loc[month],bestPar$scl[month],bestPar$shp[month])))
      }
      if (month==1){Reg.prob=probacum.month}else{
        Reg.prob=rbind(Reg.prob,probacum.month)}
    }
  } else if (distr == "GLO"){
    for (month in 1:12){
      dataset.month=dataset[which(dataset[,2]==month),1:(n.sites+2)]
      probacum.month=dataset.month
      index.flood=colMeans(dataset.month,na.rm = TRUE)
      for (site in 3:(n.sites+2)){
        dataset.month[,site]=dataset.month[,site]-index.flood[site]
        probacum.month[,site]=as.matrix(cdfglo(dataset.month[,site],
                                               c(bestPar$loc[month],bestPar$scl[month],bestPar$shp[month])))
      }
      if (month==1){Reg.prob=probacum.month}else{
        Reg.prob=rbind(Reg.prob,probacum.month)}
    }
  } else if (distr == "PE3"){
    for (month in 1:12){
      dataset.month=dataset[which(dataset[,2]==month),1:(n.sites+2)]
      probacum.month=dataset.month
      index.flood=colMeans(dataset.month,na.rm = TRUE)
      for (site in 3:(n.sites+2)){
        dataset.month[,site]=dataset.month[,site]-index.flood[site]
        probacum.month[,site]=as.matrix(cdfpe3(dataset.month[,site],
                                               c(bestPar$loc[month],bestPar$scl[month],bestPar$shp[month])))
      }
      if (month==1){Reg.prob=probacum.month}else{
        Reg.prob=rbind(Reg.prob,probacum.month)}
    }
  } else if (distr == "GNO"){
    for (month in 1:12){
      dataset.month=dataset[which(dataset[,2]==month),1:(n.sites+2)]
      probacum.month=dataset.month
      index.flood=colMeans(dataset.month,na.rm = TRUE)
      for (site in 3:(n.sites+2)){
        dataset.month[,site]=dataset.month[,site]-index.flood[site]
        probacum.month[,site]=as.matrix(cdfgno(dataset.month[,site],
                                               c(bestPar$loc[month],bestPar$scl[month],bestPar$shp[month])))
      }
      if (month==1){Reg.prob=probacum.month}else{
        Reg.prob=rbind(Reg.prob,probacum.month)}
    }
  } else if (distr == "GPA"){
    for (month in 1:12){
      dataset.month=dataset[which(dataset[,2]==month),1:(n.sites+2)]
      probacum.month=dataset.month
      index.flood=colMeans(dataset.month,na.rm = TRUE)
      for (site in 3:(n.sites+2)){
        dataset.month[,site]=dataset.month[,site]-index.flood[site]
        probacum.month[,site]=as.matrix(cdfgpa(dataset.month[,site],
                                               c(bestPar$loc[month],bestPar$scl[month],bestPar$shp[month])))
      }
      if (month==1){Reg.prob=probacum.month}else{
        Reg.prob=rbind(Reg.prob,probacum.month)}
    }
  } else{stop("parameter distr can only be GEV, GLO, PE3, GNO or GPA")}
  Reg.prob <- Reg.prob[order(Reg.prob[, 1]), ]
  return(Reg.prob)
}
