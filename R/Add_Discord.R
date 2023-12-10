#' Add_Discord
#'
#' @param dataset
#' A csv file with data from each sites.
#' The first column is the Years.
#' The second column is the Months (from 1 to 12).
#' The third column is the data coming from each site.
#' Use data(dataset) as example.
#' @return
#' Hosking and Wallis' discordancy measure
#' calculated using the additive approach as proposed in
#' MARTINS et al. (2022) https://doi.org/10.1590/1678-4499.20220061
#' @export
#' @importFrom lmomRFA regsamlmu
#' @importFrom rrcov   getDistance Cov
#' @examples
#' data(dataset)
#' Add_Discord(dataset)

Add_Discord <- function(dataset){
  #if (!require(lmomRFA)) install.packages('lmomRFA')
  #if (!require(rrcov)) install.packages('rrcov')
  #library(lmomRFA)
  #library(rrcov)
  n <- length(dataset[1,])
  n.sites <- n-2
  d <- matrix(NA,12,(n.sites)+1)
  for (month in 1:12){
    dataset.month <- dataset[which(dataset[,2]==month),3:n]
    x1.atoutset <- regsamlmu(dataset.month, lcv  =  FALSE)
    d[month,] <- c(month,sqrt(getDistance(Cov(x1.atoutset[,4:6]))))}
  sites <- c(colnames(dataset.month))
  colnames(d) <- c("Month",sites)
  return(d)
}
