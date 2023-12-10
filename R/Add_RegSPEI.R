#' Add_RegSPEI
#'
#' @param RegProb
#' The cumulative probabilities for the sites forming
#' the homogeneous region or group.
#' Calculated using the 'Add_RegProb' function.
#' @return
#' The Regional SPEI
#' @export
#' @importFrom stats qnorm
#' @examples
#' data(RegProb)
#' RegSPEI=Add_RegSPEI(RegProb)

Add_RegSPEI=function(RegProb){
  n.sites=length(RegProb[1,])-2
  if (n.sites<7){stop("The number of sites should be larger than 6.")}
  RegSPEI = RegProb
  for (site in 3:(n.sites+2)){
    RegSPEI[,site]=qnorm(RegProb[,site], mean = 0, sd = 1)
  }
  return(RegSPEI)
}
