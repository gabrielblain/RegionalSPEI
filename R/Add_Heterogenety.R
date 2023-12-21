#' Add_Heterogenety
#'
#' @param dataset
#' A csv file with data from each sites.
#' The first column is the Years.
#' The second column is the Months (from 1 to 12).
#' The third column is the data coming from each site.
#' Use data(dataset) as example.
#' @param rho
#' A 1-column csv file or matrix
#' with 12 monthly values for
#' coefficients of spatial correlation.
#' The first coefficient is for January.
#' The last coefficient is for December.
#' @param Ns
#' Number of simulated groups of series.
#' Default is 100, but at least 500 is recommended.
#' @return
#' 'Hosking' and 'Wallis'' heterogeneity measure
#' calculated using the additive approach proposed in
#' 'LOPES MARTINS et al. (2022)'
#' https://doi.org/10.1590/1678-4499.20220061
#' @export
#' @importFrom lmomRFA regsamlmu
#' @importFrom lmom pelkap pelwak quakap quawak
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm sd
#' @examples
#' data(dataset)
#' rho=as.matrix(rep(0.5,12))
#' Add_Heterogenety(dataset,rho = rho,Ns = 100)

Add_Heterogenety <- function(dataset,rho,Ns = 100){
  if (Ns<100){stop("Ns should be larger than 99.")}
  max.rho <- length(which(rho < 1))
  min.rho <- length(which(rho > (-1)))
  if (max.rho !=12 || min.rho !=12){stop("rho should largaer than -1 and smaller than 1.")}
  n.sites <- length(dataset[1,])-2
  if (n.sites<7){stop("The number of sites should be larger than 6.")}
  vetor.numerador <- matrix(NA,n.sites,1)
  V.sim <- matrix(NA,Ns,1)
  H.month <- matrix(NA,12,1)
  for (month in 1:12){
    rho.month <- rho[month,1]
    dataset.month <- dataset[which(dataset[,2]==month),3:(n.sites+2)]
    index.flood <- colMeans(dataset.month,na.rm = TRUE)
    for (site in 1:n.sites){
      dataset.month[,site] <- dataset.month[,site]-index.flood[site]
    }
    x1.atoutset <- regsamlmu(dataset.month, lcv = FALSE)
    weight <- sum(x1.atoutset[,2])
    vet.l2<- as.matrix(x1.atoutset[,2]*x1.atoutset[,4])
    vet.t3<- as.matrix(x1.atoutset[,2]*x1.atoutset[,5])
    vet.t4<- as.matrix(x1.atoutset[,2]*x1.atoutset[,6])
    vet.t5<- as.matrix(x1.atoutset[,2]*x1.atoutset[,7])
    reg.l2 <- sum(vet.l2)/weight
    reg.t3 <- sum(vet.t3)/weight
    reg.t4 <- sum(vet.t4)/weight
    reg.t5 <- sum(vet.t5)/weight
    rmom <- c(0,reg.l2,reg.t3,reg.t4,reg.t5)
    reg.par <- try(pelkap(rmom),TRUE)
    is.kappa=length(reg.par)
    if (is.kappa==1){reg.par=try(pelwak(rmom),TRUE)}
    for (v in 1:n.sites){
      vetor.numerador[v] <- x1.atoutset[v,2]*(x1.atoutset[v,4]-reg.l2)^2}
    V <- sqrt(sum(vetor.numerador)/weight)
    max.n.years <- max(x1.atoutset[,2])
    data.sim <- matrix(NA,max.n.years,n.sites)
    sigma <- matrix(rho.month,n.sites,n.sites); diag(sigma) <- 1
    V.sim <- matrix(NA,Ns,1)
    for (ns in 1:Ns){
      u.sim <- pnorm(mvrnorm(n  =  max.n.years, mu = rep(0,n.sites), Sigma = sigma, tol  =  1e-06, empirical  =  FALSE))
      for (site in 1:n.sites){
        if (is.kappa==1){
        data.sim[1:x1.atoutset[site,2],site] <- quawak(u.sim[1:x1.atoutset[site,2],site],
                                                    c(reg.par[1],reg.par[2],reg.par[3],reg.par[4],reg.par[5]))}else{
        data.sim[1:x1.atoutset[site,2],site] <- quakap(u.sim[1:x1.atoutset[site,2],site],
                                                    c(reg.par[1],reg.par[2],reg.par[3],reg.par[4]))}
                             }
      x1.sim <- regsamlmu(data.sim, lcv  =  FALSE)
      vet.l2.sim<- as.matrix(x1.atoutset[,2]*x1.sim[,4])
      reg.l2.sim <- sum(vet.l2.sim)/weight
      for (v in 1:n.sites){
        vetor.numerador[v] <- x1.atoutset[v,2]*(x1.sim[v,4]-reg.l2.sim)^2}
      V.sim[ns,1] <- sqrt(sum(vetor.numerador)/weight)
    }
    H.month[month,1] <- c((V-mean(V.sim))/sd(V.sim))
  }
  return(H.month)
}
