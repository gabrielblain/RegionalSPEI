#' Add_GoodnessFit
#'
#' @param dataset
#' A csv file with data from each sites.
#' The first column is the Years.
#' The second column is the Months (from 1 to 12).
#' The third column is the data coming from each site.
#' Use data(dataset) as example.
#' @param rho
#' Average spatial correlation (between 0 and 1).
#' defauld is 0.
#' @param Ns
#' Number of simulated groups of series.
#' Default is 100, but at least 500 is recommended.
#' @return
#' Hosking and Wallis (1991) Goodness of fit measure and
#' its bivatriate extension proposed in Kjeldsen and Prosdocimi (2015).
#' @export
#' @importFrom lmomRFA regsamlmu
#' @importFrom lmom pelwak pelgev pelglo pelpe3 pelgno pelgpa lmrgev lmrglo lmrpe3 lmrgno lmrgpa quawak
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm mahalanobis
#' @examples
#' data(dataset)
#' Add_GoodnessFit(dataset,rho=0.5,Ns=100)

Add_GoodnessFit <- function(dataset,rho=0,Ns=100){
  if (Ns<100){stop("Ns should be larger than 99.")}
  if (rho<0 || rho>1){stop("rho should Varies between 0 and 1.")}
  #if (!require(lmomRFA)) install.packages('lmomRFA')
  #if (!require(MASS)) install.packages('MASS')
  #library(lmomRFA)
  #library(MASS)
  n.sites <- length(dataset[1,])-2
  if (n.sites<7){stop("The number of sites should be larger than 6.")}
  numerador.bias <- matrix(NA,Ns,1)
  numerador.bias.quad <- matrix(NA,Ns,1)
  sknumerador.bias <- matrix(NA,Ns,1)
  sknumerador.bias.quad <- matrix(NA,Ns,1)
  kurt.sKew.numerador.bias <- matrix(NA,Ns,1)
  Z <- matrix(NA,12,5)
  Zbivar <- matrix(NA,12,5)
  rownames(Z) <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(Z) <- c("GEV","GLO","PE3","GNO","GPA")
  rownames(Zbivar) <- c("Jan","Feb","Mar","Apr","May","Jun",
                     "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(Zbivar) <- c("GEV","GLO","PE3","GNO","GPA")
  for (month in 1:12){
    dataset.month <- dataset[which(dataset[,2]==month),3:(n.sites+2)]
    index.flood <- colMeans(dataset.month,na.rm = TRUE)
    for (site in 1:n.sites){
      dataset.month[,site] <- dataset.month[,site]-index.flood[site]
    }
    x1.atoutset <- regsamlmu(dataset.month, lcv = FALSE)
    weight <- sum(x1.atoutset[,2])
    vet.l2 <- as.matrix(x1.atoutset[,2]*x1.atoutset[,4])
    vet.t3 <- as.matrix(x1.atoutset[,2]*x1.atoutset[,5])
    vet.t4 <- as.matrix(x1.atoutset[,2]*x1.atoutset[,6])
    vet.t5 <- as.matrix(x1.atoutset[,2]*x1.atoutset[,7])
    reg.l2 <- sum(vet.l2)/weight
    reg.t3 <- sum(vet.t3)/weight
    reg.t4 <- sum(vet.t4)/weight
    reg.t5 <- sum(vet.t5)/weight
    rmom <- c(0,reg.l2,reg.t3,reg.t4,reg.t5)
    reg.par <- try(pelwak(rmom))
    max.n.years <- max(x1.atoutset[,2])
    data.sim <- matrix(NA,max.n.years,n.sites)
    sigma <- matrix(rho,n.sites,n.sites); diag(sigma) <- 1
    for (ns in 1:Ns){
      u.sim <- pnorm(mvrnorm(n = max.n.years, mu=rep(0,n.sites), Sigma=sigma, tol = 1e-06, empirical = FALSE))
      for (site in 1:n.sites){
        data.sim[1:x1.atoutset[site,2],site] <- quawak(u.sim[1:x1.atoutset[site,2],site],
                                                    c(reg.par[1],reg.par[2],reg.par[3],reg.par[4],reg.par[5]))}
      x1.sim <- regsamlmu(data.sim, lcv = FALSE)
      vet.l2.sim<- as.matrix(x1.atoutset[,2]*x1.sim[,4])
      vet.t3.sim<- as.matrix(x1.atoutset[,2]*x1.sim[,5])
      vet.t4.sim<- as.matrix(x1.atoutset[,2]*x1.sim[,6])
      reg.l2.sim <- sum(vet.l2.sim)/weight
      reg.t3.sim <- sum(vet.t3.sim)/weight
      reg.t4.sim <- sum(vet.t4.sim)/weight
      rmom.sim <- c(0,reg.l2.sim,reg.t3.sim,reg.t4.sim)
      numerador.bias[ns] <- rmom.sim[4]-rmom[4]
      numerador.bias.quad[ns] <- (rmom.sim[4]-rmom[4])^2
      sknumerador.bias[ns] <- rmom[3]-rmom.sim[3]
      sknumerador.bias.quad[ns] <- (rmom[3]-rmom.sim[3])^2
      kurt.sKew.numerador.bias[ns] <- (rmom[3]-rmom.sim[3])*(rmom[4]-rmom.sim[4])
    }
    para.gev=pelgev(rmom)
    para.glo=pelglo(rmom)
    para.pe3=pelpe3(rmom)
    para.nor=pelgno(rmom)
    para.gp=pelgpa(rmom)
    B4=sum(numerador.bias)/Ns
    SD=sqrt((sum(numerador.bias.quad)-(Ns*(B4^2)))/(Ns-1))
    lmomgev=lmrgev(para = para.gev, nmom = 4)
    lmomglo=lmrglo(para = para.glo, nmom = 4)
    lmompe3=lmrpe3(para = para.pe3, nmom = 4)
    lmomnor=lmrgno(para = para.nor, nmom = 4)
    lmomgp=lmrgpa(para = para.gp, nmom = 4)
    Z.GEV=(lmomgev[4]-rmom[4]+B4)/SD
    Z.GLO=(lmomglo[4]-rmom[4]+B4)/SD
    Z.PE3=(lmompe3[4]-rmom[4]+B4)/SD
    Z.GNO=(lmomnor[4]-rmom[4]+B4)/SD
    Z.GPA=(lmomgp[4]-rmom[4]+B4)/SD
    Z[month,] <- as.matrix(cbind(Z.GEV,Z.GLO,Z.PE3,Z.GNO,Z.GPA))
    B3 <- sum(sknumerador.bias)/Ns
    SD.Skew <- sqrt(((sum(sknumerador.bias.quad))-Ns*(B3^2))/(Ns-1))
    SD.kurt.Skew <- (sum(kurt.sKew.numerador.bias)-(Ns*B3*B4))/(Ns-1)
    StR<-matrix(c(SD.Skew^2,SD.kurt.Skew,SD.kurt.Skew,SD^2),nrow=2,ncol=2)
    tb <- c(rmom[3]-B3,rmom[4]-B4)
    Tdist <- c(lmomgev[3],lmomgev[4])
    Zbivar.GEV <- mahalanobis(Tdist, tb, StR)
    Tdist <- c(lmomglo[3],lmomglo[4])
    Zbivar.GLO <- mahalanobis(Tdist, tb, StR)
    Tdist <- c(lmompe3[3],lmompe3[4])
    Zbivar.PE3 <- mahalanobis(Tdist, tb, StR)
    Tdist <- c(lmomnor[3],lmomnor[4])
    Zbivar.nor <- mahalanobis(Tdist, tb, StR)
    Tdist <- c(lmomgp[3],lmomgp[4])
    Zbivar.gp <- mahalanobis(Tdist, tb, StR)
    Zbivar[month,] <- as.matrix(cbind(Zbivar.GEV,Zbivar.GLO,Zbivar.PE3,Zbivar.nor,Zbivar.gp))
  }
  GoodnessFit <- list(Z=Z,Zbivar=Zbivar)
  return(GoodnessFit)
}
