#' Add_RegPar
#'
#' @param dataset
#' A csv file with data from each sites.
#' The first column is the Years.
#' The second column is the Months (from 1 to 12).
#' The third column is the data coming from each site.
#' See 'Dataset.csv' as example.
#' @return
#' Regional parameter for several distributions
#' GEV: Generalized extreme value distribution.
#' GLO: Generalized logistic distribution.
#' PE3: Pearson type III distribution.
#' GNO: Generalized normal distribution.
#' GPA: Generalized pareto distribution.
#' @export
#' @importFrom lmomRFA regsamlmu
#' @importFrom lmom pelgev pelglo pelpe3 pelgno pelgpa
#' @examples
#' data(dataset)
#' Add_RegPar(dataset)


Add_RegPar=function(dataset){
  #if (!require(lmomRFA)) install.packages('lmomRFA')
  #library(lmomRFA)
  n.sites=length(dataset[1,])-2
  if (n.sites<7){stop("The number of sites should be larger than 6.")}
  par.GEV=matrix(NA,12,3)
  rownames(par.GEV)=c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(par.GEV)=c("loc","scl","shp")
  par.GLO=matrix(NA,12,3)
  rownames(par.GLO)=c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(par.GLO)=c("loc","scl","shp")
  par.PE3=matrix(NA,12,3)
  rownames(par.PE3)=c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(par.PE3)=c("loc","scl","shp")
  par.GNO=matrix(NA,12,3)
  rownames(par.GNO)=c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(par.GNO)=c("loc","scl","shp")
  par.GPA=matrix(NA,12,3)
  rownames(par.GPA)=c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Ago","Sep","Oct","Nov","Dec")
  colnames(par.GPA)=c("loc","scl","shp")

  for (month in 1:12){
    dataset.month=dataset[which(dataset[,2]==month),3:(n.sites+2)]
    index.flood=colMeans(dataset.month,na.rm = TRUE)
    for (site in 1:n.sites){
      dataset.month[,site]=dataset.month[,site]-index.flood[site]
    }
    x1.atoutset=regsamlmu(dataset.month, lcv = FALSE)
    weight=sum(x1.atoutset[,2])
    vet.l2<- as.matrix(x1.atoutset[,2]*x1.atoutset[,4])
    vet.t3<- as.matrix(x1.atoutset[,2]*x1.atoutset[,5])
    vet.t4<- as.matrix(x1.atoutset[,2]*x1.atoutset[,6])
    reg.l2 <- sum(vet.l2)/weight
    reg.t3 <- sum(vet.t3)/weight
    reg.t4 <- sum(vet.t4)/weight
    rmom=c(0,reg.l2,reg.t3,reg.t4)
    par.GEV[month,]=pelgev(rmom)
    par.GLO[month,]=pelglo(rmom)
    par.PE3[month,]=pelpe3(rmom)
    par.GNO[month,]=pelgno(rmom)
    par.GPA[month,]=pelgpa(rmom)
  }
  RegPar=list(par.GEV=par.GEV,par.GLO=par.GLO,par.PE3=par.PE3,
              par.GNO=par.GNO,par.GPA=par.GPA)
  return(RegPar)
}
