#' Impute missing values by random draw from the marginal distributions
#'
#' This is a simple impute to create initial values for missing data
#'
#' @param dat a data matrix with possible missing values
#' @param miscode a set of values used for denoting missing values
#' @param nimpute the number of copies of the imputed data
#'
#'
#' @details This function imputes missing data by random draw from
#'               the marginal distribution of the respective variables
#' @return 1. Imputed data set with the last variable indexing the datasets
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  fit=mdimpute(dat=bone,miss=c(-9))
#'  fit[[1]]
#' }
#'
#'@export
#'

mdimpute=function(dat=dat,miscode=c(-9),nimpute=5){
  #1. Find the missing data indicators
  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
  } # missing data indicators

  #2. Get the conditional frequencies and impute using random draws.
  impdat=array(0,c(n*nimpute,p+1))
  for(ni in 1:nimpute){
    dstart=(ni-1)*n+1
    dend=ni*n
    impdat[dstart:dend,1:p]=as.matrix(dat)
    impdat[dstart:dend,p+1]=ni  # add a variable index the datasets
  }
  for(k in 1:p){
    obsset=subset(c(1:n),misdat[,k]==1) # observed data locations
    a=table(dat[obsset,k])
    a_cat=as.numeric(names(a))
    a_freq=as.numeric(a)
    if(sum(1-misdat[,k])>0){
      misset=subset(c(1:n),misdat[,k]==0) # subset the locations of missing values
      imp=sample(x=a_cat,prob=a_freq/sum(a_freq),size=length(misset)*nimpute,replace=TRUE)
      imp=matrix(imp, ncol=nimpute)
      for(ni in 1:nimpute){
        dstart=(ni-1)*n+1
        dend=ni*n
        impdat[(ni-1)*n+misset,k]=imp[,ni] # fill in the imputed values by random draws
      }
    }
  }

return(list(impdat))# add a data set index
}
