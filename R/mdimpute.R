#' Impute missing values by random draw from the marginal distributions
#'
#' This is a simple impute to create initial values for missing data
#'
#' @param dat a data matrix with possible missing values
#' @param miscode a set of values used for denoting missing values
#'
#'
#' @details This function imputes missing data by random draw from
#'               the marginal distribution of the respective variables
#' @return 1. One imputed data set
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  mdimpute(dat=bone,miss=c(-9))
#' }
#'
#'@export
#'

mdimpute=function(dat=dat,miscode=c(-9)){
  #1. Find the missing data indicators
  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
    } # missing data indicators

  #2. Get the conditional frequencies and impute using random draws.
  impdat=dat
  for(k in 1:p){
    a=table(dat[,k])
    a_cat=as.numeric(names(a))
    a_freq=as.numeric(a)
    if(sum(1-misdat[,k])>0){
      len=length(a_cat)
      imp=sample(x=a_cat[2:len],prob=a_freq[2:len]/sum(a_freq[2:len]),size=a_freq[1],replace=TRUE)
      misset=order(misdat[,k])[1:a_freq[1]] # subset the locations of missing values
      impdat[misset,k]=imp # fill in the imputed values by random draws
      }
#    print(table(dat[,k]))
#    print(table(impdat[,k]))
   }
return(list(impdat))
}
