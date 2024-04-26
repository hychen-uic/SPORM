#' The function finds the locations of setb elements in order of seta.
#' Both sequences are sorted from smallest to the largest
#'
#' @param seta a set of distinctive elements ordered from smallest to the largest.
#' @param setb a set of distinctive elements ordered from smallest to the largest.
#' @param method either use "nearestleft" for low bound of interval censoring
#'        or "nearestright" for upp bound of interval censoring
#'
#' @details This function identifies locations of setb in the sequence of seta
#'
#' @return a set of locations in sequence seta for the positions in setb
#'
#' @examples \dontrun{
#' seta=c(-1,0,3,5,12,17,18)
#' setb=c(-9,-7,-5,2,2,3,4,5,8,15,17,17,23,27,29)
#' insertloc(seta,setb)
#' insertloc(seta,setb,method="nearestright")
#' }
#'
#'@export
#'

insertloc=function(seta,setb,method='nearestleft'){
  # seta is a set of ordered values
  # setb is another set of ordered values
  # the function is to locate the nearest seta values for every value in setb

  loc=rep(0,length(setb))
  if(method=="nearestleft"){
    i=1
    k=1
    while(k<=length(setb)){
      if(setb[k]>max(seta)){
        loc[k]=length(seta)
        k=k+1
      }else{
        if(setb[k]<=seta[i]){
          loc[k]=i
          k=k+1
        }else{
          i=i+1
        }
      }
    }
  }else{
    i=length(seta)
    k=length(setb)
    while(k>=1){
      if(setb[k]>=max(seta)){
        loc[k]=length(seta)
        k=k-1
      }else{
        if(setb[k]>=seta[i]){
          loc[k]=i
          k=k-1
        }else{
          i=i-1
          if(i==0){loc[1:k]=1;k=0}
        }
      }
    }
  }

  list(loc)
}
