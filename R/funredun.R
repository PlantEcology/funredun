#' Functional Redundancy
#'
#' The function calculates functional redundancy as Simpson's D minus dissimilarity between species weighted by the proportions of the species.
#'
#' @param spDat Data frame with rows as sites, columns as species, and elements as counts
#' @param funDat Data frame with rows as species (same as spDat column names), columns as functional traits, elements as counts, measures, binary, etc.
#' @param method Available options include "bray", "gower", and "altGower". See \code{\link[vegan]{vegdist}} for details. Default is Gower distance.
#' @param redund Redundancy calculation as difference from Simpson's D (R = D - Q) or uniqueness (U = Q/D). Default is difference (TRUE).
#' @param funDiv Functional Diversity as Rao's Q (Botta-Dukát 2005). Default is false.
#' @return A data frame with rows as sites and a column of functional redundancy
#' @export

funredun=function(spDat,funDat,method='gower',redund=TRUE,funDiv=FALSE){
  #sets progress bar
  pb=utils::txtProgressBar(min=0,max=nrow(spDat),style=3,width=50,char="=")

  #calculates distance between species based on functional traits
  spDist=vegan::vegdist(funDat,method=method,upper=T)

  #creates a vector of Simpson's diversity for the community
  D=as.vector(vegan::diversity(spDat,index='simpson'))

  #creates empty data frame P with same columns as input community data
  P=spDat[0,]

  #loop to calculate proportion of each species at each site
  for(i in 1:nrow(spDat)){
    for(j in 1:ncol(spDat)){
      P[i,j]=spDat[i,j]/sum(spDat[i,])
    }
  }

  #empty vectors of Q and FR
  Q=vector()
  FR=vector()
  Q=rep(0,nrow(P))
  FR=rep(0,nrow(P))

  #loop to calculate Q as product of dissimilarity between species i and j, proportion of species i, proportion of species j
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      for(k in 2:ncol(P)){
        Q[i]=Q[i]+(usedist::dist_get(spDist,names(P)[j],names(P)[k])*P[i,j]*P[i,k])
      }
    }

    #calculates function redundance as D-Q or functional uniqueness Q/D
    if (redund==TRUE) {
      FR[i]=(D[i]-Q[i])
    } else {
      if (D[i]==0) {
        FR[i]=NA
      } else {
      FR[i]=Q[i]/D[i]
      }
    }
    utils::setTxtProgressBar(pb,i)
  }

  #builds output
  if (funDiv==FALSE) {
    FRoutput=data.frame(cbind(FR))
    names(FRoutput)='Func Redun'
    rownames(FRoutput)=rownames(spDat)
    return(FRoutput)
  } else {
    FRoutput=data.frame(cbind(FR,Q))
    names(FRoutput)=c('Func Redun','Func Div')
    rownames(FRoutput)=rownames(spDat)
    return(FRoutput)
  }
}
