#' Functional Redundancy
#' 
#' The function calculates functional redundancy as Simpson's D minus dissimilarity between species weighted by the proportions of the species.
#'
#' @param spDat Data frame with rows as sites, columns as species, and elements as counts
#' @param funDat Data frame with rows as species (same as spDat column names), columns as functional traits, elements as counts, measures, binary, etc.
#' @param method Default is Bray-Curtis dissimilarity. Available options include "bray", "gower", and "altGower". See \code{\link[vegan]{vegdist}} for details. 
#' @return A data frame with rows as sites and a column of functional redundancy
#' @export

funredun=function(spDat,funDat,method='bray'){

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
  
  #sets progress bar
  pb=txtProgressBar(min=0,max=nrow(P),style=3,width=50,char="=")
  
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
    
    #calculates function diversity as Simpson's D - Q
    FR[i]=D[i]-Q[i]
    setTxtProgressBar(pb,i)
  }
  FRoutput=data.frame(cbind(FR))
  names(FRoutput)='Func Redun'
  rownames(FRoutput)=rownames(spDat)
  return(FRoutput)
}