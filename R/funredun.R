funredun<-function(spDat,funDat,method='gower',redund=TRUE,funDiv=FALSE){
  #sets progress bar
  pb<-utils::txtProgressBar(min=0,max=nrow(spDat),style=3,width=50,char="▓")
  
  #calculates distance between species based on functional traits
  spDist<-vegan::vegdist(funDat,method=method,upper=T)
  
  #creates a vector of Simpson's diversity for the community
  D<-as.vector(vegan::diversity(spDat,index='simpson'))
  
  #creates empty data frame P with same columns as input community data
  P<-spDat[0,]
  
  #loop to calculate proportion of each species at each site
  for(i in 1:nrow(spDat)){
    for(j in 1:ncol(spDat)){
      P[i,j]<-spDat[i,j]/sum(spDat[i,])
    }
  }
  
  #empty vectors of Q and FR
  Q<-as.vector(rep(0,nrow(P)))
  FR<-as.vector(rep(0,nrow(P)))
  
  #loop to calculate Q as product of dissimilarity between species i and j, proportion of species i, proportion of species j
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      for(k in 2:ncol(P)){
        Q[i]<-Q[i]+(usedist::dist_get(spDist,names(P)[j],names(P)[k])*P[i,j]*P[i,k])
      }
    }
    
    #calculates function redundance as D-Q or functional uniqueness Q/D
    ifelse(redund==TRUE,FR[i]<-(D[i]-Q[i]),ifelse(D[i]==0,FR[i]<-NA,FR[i]<-Q[i]/D[i]))
    utils::setTxtProgressBar(pb,i)
  }
  
  #builds output
  if (funDiv==FALSE) {
    FRoutput<-data.frame(FR)
    names(FRoutput)<-'Func_Redun'
    rownames(FRoutput)<-rownames(spDat)
    return(FRoutput)
  } else {
    FRoutput<-data.frame(FR,Q)
    names(FRoutput)<-c('Func_Redun','Func_Div')
    rownames(FRoutput)<-rownames(spDat)
    return(FRoutput)
  }
}
