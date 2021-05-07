#process network features
library(ggplot2)
library(vegan)
library(ineq)

f_01_averageLocalClusteringCoeff <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    #return(mean(transitivity(gc, type = 'local'), na.rm = T)) 
    ##there seems to be a bug in igraph, when the graph is directed, it doesn't calculate
    ##transitivity as it described that for directed graph correctly
    return(mean(transitivity(as.undirected(gc,'collapse'), type = 'local'),na.rm=T))
    
  } else {
    return(mean(transitivity(as.undirected(g,'collapse'), type = 'local'),na.rm=T))
  }
}

f_02_averageLocalClusteringCoeffWeighted <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    #return(mean(transitivity(gc, type = 'barrat'), na.rm = T))
    ##barrat transitivity does not work for graphs with multiple and/or loop edges.
    ##for a directed graph, call as.undirected with the collapse mode first.
    return(mean(transitivity(as.undirected(gc,'collapse'), type = 'barrat'),na.rm=T))
  } else {
    #return(mean(transitivity(g, type = 'barrat'), na.rm = T))
    return(mean(transitivity(as.undirected(g,'collapse'), type = 'barrat'),na.rm=T))
  }
}


f_03_globalClusteringCoeff <- function(g, GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(transitivity(gc, type = 'global'))
  } else {
    return(transitivity(g, type = 'global'))
  }
}

f_04_gdensity <- function(g, directed=T){
  if (directed) {
    ne<-length(E(g))
    nv<-length(V(g))
    return(ne/(nv*(nv-1)))
  }else{
    return(graph.density(g))
  }
}

f_05_proportionSingletons <- function(g){
  sum(degree(g)==0)/length(V(g))
}

##useful in directed networks, as selection case has more reciprocal connections
f_06_proportionEndpoints <- function(g){
  sum(degree(g)==1)/length(V(g))
}

f_07_meanDegree <- function(g){
  mean(degree(g))
}

#
#f_08_meanDegreeNotSingletons <- function(g){
#  mean(degree(g)[degree(g)!=0])
#}

f_08_assortativityDegree<-function(g) {
  assortativity_degree(g, directed = TRUE)
}

f_09_meanStrength <- function(g){
  mean(strength(g))
}

#test the log-log pearson coefficient of degree distribution and degree
f_10_straightness <- function(g){
  a<-degree_distribution(g)
  a<-log(a[2:length(a)])
  b<-log(1:length(a))
  cor(a[a!=-Inf],b[a!=-Inf])
}
#
#f_10_meanStrengthNotSingletons <- function(g){
#  mean(strength(g)[degree(g)!=0])
#}

##entropy
f_11_entropyDegreeDistribution <-  function(g,verbose=F){
  y=degree(g)
  freq=prop.table(table(y))
  if (verbose){print(freq)}
  -sum(freq * log(freq))
}


f_12_ratioComponents <- function(g) { 
  cl <- clusters(g) 
  cl$no/length(V(g))
}

giant.component <- function(g) { 
  cl <- clusters(g) 
  induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
}

#this doesn't provide any new information as it's the reciprocal of f_12
#f_13_averageComponentSize <- function(g) { 
#  cl <- clusters(g) 
#  mean(cl$csize)
#}

###change to measure the ratio between bigger and the second largest components
f_13_giantConnectedRatio <- function(g){
  cl <- clusters(g)
  if (cl$no>1){
    d<-sort(cl$csize,decreasing=T)
    return(d[1]/d[2])
  }else{
    return(0)
  }
}

###change this to evenness
#f_14_entropyComponentSize <-  function(g,verbose=F){
#  cl <- clusters(g)
#  y=cl$csize
#  freq=prop.table(table(y))
#  if (verbose){print(freq)}
#  -sum(freq * log(freq, base = 2))
#}

f_14_evennessComponentSize <-  function(g){
  cl <- clusters(g)$csize
  if (length(cl)>1){
    Gini(cl, corr=T, na.rm=T)
  }else{
    return(1)
  }
  
}

###added 
f_15_CentralPointDominance<-function(g){
  gc <- giant.component(g)
  N<-length(V(gc))
  B<-betweenness(gc, weights=NA,normalized=T)
  1/(N-1)*(N*max(B)-sum(B))
}


f_16_meanEccentricity <- function(g){
  mean(eccentricity(g))
}

f_17_gdiameter <- function(g){
  return(diameter(g,weights=rep(1,length(E(g)))))
}

###slightly changed to include only components larger than 1
f_18_meanDiameterComponents <- function(g){
  E(g)$weight <- 1
  cl <- clusters(g)
  d <- 0
  n<-0
  for (m in 1:(cl$no)){
    if (cl$csize[m]>1) {
      d <- d + diameter(induced.subgraph(g, which(cl$membership == m)))
      n<-n+1
    }
  }
  return(d/n)
}

f_19_globalEfficiency <- function(g){
  ###Latora and Marchiori, use reciprocal of the global efficiency, eq. 15
  ###because this does not present divergence problem as in using average geodesic distance
  d_ij <- shortest.paths(g,weights=rep(1,length(E(g))))
  d_ij <- d_ij[lower.tri(d_ij)]
  d_ij <- d_ij[!is.infinite(d_ij)]
  N=length(V(g))
  (N*(N-1))*(1/sum(1/d_ij))###changed to the reciprocal of global efficiency
}

f_20_averageClosenessCentrality <- function(g){
  ###doesn't make sense for not connected graphs, therefore, calculating only the giant
  ###component
  gc <- giant.component(g)
  mean(closeness(gc, weights = NA))
}


motifsProportion <- function(g){ # Calculate the proportion of each of the 16 motifs out of the total motifs found
  motifs <- graph.motifs(g, size = 3)
  motifs.prop <- motifs/sum(motifs, na.rm = T)
  #names(motifs.prop) <- paste('motif',1:16,sep='')
  return(motifs.prop[-c(1,2,4,12)]) # Motifs 1,2,4 are constantly (across all cutoffs) NA and 12 is constantly 0 (I tested it)
}


##reciprocity calculates how in and out degrees reciprocate
f_33_reciprocity<-function(g){
  return(reciprocity(g))
}

##correlation between in and out number of links
f_34_inoutCorrelation<-function(g){
  return(cor(degree(g,mode='in'),degree(g,mode='out')))
}

####modules, module evenness, and module FST
##module FST calculation
communityDistances<-function(g,ceb, dis = T){
  mat<-get.adjacency(g,attr="weight") 
  if (dis == F) {
    #note that here, weight should be the distance between edges
    #if input g is similarity matrix, then this should be reverted
    mat<-1-mat
  }
  memberMap<-membership(ceb)
  comNumber<-length(ceb)
  if (comNumber<2){
    print("no sub communities detected")
    return(matrix(NA,comNumber,comNumber))
  }
  comSize<-sizes(ceb)
  withinDiv<-c()
  outMat<-matrix(NA,comNumber,comNumber)
  indList<-list()
  #first calculate all within diversities
  for (i in 1:comNumber) {
    indOut<-which(memberMap==i,arr.ind=T)
    indList[[i]]<-indOut
    if (comSize[i]>1) {
      focalMat<-mat[indOut,indOut]
      piWithin = sum(focalMat[upper.tri(focalMat)])+
        sum(focalMat[lower.tri(focalMat)])
      withinDiv<-c(withinDiv,piWithin)
    }else{
      withinDiv<-c(withinDiv,NA)
    }
  }
  ##then calculate all between diversities
  for (i in 1:(comNumber-1)) {
    if (comSize[i]==1) {
      next
    }
    for (j in (i+1):comNumber) {
      if (comSize[j]==1) {
        next
      }
      
      piBetween = (sum(mat[indList[[i]],indList[[j]]])+sum(
        mat[indList[[j]],indList[[i]]]))/2/comSize[i]/comSize[j]
      outMat[i,j]<-1-(withinDiv[i]+withinDiv[j])/
        (comSize[i]*(comSize[i]-1)+comSize[j]*(comSize[j]-1))/piBetween
      
    }
  }
  return(outMat)
}

calCommunity<-function(go,useOriginal, cutoff=0.99, similarity=T){
  if(useOriginal){
    g<-go
  }else{
    if (similarity){
      g<-delete.edges(go,which(E(go)$weight<cutoff))    
    }else{
      g<-delete.edges(go,which(E(go)$weight<quantile(E(go)$weight,cutoff)))    
    }
  }
  
  E(g)$weight<-(1-E(g)$weight)*99+1
  ceb<-cluster_edge_betweenness(g)
  outMat<-communityDistances(go,ceb, F)
  #hist(outMat,40)
  #image(outMat,col=heat.colors(20),legend=T)  
  #plotMember(g,ceb)
  return(list(ceb,outMat))
}

##community size evenness, gini index, FST###

f_35_communitySizeEvenness<-function(ceb){
  if (length(ceb)>1) {
    return(Gini(sizes(ceb),corr=T, na.rm=T))
  }else{
    return(1)
  }
}

f_36_numberCommonCommunities<-function(ceb){
  return(length(sizes(ceb)[sizes(ceb)>1]))
}

f_37_CommunityRatio<-function(ceb){
  if (length(ceb)>1) {
    d<-sort(sizes(ceb),decreasing=T)
    return(d[1]/d[2])
  }else{
    return(0)
  }
}

f_38_modularity<-function(ceb){
  return(modularity(ceb))
}

f_39_meanFST<-function(mat){
  if (length(mat)>1) {
    return(mean(mat,na.rm=T))
  }else{
    return(0)
  }
}

f_40_maxFST<-function(mat){
  if (length(mat)>1) {
    return(max(mat,na.rm=T))
  }else{
    return(0)
  }
}

f_41_minFST<-function(mat){
  if (length(mat)>1) {
    return(min(mat,na.rm=T))
  }else{
    return(0)
  }
}


calculateFeatures <- function(OriginalG, useOriginal = T, cutoff = 0.99, similarity=F){
  if (useOriginal){
    g<-OriginalG
  }else{
    if (similarity){ #if the cutoff means similarity cutoff
      g<-delete.edges(OriginalG,which(E(OriginalG)$weight<cutoff))
    }else{ #cutoff means quantile cutoff
      g<-delete.edges(OriginalG,which(E(OriginalG)$weight<quantile(E(OriginalG)$weight,cutoff)))
    }
    
  }
  
  featureVector <- vector(length=41)
  # Diagnostics of transitivity
  featureVector[1] <- f_01_averageLocalClusteringCoeff(g,F)    # Clustering coefficient averaged across all nodes
  featureVector[2] <- f_02_averageLocalClusteringCoeffWeighted(g,F)    # Barrat's clustering coefficient averaged across all nodes
  featureVector[3] <- f_03_globalClusteringCoeff(g,F)
  # Diagnostics of degree/sterngth
  featureVector[4] <- f_04_gdensity(g,T)                    # Graph density
  featureVector[5] <- f_05_proportionSingletons(g)             # Proportion of nodes with degree 0 of all the nodes
  featureVector[6] <- f_06_proportionEndpoints(g)              # Proportion of nodes with degree 1 of all the nodes
  featureVector[7] <- f_07_meanDegree(g)                       # Average degree
  featureVector[8] <- f_08_assortativityDegree(g)
  featureVector[9] <- f_09_meanStrength(g)                     # Average strength
  featureVector[10] <- f_10_straightness(g)                    # power-law test
  featureVector[11] <- f_11_entropyDegreeDistribution(g)       # Average measurement of the heterogeneity of the network. See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  featureVector[12] <- f_12_ratioComponents(g)                 # Number of components relative to networks size
  featureVector[13] <- f_13_giantConnectedRatio(g)
  featureVector[14] <- f_14_evennessComponentSize (g)
  #Diagnostics of shortest-paths
  featureVector[15] <- f_15_CentralPointDominance(g)             
  featureVector[16] <- f_16_meanEccentricity(g)                 # Eccentricity is the maximum shortest distance from a node to all other nodes. This is averaged across all nodes
  featureVector[17] <- f_17_gdiameter(g)                        # length of the longest geodesic
  featureVector[18] <- f_18_meanDiameterComponents(g)
  featureVector[19] <- f_19_globalEfficiency(g)                # Latora and Marchiori; See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  featureVector[20] <- f_20_averageClosenessCentrality(g)      #
  # Diagnostics of motifs
  featureVector[21:32] <- motifsProportion(g)
  # Diagnostics of reciprocity
  featureVector[33]<-f_33_reciprocity(g)
  featureVector[34]<-f_34_inoutCorrelation(g)                  #correlation between in and out edges per node
  # Diagnostics of communities
  commun<-calCommunity(OriginalG, useOriginal, cutoff,F)
  ceb<-commun[[1]]
  mat<-commun[[2]]
  featureVector[35]<-f_35_communitySizeEvenness(ceb)
  featureVector[36]<-f_36_numberCommonCommunities(ceb)
  featureVector[37]<-f_37_CommunityRatio(ceb)
  featureVector[38]<-f_38_modularity(ceb)
  featureVector[39]<-f_39_meanFST(mat)
  featureVector[40]<-f_40_maxFST(mat)
  featureVector[41]<-f_41_minFST(mat)
  return(featureVector)
}

calculateDegreeDist <- function(OriginalG, useOriginal = T, cutoff = 0.99, similarity=F){
  if (useOriginal){
    g<-OriginalG
  }else{
    if (similarity){ #if the cutoff means similarity cutoff
    g<-delete.edges(OriginalG,which(E(OriginalG)$weight<cutoff))
  }else{ #cutoff means quantile cutoff
    g<-delete.edges(OriginalG,which(E(OriginalG)$weight<quantile(E(OriginalG)$weight,cutoff)))
  }
  }
  
  dd_in<-degree_distribution(g, mode= "in")
  dd_out<-degree_distribution(g, mode= "out")
  return(list(dd_in, dd_out))
}



