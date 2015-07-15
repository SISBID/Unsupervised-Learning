#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Clustering Demos for use in lecture
############################################################



########################################################
#Data set 1 - Simulated Data
#Small simulated data set to demonstrate concepts with k-means clustering
##########################################################

#simulate data
n = 300
mu1 = c(3,3); mu2 = c(7,4); mu3 = c(6,5.5); 
Sig = matrix(c(1,.5,.5,1),2,2)
x1 = t(matrix(mu1,2,n/3)) + matrix(rnorm(n),n/3,2)
xx = matrix(rnorm(n*2/3),n/3,2)
x2 = t(matrix(mu2,2,n/3)) + xx%*%chol(Sig)
xx = matrix(rnorm(n*2/3),n/3,2)
x3 = t(matrix(mu3,2,n/3)) + xx%*%chol(Sig)
x = rbind(x1,x2,x3)
Y = c(rep(1,n/3),rep(2,n/3),rep(3,n/3))
y = factor(Y)

plot(x[,1],x[,2],col=as.numeric(y),pch=16)


#try changing k - which clustering looks best?
k = 3
km = kmeans(x,centers=k)
plot(x[,1],x[,2],col=km$cluster,pch=16)
cens = km$centers
points(cens[,1],cens[,2],col=1:k,pch=16,cex=3)


#code to understand K-means algorithm
require("animation")
mv.kmeans = function(x,k,cens=NULL){
  n = nrow(x)
  if(is.null(cens)){
      cens = x[sample(1:n,k),]
    }
  plot(x[,1],x[,2],pch=16)
  points(cens[,1],cens[,2],col=1:k,pch=16,cex=3)
  thr = 1e-6; ind = 1; iter = 1;
  while( ind>thr)
    {
      oldcen = cens
      km = kmeans(x,centers=cens,iter.max=1,nstart=1,algorithm="MacQueen")
      plot(x[,1],x[,2],col=km$cluster,pch=16)
      points(cens[,1],cens[,2],col=1:k,pch=16,cex=3)
      cens = km$centers
      #print(cens)
      plot(x[,1],x[,2],col=km$cluster,pch=16)
      points(cens[,1],cens[,2],col=1:k,pch=16,cex=3)
      ind = sum(diag((oldcen-cens)%*%t(oldcen-cens)))
      #print(ind)
    } 
}

#watch K-means algorithm movie
#start from random starting points
saveHTML(mv.kmeans(x,3,cens=NULL),img.name="km1")


########################################################
#Data set 2 - NCI Microarray data
#Apply K-means to cluster a high-dimensional data set.
#Apply hierarchical clustering & try out different linkages.
#Apply biclustering (Cluster heatmap) to visualize data.
##########################################################

require("ISLR")

ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

dim(ncidat)
unique(colnames(ncidat))

############
#apply K-means
K = 9
km = kmeans(t(ncidat),centers=K)

#how do we visualize K-means results?

#PCA - take SVD to get solution
#center genes, but don't scale
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = svd(t(X));
U = sv$u
V = sv$v
D = sv$d
Z = t(X)%*%V;

plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

#Re-run and see if solution changes
K = 9
km = kmeans(t(ncidat),centers=K)
plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

#try different K
K = 5
km = kmeans(t(ncidat),centers=K)
plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)


##################
#hierarchical clustering

require("ISLR")

ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

dim(ncidat)
unique(colnames(ncidat))

#complete linakge - Euclidean distance
cols = as.numeric(as.factor(colnames(ncidat)))
Dmat = dist(t(ncidat))
com.hclust = hclust(Dmat,method="complete")
plot(com.hclust,cex=.7,main="Complete Linkage")

#single linakge
dev.new()
sing.hclust = hclust(Dmat,method="single")
plot(sing.hclust,cex=.7,main="Single Linkage")

#average linakge
dev.new()
ave.hclust = hclust(Dmat,method="average")
plot(ave.hclust,cex=.7,main="Average Linkage")

#Ward's linakge
dev.new()
ward.hclust = hclust(Dmat,method="ward.D")
plot(ward.hclust,cex=.7,main="Ward's Linkage")

#complete linkage with different distances
dev.new()
Dmat = dist(t(ncidat),method="manhattan") #L1 distance
com.hclust = hclust(Dmat,method="complete")
plot(com.hclust,cex=.7,main="Complete Linkage - L1 Dist")


##########
#Biclustering - Cluster Heatmap 

require("ISLR")
ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

#filter genes using PCA
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = svd(t(X));
V = sv$v

#PC loadings - visualize data by limiting to top genes in magnitude in the PC loadings 
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 2
ord = order(abs(V[,j]),decreasing=TRUE)
x = as.matrix(X[ord[1:250],])

#cluster heatmap - uses Ward's linkage (complete is default)
heatmap(x,col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))

######################################################
