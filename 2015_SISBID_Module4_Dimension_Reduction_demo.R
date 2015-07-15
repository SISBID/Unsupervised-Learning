#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Dimension Reduction Demos for use in lecture
############################################################



########################################################
#Data set 1 - College Data
#Small data set to understand R's built in PCA functions
#princomp & prcomp
##########################################################

#read in data
library(ISLR)
data(College)
cdat = College[,2:18]

dim(cdat)
names(cdat)

#PCA
pc = princomp(cdat) #default - centers and scales

#default R plots with princomp
biplot(pc,cex=.7)
screeplot(pc)

#scatter plots - patterns among observations
i = 1; j = 2;
plot(pc$scores[,i],pc$scores[,j],pch=16,cex=.2)
text(pc$scores[,i],pc$scores[,j],rownames(cdat),cex=.6)

#look at a particular college
ind = match("Harvard University",rownames(cdat))
text(pc$scores[ind,i],pc$scores[ind,j],rownames(cdat)[ind],cex=.7,col=2)

#loadings - variables that contribute to these patterns
par(mfrow=c(2,1))
barplot(pc$loadings[,1],cex.names=.6,main="PC 1 Loadings")
barplot(pc$loadings[,2],cex.names=.6,main="PC 2 Loadings")

#variance explained
screeplot(pc)

varex = 100*pc$sdev^2/sum(pc$sdev^2)
plot(varex,type="l",ylab="% Variance Explained",xlab="Component")

#cumulative variance explained
cvarex = NULL
for(i in 1:ncol(cdat)){
  cvarex[i] = sum(varex[1:i])
}
plot(cvarex,type="l",ylab="Cumulative Variance Explained",xlab="Component")

######
#sparse PCA
library(PMA)

spc = SPC(scale(cdat),sumabsv=2,K=3)

spcL = spc$v
rownames(spcL) = names(cdat)

#scatterplots of Sparse PCs
i = 1; j = 2;
plot(spc$u[,i],spc$u[,j],pch=16,cex=.2)
text(spc$u[,i],spc$u[,j],rownames(cdat),cex=.6)

#loadings 
par(mfrow=c(2,1))
barplot(spc$v[,1],names=names(cdat),cex.names=.6,main="SPC 1 Loadings")
barplot(spc$v[,2],names=names(cdat),cex.names=.6,main="SPC 2 Loadings")

#variance explained
spc$prop.var.explained


##########################################################
#Dataset 2 - NCI Microarray Data
#Understand PCA and Sparse PCA
#PCA solution via the SVD
###########################################################

require("ISLR")

ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

dim(ncidat)
unique(colnames(ncidat))

#PCA - take SVD to get solution
#center genes, but don't scale
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = svd(t(X));
U = sv$u
V = sv$v
D = sv$d
Z = t(X)%*%V;

#PC scatterplots
cols = as.numeric(as.factor(colnames(ncidat)))
K = 3
pclabs = c("PC1","PC2","PC3","PC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(U[,i],U[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(U[,i],U[,j],colnames(X),col=cols)
}

#PC loadings - visualize data by limiting to top genes in magnitude in the PC loadings 
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 2
ord = order(abs(V[,j]),decreasing=TRUE)
x = as.matrix(X[ord[1:250],])
heatmap(x,col=gcol2)

#Variance Explained
varex = 0
cumvar = 0
denom = sum(D^2)
for(i in 1:64){
  varex[i] = D[i]^2/denom
  cumvar[i] = sum(D[1:i]^2)/denom
}

#screeplot
par(mfrow=c(1,2))
plot(1:64,varex,type="l",lwd=2,xlab="PC",ylab="% Variance Explained")
plot(1:64,cumvar,type="l",lwd=2,xlab="PC",ylab="Cummulative Variance Explained")

#######
#Sparse PCA
require("PMA")

spc = SPC(t(X),sumabsv=10,K=4)

#how many genes selected?
apply(spc$v!=0,2,sum)

#PC scatterplots
cols = as.numeric(as.factor(colnames(ncidat)))
K = 3
pclabs = c("SPC1","SPC2","SPC3","SPC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(spc$u[,i],spc$u[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(spc$u[,i],spc$u[,j],colnames(X),col=cols)
}

#SPC loadings - visualize data by limiting to gene selected by the sparse PC loadings
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 1
ind = which(spc$v[,j]!=0)
x = as.matrix(X[ind,])
heatmap(x,col=gcol2)

#variance explained
spc$prop.var.explained


##########################################################
#Dataset 3 - Digits Data
#Here only use 3's to compare and contrast PCA, NMF and ICA
###########################################################

load("UnsupL.Rdata")

#pull out 3's
dat3 = digits[which(rownames(digits)==3),]

#visulaize
par(mfrow=c(3,4))
for(i in 1:12){
  imagedigit(dat3[i,])
}

#PCA - take SVD to get solution
#don't center and scale to retain interpretation as images
svd3 = svd(dat3)
U = svd3$u
V = svd3$v #PC loadings
D = svd3$d
Z = dat3%*%V #PCs

#PC scatterplot
par(mfrow=c(1,1))
plot(Z[,2],Z[,3],pch=16)

#PC loadings
par(mfrow=c(1,4))
for(i in 1:4){
  imagedigit(V[,i])
}

#Variance Explained
varex = 0
cumvar = 0
denom = sum(D^2)
for(i in 1:256){
  varex[i] = D[i]^2/denom
  cumvar[i] = sum(D[1:i]^2)/denom
}

#screeplot
par(mfrow=c(1,2))
plot(1:256,varex,type="l",lwd=2,xlab="PC",ylab="% Variance Explained")
plot(1:256,cumvar,type="l",lwd=2,xlab="PC",ylab="Cummulative Variance Explained")


cumvar[25] #first 25 PCs explain over 90% of variance
pdat3 = dat3%*%V[,1:25] #projected data - a tenth of the original size


#######
#now NMF

require("NMF")

K = 10
nmffit = nmf(dat3+1,rank=K)
W = basis(nmffit)
H = coef(nmffit)

#plot archetypes - try changing K
par(mfrow=c(3,5))
for(i in 1:K){
  imagedigit(H[i,])
}


###########
#now ICA

require("fastICA")

K = 10
icafit = fastICA(t(dat3),n.comp=K)

#plot independent source signals - try changing K
par(mfrow=c(3,5))
for(i in 1:K){
  imagedigit(icafit$S[,i])
}

#################################################
