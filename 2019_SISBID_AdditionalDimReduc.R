######################################################
## Dimension Reduction
######################################################
library(MASS)
library(ISLR)
library(cluster)

rm(list=ls())

##
## PCA
##

## PCA for USArrests
data(USArrests)
pairs(USArrests, col=4, cex=0.7, upper.panel=panel.smooth)

states <- row.names(USArrests)
states[1:5]

apply(USArrests, 2, mean)
apply(USArrests, 2, var)

pc.out <- prcomp(USArrests, scale=TRUE)
print(pc.out$rot)

pc.out <- prcomp(USArrests, scale=TRUE, retx=TRUE)
names(pc.out)

## plotting the observations and variables in the PC space
plot(pc.out$x[,1], pc.out$x[,2], type="n", xlab="1st PC", ylab="2nd PC")
text(pc.out$x[,1], pc.out$x[,2], labels=states, cex=0.7, col=4)

plot(pc.out$rot, type="n", xlim=c(-0.6, -0.2))
text(pc.out$rot, names(USArrests), col=4)

## biplot
biplot(pc.out)

## scree plots & PVE
screeplot(pc.out)

par(mar=c(2.5,4.5,.5,.5))
plot(pc.out, col='orange', main='', type='l')
axis(1, at=c(0.7, 1.8, 3.2, 4.5), labels=paste('PC', 1:4))

print(pc.out$sdev^2)
print((pc.out$sdev^2)/sum(pc.out$sdev^2))

par(mfrow=c(1,2), mar=c(4.5,4.5,.5,.5))
plot((pc.out$sdev^2)/sum(pc.out$sdev^2), xlab="PC", ylab="PVE", 
	col=4, ylim=c(0, 1), type='b')
plot(cumsum(pc.out$sdev^2)/sum(pc.out$sdev^2), xlab="PC", col=4, 
	ylab="Cumulative PVE", ylim=c(0, 1), type='b')


## PCA for College data
data(College)
cdat = College[,2:18]
dim(cdat)
names(cdat)

## PCA
pc.col <- princomp(cdat) #default - centers and scales

#default R plots with princomp
biplot(pc.col, cex=.7)
screeplot(pc.col)

# #scatter plots - patterns among observations
i = 1; j = 2;
plot(pc.col$scores[,i],pc.col$scores[,j],pch=16,cex=.2)
text(pc.col$scores[,i],pc.col$scores[,j],rownames(cdat),cex=.6)

#look at a particular college
ind = match("Harvard University",rownames(cdat))
text(pc.col$scores[ind,i],pc.col$scores[ind,j],rownames(cdat)[ind],cex=.7,col=2)

#loadings - variables that contribute to these patterns
par(mfrow=c(2,1))
barplot(pc.col$loadings[,1],cex.names=.6,main="PC 1 Loadings")
barplot(pc.col$loadings[,2],cex.names=.6,main="PC 2 Loadings")

#variance explained
screeplot(pc.col)

varex <- 100*pc.col$sdev^2/sum(pc.col$sdev^2)
plot(varex,type="l",ylab="% Variance Explained",xlab="Component")

#cumulative variance explained
cvarex = NULL
for(i in 1:ncol(cdat)){
  cvarex[i] = sum(varex[1:i])
}
plot(cvarex,type="l",ylab="Cumulative Variance Explained",xlab="Component")

######################################################


##
## sparse PCA
##
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

######################################################


##
##MDS
##
## US arrests data
dist.1 <- daisy(USArrests)
dist.2 <- daisy(USArrests, metric="manhattan")

mds.1 <- cmdscale(dist.1)
mds.2 <- cmdscale(dist.2)
mds.3 <- isoMDS(dist.1)

#plotting the MDS resutls
par(mfrow=c(1,3), mar=c(4.5,4.5,2.5,.5))
x <- mds.1[,1]
y <- mds.1[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Euclidean Dist')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)

x <- mds.2[,1]
y <- mds.2[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Manhattan Dist')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)

x <- mds.3$points[,1]
y <- mds.3$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "n", main='Non-Metric MDS')
text(x, y, labels = rownames(USArrests), col='blue', cex=0.7)

######################################################


##
## Case Study: Analysis of NCI60 data
##

## PCA
library(ISLR)
data(NCI60)
dat <- NCI60$data
#dim(dat)

labs <- NCI60$labs
#table(labs)

##remove subtypes with only one member
indx <- names(table(labs))[table(labs) == 1]
dat <- dat[!(labs %in% indx),]
labs2 <- labs[!(labs %in% indx)]

dat2 <- t(dat)

##select subset of data and order the samples
ord <- order(labs2)
dat <- dat[ord,]
labs2 <- labs2[ord]

# incorrect PCA (if the goal is finding patterns in samples)
pc.nci.1 <- prcomp(dat, scale=TRUE, retx=TRUE)
# summary(pc.nci)

# correct PCA (if the goal is finding patterns in samples)
pc.nci.2 <- prcomp(dat2, scale=TRUE, retx=TRUE)
# summary(pc.nci)

#scree plot
par(mfrow=c(1,2), mar=c(4.5, 4.5, 1, .5))
plot((pc.nci.1$sdev^2)/sum(pc.nci.1$sdev^2), xlab="PC", ylab="PVE", 
	col=4, ylim=c(0, 1), type='b', cex=0.5)
#plot(cumsum(pc.nci$sdev^2)/sum(pc.nci$sdev^2), xlab="PC", col=4, 
#	ylab="Cumulative PVE", ylim=c(0, 1), type='b')
screeplot(pc.nci.1, type='l', col='blue', main='') 

#plotting the observations in the PC space
mycol <- as.integer(as.factor(labs2))
par(mfrow=c(1,3), mar=c(4.5,4.5,.5,.5))
plot(pc.nci.1$x[,1], pc.nci.1$x[,2], type="p", xlab="1st PC", ylab="2nd PC",
    cex=1.2, col=mycol, pch=mycol)
#legend('topright', unique(labs2), col=unique(mycol), pch=unique(mycol), cex=0.6)
plot(pc.nci.1$x[,1], pc.nci.1$x[,3], type="p", xlab="1st PC", ylab="3rd PC",
    cex=1.2, col=mycol, pch=mycol)
plot(pc.nci.1$x[,2], pc.nci.1$x[,3], type="p", xlab="2nd PC", ylab="3rd PC",
    cex=1.2, col=mycol, pch=mycol)


## PCA - direct calculation with SVD

# center genes, but don't scale
X = t(scale(dat,center=TRUE,scale=FALSE))
sv = svd(t(X));
U = sv$u
V = sv$v
D = sv$d
Z = t(X)%*%V;

## PC scatterplots
cols = as.integer(as.factor(labs2))
K = 3
pclabs = c("PC1","PC2","PC3","PC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(U[,i],U[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(U[,i],U[,j],colnames(X),col=cols)
}

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

spc = SPC(t(X),sumabsv=7,K=4)

#how many genes selected?
apply(spc$v!=0,2,sum)

#PC scatterplots
K = 3
pclabs = c("SPC1","SPC2","SPC3","SPC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(spc$u[,i],spc$u[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(spc$u[,i],spc$u[,j],colnames(X),col=cols)
}

## SPC loadings - visualize data by limiting to gene selected by the sparse PC loadings
aa = grep("grey",colors())
bb = grep("blue",colors())
cc = grep("yellow",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 1
ind = which(spc$v[,j]!=0)
x = as.matrix(X[ind,])
heatmap(x,col=gcol2)

#variance explained
spc$prop.var.explained


## MDS

dist.1 <- daisy(dat)
dist.2 <- daisy(dat, metric="manhattan")

mds.1 <- cmdscale(dist.1)
mds.2 <- cmdscale(dist.2)
mds.3 <- isoMDS(dist.1)

par(mfrow=c(1,3), mar=c(4.5,4.5,2.5,.5))
x <- mds.1[,1]
y <- mds.1[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "p", main='Euclidean Dist', pch=mycol,col=mycol)
legend('topright', unique(labs2), col=unique(mycol), pch=unique(mycol), cex=0.6)

x <- mds.2[,1]
y <- mds.2[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "p", main='Manhattan Dist', pch=mycol,col=mycol)

x <- mds.3$points[,1]
y <- mds.3$points[,2]
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
	xlim = range(x)*1.2, type = "p", main='Non-Metric MDS', pch=mycol,col=mycol)
######################################################
