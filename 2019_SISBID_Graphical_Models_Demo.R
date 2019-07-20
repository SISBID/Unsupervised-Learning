#############################################################
#2019 SISBID Module 3 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu & Ali Shojaie 
#Clustering Lab
############################################################

require("igraph")
require("XMRF")
require("huge")
require("glasso")
require("spacejam")
require("WGCNA")
require("glmnet")

rm(list=ls())

##
## Read the Sachs et al data
##
#setwd('/Users/ashojaie/Dropbox/Teaching/shortcourse_NET/Ali/')

sachscov <- as.matrix(read.table("sachscov.txt"))
dim(sachscov)

sachscor <- cov2cor(sachscov)

sachsdat <- as.matrix(read.table("sachs.data"))
dim(sachsdat)
head(sachsdat)
ps <- c("praf","pmek","plcg","PIP2","PIP3","P44","pakts","PKA","PKC","P38","pjnk")
colnames(sachsdat) <- ps
#X2 <- as.data.frame(scale(sachsdat))

p <- ncol(sachsdat)
n <- nrow(sachsdat)

load("UnsupL_SISBID_2019.Rdata")

############
#Data set - Author Data. This data set consists of word counts from chapters written by four British authors.
#############

dim(author)
colnames(author)
unique(rownames(author))
TrueAuth = as.factor(rownames(author))

X = log(1 + author[,1:69]) #log transform data


##
## coexpression network
##
tau <- 0.25
A1 <- 1 * (round(abs(sachscor),3) > tau); diag(A1) <- 0
sum(A1)/2

##
## WGCNA
##
A2 <- round(WGCNA::adjacency(sachsdat), 3); diag(A2) <- 0
A2[1:5,1:5]
sum((A2 > 0))/2

## co-expression network

## method 1: simple thresholding of correlations, at a cutoff chosen to give similar number of 
##			 edges to partial correlation methods

## a randomly chosen threshold, to get 50 nonzero values in the adjacency matrix (25 edges)
tau <- 0.22
A2 <- abs(sachscor) > tau; diag(A2) <- 0
sum(A2)/2

## method 2: testing for nonzero correlations

## testing for nonzero correlation, using Fisher Z-transform
fisherzs <- atanh(sachscor)
fisherps  <- 2*pnorm(abs(fisherzs), 0, 1/sqrt(n-3), lower.tail=FALSE)
A3 <- fisherps < (0.01/(p*(p-1)/2)); diag(A3) <- 0
sum(A3)/2

## 
## plot the three networks
##
g1 <- graph.adjacency(A1, mode="undirected")
g2 <- graph.adjacency(A2, mode="undirected")
g3 <- graph.adjacency(A3, mode="undirected")
g4 <- graph.adjacency(A4, mode="undirected")

pdf('plot.pdf', width=9, height=3)
par(mfrow = c(1,3))
plot(g1,layout=layout.circle(g1), main='A1')
plot(g2,layout=layout.circle(g2), main='A2')
plot(g3,layout=layout.circle(g3), main='A3')
plot(g4,layout=layout.circle(g4), main='A4')
dev.off()
g0 <- g2

##
## Partial correlation networks
##
abs(round(sachscor,3))
abs(round(solve(sachscor),3))

invcov <- abs(round(solve(sachscor),3)); 
invcov <- cov2cor(invcov)
invcov <- 1*(invcov > 0.39); diag(invcov) <- 0; sum(invcov)/2
g0 <- graph.adjacency(invcov, mode="undirected")
plot(g0)

## calculate lambda, based on formula in the slides
alpha = 0.01
num = qt(p=alpha/(2*(p^2)),df=n-2, lower.tail=F)
lambda = num / sqrt(n-2 + num)

##
## glasso
##
glasso.est <- glasso(s=sachscor, rho=lambda*4.2, approx=FALSE, penalize.diagonal=FALSE)
A1 <- abs(glasso.est$wi) > 1E-16; diag(A1) <- 0
g1 <- graph.adjacency(A1, mode="undirected")

##
## neighborhood selection
##
ns.est <- glasso(s=sachscor, rho=lambda, approx=TRUE, penalize.diagonal=FALSE)
A2 <- abs(ns.est$wi) > 1E-16; diag(A2) <- 0
g2 <- graph.adjacency(A2, mode="undirected")


#Gaussian graphical models (via neighborhood selection)
#huge package

neth = huge(X,method="mb")
plot(neth)

dev.off()

#larger lambda
mat = neth$path[[3]]
neti = as.undirected(graph_from_adjacency_matrix(mat))
plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_with_kk)

#smaller lambda
mat = neth$path[[4]]
neti = as.undirected(graph_from_adjacency_matrix(mat))
plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_with_kk)

##
## nonparanormal
##
scor <- cor(sachsdat,method='spearman')
scor <- 2*sin(scor*pi/6)
npn.est <- glasso(s=scor, rho=lambda, approx=FALSE, penalize.diagonal=FALSE)
A3 <- abs(npn.est$wi) > 1E-16; diag(A3) <- 0
g3 <- graph.adjacency(A3, mode="undirected")

##
## nonparanormal -- alternative estiamtion 
##
npn.cor <- huge.npn(x=sachsdat, npn.func="skeptic", npn.thresh=NULL, verbose=FALSE)
npn.est <- glasso(s=npn.cor, rho=lambda, penalize.diagonal=FALSE)
A4 <- abs(npn.est$wi) > 1E-16; diag(A4) <- 0
g4 <- graph.adjacency(A4, mode="undirected")

##
## binary network estimation
##
#head(sachsdat)

sachsbin <- 1*(sachsdat > 0) + -1*(sachsdat <= 0)
head(sachsbin)

bin.est <- matrix(0,p,p)
## estiamte the neighborhood for each node 
for(j in 1:p){
  ## this is the same method used in neighborhood selection, the only difference is 'family'
  nbr <- glmnet(x=sachsbin[,-j], y=sachsbin[,j], family='binomial', lambda=lambda) 
  bin.est[j,-j] <- 1*(abs(as(nbr$beta,"matrix")) > 0)	#store the estimates in jth row of matrix
}
A6 <- bin.est; diag(A6) <- 0
sum(A6)/2
g6 <- graph.adjacency(A6, mode="undirected")

## 
## plot the networks
##
pdf('plot.pdf', width=9, height=6)
par(mfrow = c(2,3), mar=c(1,1,4,1))
#plot(g0,layout=layout.circle(g0), main='co-expression')
#plot(g0,layout=layout.circle(g0), main='direct inv')
plot(g1,layout=layout.circle(g1), main='glasso')
plot(g2,layout=layout.circle(g2), main='NS')
plot(g3,layout=layout.circle(g3), main='nonparanormal')
plot(g4,layout=layout.circle(g4), main='nonparanormal - v2')
plot(g5,layout=layout.circle(g5), main='spacejam')
plot(g6,layout=layout.circle(g3), main='Binary')
dev.off()

##
## Poisson Graphical Models via XMRF package
##
lam = lambdaMax(X)*sqrt(log(ncol(X))/nrow(X))*0.02
net = XMRF(t(X),method="LPGM",lambda.path=lam,N=1,th=.001)
mat = net$network[[1]]
neti = as.undirected(graph_from_adjacency_matrix(mat))
plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_with_kk)

plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_in_circle)

##
## spacejam
##
spacejam.est <- SJ(sachsdat, lambda=0.5)
A5 <- 1*(spacejam.est$G)[,,1]
g5 <- graph.adjacency(A5, mode="undirected")
