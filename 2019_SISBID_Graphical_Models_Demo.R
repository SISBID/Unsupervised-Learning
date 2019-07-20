#############################################################
#2019 SISBID Module 3 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu & Ali Shojaie 
#Clustering Lab
############################################################

############
#Data set - Author Data. This data set consists of word counts from chapters written by four British authors.
#############

require("igraph")
require("XMRF")
require("huge")


load("UnsupL_SISBID_2019.Rdata")

dim(author)
colnames(author)
unique(rownames(author))
TrueAuth = as.factor(rownames(author))

X = log(1 + author[,1:69]) #log transform data


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


#Poisson Graphical Models via XMRF package

lam = lambdaMax(X)*sqrt(log(ncol(X))/nrow(X))*0.02
net = XMRF(t(X),method="LPGM",lambda.path=lam,N=1,th=.001)
mat = net$network[[1]]
neti = as.undirected(graph_from_adjacency_matrix(mat))
plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_with_kk)

plot(neti,vertex.label=colnames(X),vertex.size=2,vertex.label.cex=1.2,vertex.label.dist=1,layout=layout_in_circle)


