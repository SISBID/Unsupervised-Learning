#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Unsupervised Lab
############################################################


####################
#Data Description
#############
#gdat - Gene Expression Data - n = 445 patients x p = 353 genes
#Only 353 genes with somatic mutations from COSMIC are retained
#Data is Level III TCGA BRCA RNA-Sequencing gene expression data that have already been pre-processed according to the following steps:
#1) reads normalized by RPKM
#2) samples are quantile normalized
#3) corrected for overdispersion by a power-transformation (using a Kolmogorov-Smirnoff goodness of fit test to a Poisson distribution)
#Short gene name labels are given as the column names
###############
#cdat - Clinical Data - n = 445 patients x q = 6 clinical features
#1) Subtype - denotes 5 PAM50 subtypes including Basal-like, Luminal A, Luminal B, HER2-enriched, and Normal-like
#2) ER-Status - estrogen-receptor status
#3) PR-Status - progesterone-receptor status
#4) HER2-Status - human epidermal growth factor receptor 2 status
#5) Node - number of lymph nodes involved
#6) Metastasis - indicator for whether the cancer has metastasized
###############


############
#Problem 1 - Dimension reduction 
############
#Problem 1a - Apply PCA, NMF, ICA and MDS to this dataset. Compare and contrast the results using these methods.  

#Problem 1b - Relate the dimension reduction results with the clinical data. Is any clinical information reflected in the lower dimensional spaces?

#Problem 1c - Overall, which dimension reduction method do you recommend for this data set and why?

############
#Problem 2 - Clustering 
############
#Problem 2a - Apply various clustering algorithms such as K-means (explore different K), hierarchical clustering (explore different linkages), NMF, and biclustering. Compare the clustering results using these methods.

#Problem 2b - Relate the clustering results with the clinical data. Can the clustering algorithm recover some of the clinical information such as cancer subtypes?

#Problem 2c - Overall, which clustering method(s) do you recommend for this data set and why?

############
#Problem 3 - Multiple comparisons 
############
#Problem 3a - Identify important genes to differetiate ER postive and negative, PR postive and negative, HER2 postive and negative, metastasis status.

#Problem 3b - Try different procedures to adjust for multiple comparisons.

#Problem 3c - Examine the lists of genes identified using different methods for each clinical response. 

############
#Problem 4 - Graphical models
############

# Use graphical models to explore the interaction among genes.




###############################################################
###############################################################
#R scripts to help out with the BRCA case study Lab
#Don't peek at this if you want to practice coding on your own!!
##################################################################

#load data
load("brcadat.Rdata")

#explore data
dim(gdat)
dim(cdat)

table(cdat$Subtype)
table(cdat$ER)
table(cdat$PR)
table(cdat$HER2)
table(cdat$Node)
table(cdat$Metastasis)
table(cdat$ER,cdat$PR)

hist(c(gdat),breaks=100)

################
#visualize data

#cluster heatmap - biclustering
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:2],bb[1:25],cc[1:50])]

# without scaling
heatmap(gdat,col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))
# with scaling
heatmap(scale(gdat),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))

Cols=function(vec){cols=rainbow(length(unique(vec)))
                   return(cols[as.numeric(as.factor(vec))])}

heatmap(scale(gdat),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"),labRow=cdat$Subtype,RowSideColors=Cols(cdat$Subtype))


####################
#pattern recognition
#dimension reduction

Cols=function(vec){cols=rainbow(length(unique(vec)))
                   return(cols[as.numeric(as.factor(vec))])}


#PCA
sv = svd(gdat)
V = sv$v
Z = gdat%*%V

K = 3
pclabs = c("PC1","PC2","PC3","PC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(Z[,i],Z[,j],pch=16,xlab=pclabs[i],ylab=pclabs[j],col=Cols(cdat$Subtype))
}
legend(-30,-20,pch=16,col=rainbow(5),levels(cdat$Subtype))

#MDS 
dev.off()
Dmat = dist(gdat)
mdsres = cmdscale(Dmat,k=2)
plot(mdsres[,1],mdsres[,2],pch=16,col=Cols(cdat$Subtype))
legend(-20,30,pch=16,col=rainbow(5),levels(cdat$Subtype))


#NMF
require("NMF")

K = 4
nmffit = nmf(gdat,rank=K) 
W = basis(nmffit)
H = coef(nmffit)

kk = 3
pclabs = c("NMF1","NMF2","NMF3","NMF4")
par(mfrow=c(1,kk))
for(i in 1:kk){
  j = i+1
  plot(W[,i],W[,j],pch=16,xlab=pclabs[i],ylab=pclabs[j],col=Cols(cdat$Subtype))
}
legend(0,65,pch=16,col=rainbow(5),levels(cdat$Subtype))

par(mfrow=c(1,2))
basismap(nmffit,annRow=cdat$Subtype,scale="col",legend=FALSE)
coefmap(nmffit,annCol=colnames(gdat),scale="col",legend=FALSE)


#ICA
require("fastICA")

K = 4
icafit = fastICA(gdat,n.comp=K)

kk = 3
pclabs = c("ICA1","ICA2","ICA3","ICA4")
par(mfrow=c(1,kk))
for(i in 1:kk){
  j = i+1
  plot(icafit$A[i,],icafit$A[j,],pch=16,xlab=pclabs[i],ylab=pclabs[j],col=Cols(cdat$Subtype))
}
legend(-4,-2,pch=16,col=rainbow(5),levels(cdat$Subtype))

#######################
#clustering

#K-means
K = 5
km = kmeans(gdat,centers=K)
table(km$cluster,cdat$Subtype)


#hierarchical
#which linakges is the best?
#which distance metric is the best?

Dmat = dist(gdat)
com.hc = hclust(Dmat,method="ward.D")

dev.off()
plot(com.hc,labels=cdat$Subtype,cex=.5)

res.com = cutree(com.hc,5)
table(res.com,cdat$Subtype)


#################
#genes significantly associated with ER or PR Status, etc

x = gdat[cdat$ER=="Positive" | cdat$ER=="Negative",]
y.er = cdat$ER[cdat$ER=="Positive" | cdat$ER=="Negative"]
y.label<-rep(1, length(y.er))
y.label[y.er == "Positive"]=2


ps <- NULL
for(i in 1:ncol(gdat)) ps <- c(ps,
 t.test(x[y.label==1,i],x[y.label==2,i])$p.value)
fdrs.bh <- p.adjust(ps, method="BH")

cat("Number of Tests significant with alpha=0.1 using Bonferroni correction:",
sum(ps<0.1/length(y.label)), fill=TRUE)

cat("Number of Tests with FDR below 0.1:",
sum(fdrs.bh<0.1), fill=TRUE)

plot(sort(ps,decreasing=FALSE),ylab="P-Values")
#BH procedure
abline(a=0, b=0.1/length(y.label),col="red")
#Bonferroni
abline(a=0.1/length(y.label), b=0,col="blue")


#SAM
library(samr)

# samfit x: p by n
samfit<-SAM(t(x),y.label,resp.type="Two class unpaired", fdr.output=0.1)    

# examine significant gene list 
 
print(samfit)
 # plot results
plot(samfit)      


#######################
#graphical models - how are genes related?

#this takes a bit of time
require("igraph")
require("XMRF")

lam = lambdaMax(gdat)*sqrt(log(ncol(gdat))/nrow(gdat))*0.01
net = XMRF(t(gdat),method="LPGM",lambda.path=lam,N=1,th=.0025)
plot(net)






