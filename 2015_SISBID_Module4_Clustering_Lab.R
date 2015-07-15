#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Clustering Lab
############################################################

############
#Data set - Author Data. This data set consists of word counts from chapters written by four British authors.
#This lab will put together concepts from both dimension reduction and clustering.
#There are ultimately 3 goals to this lab:
#1) Correctly cluster author texts in an unsupervised manner.
#2) Determine which words are responsible for correctly separating the author texts.
#3) Visualize the author texts, words and the results of your analysis
#############


#############
#Problem 1 - Visualization
#############
#Problem 1a - We wish to plot the author texts as well as the words via a 2D scatterplot.  Which method would be best to use?  Why?

#Problem 1b - Apply PCA to visualize the author texts.  Explain the results.

#Problem 1c - Apply MDS to visualize the author texts.  Interpret the results.

#Problem 1d - Can you use MDS to help determine which distance is appropriate for this data?  Which one is best and why?

#Problem 1e - Apply MDS with your chosen distance to visualize the words.  Interpret the results.

##########
#Problem 2 - K-means
##########
#Problem 2a - Apply K-means with K=4 to this data.

#Problem 2b - How well does K-mean do at separating the authors?

#Problem 2c - Is K-means an appropriate clustering algorithm for this data?  Why or Why nor?

#############
#Problem 3 - Hierarchical Clustering
#############
#Problem 3a - Apply hierarchical clustering to this data set.

#Problem 3b - Which distance is best to use?  Why?

#Problem 3c - Which linkage is best to use?  Why?

#Problem 3d - Do any linkages perform particularly poorly?  Explain this result.

#Problem 3e - Visualize your hierarchical clustering results.

###########
#Problem 4 - Biclustering
###########

#Problem 4a - Apply the cluster heatmap method to visualize this data.  Which distance and linkage functions did you use?

#Problem 4b - Interpret the cluster heatmap.  Which words are important for distinguishing author texts?


###########
#Problem 5 - NMF
###########

#Problem 5a - Apply NMF with K = 4 and use W to assign cluster labels to each observation.

#Problem 5b - How well does NMF perform?  Interpret and explain this result.

#Problem 5c - Can you use the NMF to determine which words are important for distinguishing author texts?  How?  What did you find?

#############
#Problem 6 - Wrap-up
############
#Problem  6a - Overall, which method is the best at clustering the author texts?  Why is this the case?

#Problem 6b - Which words are key for distinguishing the author texts?  How did you determine these?

#Problem 6c - Overall, which is the best method for providing a visual summary of the data?

#######################


###############################################################
###############################################################
#R scripts to help out with the Clustering Lab
#Don't peek at this if you want to practice coding on your own!!
##################################################################


#######################
#author data

load("UnsupL.Rdata")

#understand the data a bit
dim(author)
colnames(author)
unique(rownames(author))
TrueAuth = as.factor(rownames(author))

par(mfrow=c(2,2)) 
hist(author[,colnames(author)=="the"],breaks=25)
hist(author[,colnames(author)=="a"],breaks=25)
hist(author[,colnames(author)=="and"],breaks=25)
hist(author[,colnames(author)=="things"],breaks=25)

#take out bookID 
X = author[,1:69] 


#############
#Visulaizing data - how to visulaize texts?  words? in 2-dimensions

#trying PCA
sv = svd(X);
V = sv$v
Z = X%*%V;
plot(Z[,1],Z[,2],type="n")
text(Z[,1],Z[,2],rownames(X),col=as.numeric(TrueAuth),cex=.5)
#why doesn't this work well?

########
#trying MDS (classical)
#can you use MDS to decide which distance is best to understand this data? 

#visualizing author texts
Dmat = dist(X,method="canberra")
mdsres = cmdscale(Dmat,k=2)
plot(mdsres[,1],mdsres[,2],type="n")
text(mdsres[,1],mdsres[,2],rownames(X),col=as.numeric(TrueAuth),cex=.5)

#visulaizing words
Dmat = dist(t(X),method="canberra")
mdsresW = cmdscale(Dmat,k=2)
plot(mdsresW[,1],mdsresW[,2],type="n")
text(mdsresW[,1],mdsresW[,2],colnames(X))

##############
#K- means
K = 4
km = kmeans(X,centers=K)
table(km$cluster,TrueAuth)

plot(mdsres[,1],mdsres[,2],type="n")
text(mdsres[,1],mdsres[,2],rownames(X),col=km$cluster,cex=.5)

###############
#hierarchical clustering
#which distance is appropraite?  Why?  

Dmat = dist(X,method="canberra")
com.hc = hclust(Dmat,method="complete")
res.com = cutree(com.hc,4)
table(res.com,TrueAuth)

plot(com.hc,cex=.5)

#which linkage is best?  Why?

Dmat = dist(X,method="canberra")
com.hc = hclust(Dmat,method="ward.D")
res.com = cutree(com.hc,4)
table(res.com,TrueAuth)

plot(com.hc,cex=.5)

#do any preform terribly?  Why?

#visualize hierarchical clustering reuslts using MDS
table(res.com,TrueAuth)
plot(mdsres[,1],mdsres[,2],type="n")
text(mdsres[,1],mdsres[,2],rownames(X),col=res.com,cex=.5)


#############
#cluster heatmap

heatmap(X,distfun=function(x)dist(x,method="canberra"),hclustfun=function(x)hclust(x,method="ward.D"))

heatmap(scale(X),distfun=function(x)dist(x,method="canberra"),hclustfun=function(x)hclust(x,method="ward.D"),cex=1.5)

#############
#NMF

require("NMF")

K = 4
nmffit = nmf(X,rank=K) 
W = basis(nmffit)
H = coef(nmffit)

cmap = apply(W,1,which.max)
table(cmap,TrueAuth)

par(mfrow=c(1,2))
basismap(nmffit,annRow=rownames(X),scale="col",legend=FALSE)
coefmap(nmffit,annCol=colnames(X),scale="col",legend=FALSE)

par(mfrow=c(1,1))
coefmap(nmffit,annCol=colnames(X),scale="col",legend=FALSE)

#which words are most important for distinguishing authors?

########################################################
