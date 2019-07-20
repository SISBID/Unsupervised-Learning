#############################################################
#2019 SISBID Module 3 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu & Ali Shojaie 
#Dimension Reduction Lab
############################################################

############
#Data set - Digits Data.
#Either use all digits or choose 2-3 digits if computational speed is a problem. Looking at 3's, 8's and 5's are interesting. Note that NMF takes quite a while to run, so you may want to limit the digits considered for that problem.
############

############
#Problem 1 - PCA
############
#Problem 1a - Apply PCA to this data.

#Problem 1b - Do the first several PCs well separate different digits?  Why or why not?

#Problem 1c - Use the first several PCs and PC loadings to evaluate the major patterns in the digits data.  Can you come up with a description of the pattern found by each of the first five PCs?


#Problem 1d - How many PCs are needed to explain 95% of the variance?  You must decide how many PCs to retain.  Which do you pick and why?


#########
#Problem 2 - MDS
#########

#Problem 2a - Apply MDS (classical or non-metric) to this data.  Try out several distance metrics and different numbers of MDS components.

#Problem 2b - Which distance metric is best for this data?  Which one reveals the most separation of the digits?

#Problem 2c - Compare and contrast the MDS component maps to the dimension reduction of PCA.  Which is preferable?


#############
#Problem 3 - NMF.
#############

#Problem 3a - Apply NMF to this data.

#Problem 3b - Which value of K did you use?  Why?  What happens when you slightly change your chosen K?

#Problem 3c - Interpret the archetypes found.  Do any of them accurately reflect the different digits?  Which ones?

#Problem 3d - Plot NMF basis scatterplots of the factors associated with differences between the digits from 2c.  Do these scatterplots well separate the different digits?  Why or why not?

#############
#Problem 4 - ICA.
#############

#Problem 4a - Apply ICA to this data set.

#Problem 4b - Which value of K did you use?  Why?  What happens when you slightly change your chosen K?

#Problem 4c - Interpret the independent image signals found.  Do any other them accurately reflect the different digits?  Which ones?

###############
#Problem 5 - Comparisons.
##############

#Problem 5a - Compare and contrast PCA, MDS, NMF, and ICA on this data set.  Which one best separates the different digits?  Which one reveals the most interesting patterns?

#Problem 5b - Overall, which method do you recommend for this data set and why?

################################################


##########
#Additional Data set - NCI Microarray data
#(If you have time - take a further look at this data set using various methods for dimension reduction.  Also you may be interested in trying MDS to visualize this data.)
###########



###############################################################
###############################################################
#R scripts to help out with the Dimension Reduction Lab
#Don't peek at this if you want to practice coding on your own!!
##################################################################
#code for digits - ALL
load("UnsupL_SISBID_2019.Rdata")

#visulaize
par(mfrow=c(4,8))
for(i in 1:32){
  imagedigit(digits[i,])
}

########Problem 1 - PCA
#PCA - take SVD to get solution
#don't center and scale to retain interpretation as images
svdd = svd(digits)
U = svdd$u
V = svdd$v #PC loadings
D = svdd$d
Z = digits%*%V #PCs

#PC scatterplot
i = 1; j = 2;
plot(U[,i],U[,j],type="n")
text(U[,i],U[,j],rownames(digits),col=rownames(digits),cex=.7)

#PC loadings
par(mfrow=c(3,5))
for(i in 1:15){
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


###############
#Problem 2 - MDS
#classical MDS
#(Note, this may take some time - try only on 3's and 8's)

dat38 = rbind(digits[which(rownames(digits)==3),],digits[which(rownames(digits)==8),])

#PCA for comparison
svdd = svd(dat38)
U = svdd$u
V = svdd$v #PC loadings
D = svdd$d
Z = digits%*%V #PCs

#MDS
Dmat = dist(dat38,method="manhattan") #Manhattan (L1) Distance
mdsres = cmdscale(Dmat,k=10)

i = 1; j = 2;
par(mfrow=c(1,2))
plot(mdsres[,i],mdsres[,j],type="n")
text(mdsres[,i],mdsres[,j],rownames(dat38),col=rownames(dat38))

plot(U[,i],U[,j],type="n",xlab="PC1",ylab="PC2")
text(U[,i],U[,j],rownames(dat38),col=rownames(dat38))


#########Problem 3 - NMF
#NMF

require("NMF")

dat38 = rbind(digits[which(rownames(digits)==3),],digits[which(rownames(digits)==8),])

K = 20
nmffit = nmf(dat38+1,rank=K) #note - this takes a while
W = basis(nmffit)
H = coef(nmffit)

#plot archetypes
par(mfrow=c(4,5))
for(i in 1:K){
  imagedigit(H[i,])
}

#plot scaterrplots of W for interesting archetypes
i = 1; j = 2;
par(mfrow=c(1,1))
plot(W[,i],W[,j],type="n")
text(W[,i],W[,j],rownames(dat38),col=rownames(dat38),cex=.7)

##################
#Problem 4 - ICA

require("fastICA")

K = 20
icafit = fastICA(t(digits),n.comp=K)

#plot independent source signals 
par(mfrow=c(4,5))
for(i in 1:K){
  imagedigit(icafit$S[,i])
}

#################################################################
