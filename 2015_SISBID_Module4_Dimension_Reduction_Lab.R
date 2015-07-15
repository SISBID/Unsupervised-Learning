#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
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

#############
#Problem 2 - NMF.
#############

#Problem 2a - Apply NMF to this data.

#Problem 2b - Which value of K did you use?  Why?  What happens when you slightly change your chosen K?

#Problem 2c - Interpret the archetypes found.  Do any of them accurately reflect the different digits?  Which ones?

#Problem 2d - Plot NMF basis scatterplots of the factors associated with differences between the digits from 2c.  Do these scatterplots well separate the different digits?  Why or why not?

#############
#Problem 3 - ICA.
#############

#Problem 3a - Apply ICA to this data set.

#Problem 3b - Which value of K did you use?  Why?  What happens when you slightly change your chosen K?

#Problem 3c - Interpret the independent image signals found.  Do any other them accurately reflect the different digits?  Which ones?

###############
#Problem 4 - Comparisons.
##############

#Problem 4a - Compare and contrast PCA, NMF, and ICA on this data set.  Which one best separates the different digits?  Which one reveals the most interesting patterns?

#Problem 4b - Overall, which method do you recommend for this data set and why?

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
load("UnsupL.Rdata")

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
i = 4; j = 3;
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


#########Problem 2 - NMF
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
i = 11; j = 15;
par(mfrow=c(1,1))
plot(W[,i],W[,j],type="n")
text(W[,i],W[,j],rownames(dat38),col=rownames(dat38),cex=.7)

##################
#Problem 3 - ICA

K = 20
icafit = fastICA(t(digits),n.comp=K)

#plot independent source signals 
par(mfrow=c(4,5))
for(i in 1:K){
  imagedigit(icafit$S[,i])
}

#################################################################
