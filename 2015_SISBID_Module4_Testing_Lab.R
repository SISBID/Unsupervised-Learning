#############################################################
#2015 SISBID Module 4 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Testing Lab
############################################################

############
#Data set - Prostate Data (Singh et al. 2002). This data set consists of gene
#expression levels for 6033 genes among 102 men. 
# The dataset is available from the R package "sda"

#The ultimate goal of this lab is to perform multiple testing to identify differentially expressed genes.


#############

#############
#Problem 1 - We wish to identify important genes to differetiate cancer or healthy patients. What kind of tests are reasonable?

#Problem 2 - In order to adjust for multiple comparisons, which procedures should one use?

#Problem 3 - Examine the list of genes identified. 

#######################


###############################################################
###############################################################
#R scripts to help out with the Clustering Lab
#Don't peek at this if you want to practice coding on your own!!
##################################################################



########################################################





install.packages("sda")
library(sda)

## import data
data(singh2002)
x = singh2002$x
y = singh2002$y

n1 = sum(y == "healthy")
n2 = length(y) - n1


ps<-NULL
for(i in 1:ncol(x)) ps <- c(ps,
t.test(x[1:n1,i], x[(n1+1):(n1+n2),i])$p.value)

## ordered p-values
names(ps)<-seq(1,ncol(x),1)
p1 =sort (ps)
p1[1:50]
## plot ordered p-values
plot(p1[1:100], pch=rep('*',100),ylim=c(0,0.003), ylab="ordered p-values")

## rejection boundry of Benjamini-Hochberg's procedure
abline(a=0, b=0.1/ncol(x), col="red")

##rejection boundary of Bonferroni at 0.1
abline(a=0.1/ncol(x), b=0, col="blue",lty=5)

cat("Compute the no. rejection by Bonferroni:",
max(which(sort(ps,decreasing=FALSE) < .1/ncol(x))),
                                        fill=TRUE)
#6
cat("Compute the BH FDR Directly:",
max(which(sort(ps,decreasing=FALSE) < .1*(1:ncol(x))/ncol(x))),
                                        fill=TRUE)
#57

arrows(x0 = 61, y0 = 0.00085, x1 = 58, y1 = p1[57], length = 0.1)
text(63.5, 0.00085, labels="imax = 57", cex=.8, pos=4, col="black")
legend("topleft",legend=c("BH's Procedure", "Bonferroni", "Ordered p-values"),
       lty=c(1, 5, NA), col=c("red","blue", "black"), pch = c(NA, NA, '*'))
  
  
### SAM procedure
library(samr)
ylab<-rep(2,length(y))
ylab[y == "healthy"]=1

samfit<-SAM(t(x),ylab,resp.type="Two class unpaired")    

# examine significant gene list 
 
print(samfit)
 # plot results
plot(samfit)      


