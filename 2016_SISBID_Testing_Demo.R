#############################################################
#2016 SISBID Module 5 - Unsupervised Learning
#Genevera I. Allen & Yufeng Liu
#Testing Demos for use in lecture
############################################################



########################################################
#Data set 1 - Simulated Data
#Small simulated data set to demonstrate multiple testing 
#when all null hypthesis hold 
##########################################################

#simulate data
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
ps <- NULL
for(i in 1:1000) ps <- c(ps,
 t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Around 5% of p-values are below 0.05:",
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)

#SAM
library(samr)
ylab<-y+1

samfit<-SAM(x,ylab,resp.type="Two class unpaired", fdr.output=0.1)    

# examine significant gene list 
 
print(samfit)
 # plot results
plot(samfit)      


#Note: function p.adjust returns p-values adjusted using one
#of several methods for a given set of p-values.
# For example, for the Bonferroni correction, p-values are
# multiplied by the number of comparisons.

########################################################
#Data set 2 - Simulated Data
#Small simulated data set to demonstrate multiple testing 
#when not all null hypthesis hold 
##########################################################

#simulate data
x <- matrix(rnorm(1000*50),ncol=50)
y <- sample(c(0,1),50,rep=TRUE)
x[1:100,y==0] <- x[1:100,y==0] + 1
ps <- NULL
for(i in 1:1000) ps <- c(ps,
t.test(x[i,y==0],x[i,y==1])$p.value)
cat("Way more than 5% of p-values are below 0.05:",
mean(ps<.05),fill=TRUE)
fdrs.bh <- p.adjust(ps, method="BH")
plot(ps,fdrs.bh)
plot(fdrs.bh)
cat("Number of Tests with FDR below 0.4:",
sum(fdrs.bh<0.4), fill=TRUE)
cat("Compute the BH FDR Directly:",
max(which(sort(ps,decreasing=FALSE) < .4*(1:1000)/1000)),
                                        fill=TRUE)
plot(sort(ps,decreasing=FALSE),ylab="P-Values")
#BH procedure
abline(a=0, b=0.4/1000,col="red")
#Bonferroni
abline(a=0.4/1000, b=0,col="blue")

#SAM
library(samr)
ylab<-y+1

samfit<-SAM(x,ylab,resp.type="Two class unpaired", fdr.output=0.1)    

# examine significant gene list 
 
print(samfit)
 # plot results
plot(samfit)      


