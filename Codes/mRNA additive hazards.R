
#########################################################################
#		Install requiered Packages						#
#########################################################################

install.packages("devtools")

install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("Biobase")
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
BiocManager::install("S4Vectors")
install_github("binderh/CoxBoost")
install_github("fbertran/peperr")
install.packages("locfit")
install.packages("EnrichmentBrowser")
library("devtools")
library("S4Vectors")
library('BiocManager')
library("Biobase")
library("EnrichmentBrowser")
library("matrixStats")
library("Hmisc") 
library("lattice")
library("survival")
library("Formula")
library("ggplot2")
library(utils)
library("Biobase")
#library('CoxBoost')
library('peperr')
library('locfit')
library('Hmisc')
library('rrp')
library(erer)
library(pec)
library(preprocessCore)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cqn")
library(cqn)
library(scales)
library(bestNormalize)

nn <- cqn(t(z),lengthMethod = "fixed")


#########################################################################
#		Change directory								#
#########################################################################


setwd("E:\\students\\Salimi")


#########################################################################
#		read data from directory						#
#########################################################################


X=read.table("mRNA.txt", sep='\t', header=T,fill=T)
dim(X)
X= X[,-1] # remove sample names

X[,1]
X[,2]
head(X)

colnames(X)
rownames(X)
dim(X) ##286 22283

#########################################################################
#		read clinical data 							#
#########################################################################


S.time= X[,1]
S.status = X[,2]
X= X[,-(1:2)] ##mRNAs
dim(X)
#p values and hazard ratios using single variate COxPH

#########################################################################
#		Normalization of data 							#
#########################################################################


#z <- (log2(X))
z <- X

new.z <- matrix(0, ncol = ncol(z), nrow = nrow(z))
dim(new.z)

for (jj in 1:ncol(z)){
bn <- bestNormalize(z[,jj])
#hist(z[,1])
new.z[,jj] <- bn$x.t
}

colnames(new.z) <- colnames(z)

write.table(new.z, file='E:/students/Salimi/normalizedmRNA.txt', sep=',')


cft <- S.time +runif(nrow(z))		# censored failure time
del <- S.status			# failure status

Z=read.table("Payam_normalized_data.txt", sep=',', header=T,fill=T)
dim(Z)
colnames(Z)
#summary(Z[,1])
#summary(Z[,10])
new.z=Z
###################Screening by adjusted univariate P-values
library(addhazard)

n.data <- data.frame(cft,del,new.z)
dim(n.data)
colnames(n.data)
my.pval <- c()

for (kk in 1:ncol(new.z)){
fit3 <- ah(Surv(cft ,del == 1) ~ new.z[,kk] , ties = FALSE, data = n.data, robust = TRUE)
a <- summary(fit3)

my.pval[kk] <- a$coefficients[,6]
}

sig.z <- new.z[,which(my.pval<0.05)]
colnames(sig.z )
write.table(sig.z, file='E:/students/Salimi/reduceddata.txt', sep=',')
reduceddata = sig.z
library(foreign)

reduceddata <- read.table(file.choose(), header = TRUE, dec =".",sep =",")
head(reduceddata)
dim(reduceddata)



#read.table(reduceddata.txt)
#reduceddata
#########################################################################
#		Load R codes for additive hazards					#
#########################################################################


############## Run codes 
##load dll file

dyn.load("D:\\cr.dll","arsim")
is.loaded("arsim")

dyn.load("D:\\ahfun.dll","arsim")
is.loaded("arsim")
source("ahfun.R")

############## 





#########################################################################
#		Analysis on reduced data all 286 patients				#
#########################################################################



z.reduceddata =as.matrix(reduceddata)
################################total data

n.data <- data.frame(cft,del,z.reduceddata)
dim(n.data)
colnames(n.data)
my.pval <- c()

for (kk in 1:ncol(z.reduceddata)){
fit3 <- ah(Surv(cft ,del == 1) ~ z.reduceddata[,kk] , ties = FALSE, data = n.data, robust = TRUE)
a <- summary(fit3)

my.pval[kk] <- a$coefficients[,6]
}

sig.z <- z.reduceddata[,which(my.pval<0.001)]
colnames(sig.z )



write.table(sig.z, file='C:\\Users\\Dear User\\Dropbox\\AminiTapak collaboration\\2 Additive hazards\\New Analysis 1401-2-23\\reduceddata0.001.txt', sep=',')
reduceddata = sig.z
library(foreign)

reduceddata <- read.table(file.choose(), header = TRUE, dec =".",sep =",")
head(reduceddata)
dim(reduceddata)


z.reduceddata = data.matrix(reduceddata)




##variable selection

sica.beta.tot <- matrix(0,ncol = ncol(z.reduceddata), nrow = 100)
dim(sica.beta.tot)
##variable selection
for (kkk in 1:100){
	# solve by coordiante descent
	ans <- cd(cft, del, z.reduceddata, "SICA",nlam=100 )
	sol <- ans$sol; lam <- ans$lam
	#lam <- lam[-c(1,2)]


	# cross-validate
	ans <- cv(cft, del, z.reduceddata, "SICA", lam = lam)
	i <- ans$i; j <- ans$j
	i;j
	#plot(lam, ans$cv.err[,1])

	#sol[, i]
	#bet.hat <- sol[,length(lam)/2 ]
	bet.hat <- sol[,i,j  ]
	#bet.hat
	sica.beta.tot[kkk,] <- bet.hat
	print(kkk)
#which(bet.hat!=0)
}


sica.frec.tot <- matrix(0,ncol=ncol(z.reduceddata), nrow=100)
for(jjj in 1:100){
 sica.frec.tot[jjj,] <- ifelse(sica.beta.tot[jjj,] !=0,1,0)

}


a.sica.tot <- colSums(sica.frec.tot)
a.sica.tot[which(a.sica.tot != 0)]
sica.100.names.tot <- matrix(colnames(z.reduceddata[,which(a.sica.tot!=0)]),ncol=1)
a1=cbind(sica.100.names.tot,a.sica.tot[which(a.sica.tot != 0)])

write.csv(a1, file = "sica.additive.csv")



#matrix(colnames(z.reduceddata[,which(bet.hat!=0)]),ncol=1)



#######SCAD
scad.beta.tot <- matrix(0,ncol = ncol(z.reduceddata), nrow = 100)
dim(scad.beta.tot)
##variable selection
for (kkk in 1:100){
	# solve by coordiante descent
	ans <- cd(cft, del, z.reduceddata, "SCAD",nlam=100 )
	sol <- ans$sol; lam <- ans$lam
	#lam <- lam[-c(1,2)]


	# cross-validate
	ans <- cv(cft, del, z.reduceddata, "SCAD", lam = lam)
	i <- ans$i; j <- ans$j
	i;j
	#plot(lam, ans$cv.err[,1])

	#sol[, i]
	#bet.hat <- sol[,length(lam)/2 ]
	bet.hat <- sol[,i  ]
	#bet.hat
	scad.beta.tot[kkk,] <- bet.hat
	print(kkk)
#which(bet.hat!=0)
}


scad.frec.tot <- matrix(0,ncol=ncol(z.reduceddata), nrow=100)
for(jjj in 1:100){
 scad.frec.tot[jjj,] <- ifelse(scad.beta.tot[jjj,] !=0,1,0)

}


a.scad.tot <- colSums(scad.frec.tot)
a.scad.tot[which(a.scad.tot != 0)]
scad.100.names.tot <- matrix(colnames(z.reduceddata[,which(a.scad.tot!=0)]),ncol=1)
a2=cbind(scad.100.names.tot,a.scad.tot[which(a.scad.tot != 0)])
write.csv(a2, file = "scad.additive.csv")

##############################MCP


#######mcp
mcp.beta.tot <- matrix(0,ncol = ncol(z.reduceddata), nrow = 100)
dim(mcp.beta.tot)
##variable selection
for (kkk in 1:100){
	# solve by coordiante descent
	ans <- cd(cft, del, z.reduceddata, "MCP",nlam=100 )
	sol <- ans$sol; lam <- ans$lam
	#lam <- lam[-c(1,2)]


	# cross-validate
	ans <- cv(cft, del, z.reduceddata, "MCP", lam = lam)
	i <- ans$i; j <- ans$j
	i;j
	#plot(lam, ans$cv.err[,1])

	#sol[, i]
	#bet.hat <- sol[,length(lam)/2 ]
	bet.hat <- sol[,i  ]
	#bet.hat
	mcp.beta.tot[kkk,] <- bet.hat
	print(kkk)
#which(bet.hat!=0)
}


mcp.frec.tot <- matrix(0,ncol=ncol(z.reduceddata), nrow=100)
for(jjj in 1:100){
 mcp.frec.tot[jjj,] <- ifelse(mcp.beta.tot[jjj,] !=0,1,0)

}


a.mcp.tot <- colSums(mcp.frec.tot)
a.mcp.tot[which(a.mcp.tot != 0)]
mcp.100.names.tot <- matrix(colnames(z.reduceddata[,which(a.mcp.tot!=0)]),ncol=1)
a3=cbind(mcp.100.names.tot,a.mcp.tot[which(a.mcp.tot != 0)])
write.csv(a3, file = "mc.additive.csv")



#######Lasso
Lasso.beta.tot <- matrix(0,ncol = ncol(z.reduceddata), nrow = 100)
dim(Lasso.beta.tot)
##variable selection
for (kkk in 1:100){
	# solve by coordiante descent
	ans <- cd(cft, del, z.reduceddata, "Lasso",nlam=100 )
	sol <- ans$sol; lam <- ans$lam
	#lam <- lam[-c(1,2)]
	# cross-validate
	ans <- cv(cft, del, z.reduceddata, "Lasso", lam = lam)
	i <- ans$i; j <- ans$j
	i;j
	#plot(lam,ans$cv.err[,1])

	#sol[, i]
	#bet.hat <- sol[,length(lam)/2 ]
	bet.hat  <- sol[,i  ]
	#bet.hat
	Lasso.beta.tot[kkk,] <- bet.hat
	print(kkk)
#which(bet.hat!=0)
}


Lasso.frec.tot <- matrix(0,ncol=ncol(z.reduceddata), nrow=100)
for(jjj in 1:100){
 Lasso.frec.tot[jjj,] <- ifelse(Lasso.beta.tot[jjj,] !=0,1,0)

}


a.Lasso.tot <- colSums(Lasso.frec.tot)
a.Lasso.tot[which(a.Lasso.tot != 0)]
Lasso.100.names.tot <- matrix(colnames(z.reduceddata[,which(a.Lasso.tot!=0)]),ncol=1)
a4=cbind(Lasso.100.names.tot,a.Lasso.tot[which(a.Lasso.tot != 0)])
write.csv(a4, file = "Lasso.additive.csv")

########################Enet

Enet.beta.tot <- matrix(0,ncol = ncol(z.reduceddata), nrow = 100)
dim(Enet.beta.tot)
##variable selection
for (kkk in 21:100){
	# solve by coordiante descent
	ans <- cd(cft, del, z.reduceddata, "Enet",nlam=100 )
	sol <- ans$sol; lam <- ans$lam
	#lam <- lam[-c(1,2)]


	# cross-validate
	ans <- cv(cft, del, z.reduceddata, "Enet", lam = lam)
	i <- ans$i; j <- ans$j
	i;j
	#plot(lam, ans$cv.err[,1])

	#sol[, i]
	#bet.hat <- sol[,length(lam)/2,j ]
	bet.hat <- sol[,i ,j ]
	#bet.hat
	Enet.beta.tot[kkk,] <- bet.hat
	print(kkk)
#which(bet.hat!=0)
}


Enet.frec.tot <- matrix(0,ncol=ncol(z.reduceddata), nrow=20)
for(jjj in 1:20){
 Enet.frec.tot[jjj,] <- ifelse(Enet.beta.tot[jjj,] !=0,1,0)

}


a.Enet.tot <- colSums(Enet.frec.tot)
a.Enet.tot[which(a.Enet.tot != 0)]
Enet.100.names.tot <- matrix(colnames(z.reduceddata[,which(a.Enet.tot!=0)]),ncol=1)
a5= cbind(Enet.100.names.tot,a.Enet.tot[which(a.Enet.tot != 0)])
write.csv(a5, file = "Enet.additive.csv")
write.table(cbind(Enet.100.names.tot,a.Enet.tot[which(a.Enet.tot != 0)]), file='C:/Users/Dear User/Dropbox/AminiTapak collaboration/CoxBoost/Enet.txt', sep=',')


write.xlsx(final, file='C:/Users/Dear User/Dropbox/AminiTapak collaboration/CoxBoost/selected_genes.xlsx', 



a6 = unique(c(a1[,1],a2[,1],a3[,1],a4[,1],a5[,1]))
write.csv(a6, file = "unique.all.additive.csv")


#######################################################
s1 <- a2
l1 <- a4
si1 <- a1
m1 <- a3


gene <- unique(c(s1[,1],l1[,1],si1[,1],m1[,1]))
final <- data.frame(matrix(0,ncol = 5, nrow = length(gene)))
colnames(final) <- c("gene name","SCAD", "Lasso", "SICA" , "MCP")
final[,1] <- gene
final[which(match(final[,1],s1)!="NA"),2]  <- as.numeric(s1[,2])
final[which(match(final[,1],l1)!="NA"),3]  <- as.numeric(l1[,2])
final[which(match(final[,1],si1)!="NA"),4] <- as.numeric(si1[,2])
final[which(match(final[,1],m1)!="NA"),5]  <- as.numeric(m1[,2])



dim(final)
library("xlsx")
write.table(final, file='E:/students/Salimi/selected_genes.txt', sep="\t")
write.xlsx(final, file='E:/students/Salimi/selected_genes.xlsx', 
sheetName = "Sheet1",col.names = TRUE, row.names = TRUE, append = FALSE)


########################P-value calculation
library(addhazard)
dim(z.reduceddata[,gene])

n.data <- data.frame(cft,del,z.reduceddata[,gene])
zzz <- z.reduceddata[,gene]
dim(n.data)
colnames(n.data)
my.pval <- c()
coef <- matrix(0,ncol = 6, nrow = length(gene))
for (kk in 1:ncol(z.reduceddata[,gene])){
fit3 <- ah(Surv(cft+runif(length(cft)) ,del == 1) ~ zzz[,kk] , ties = FALSE, data = n.data, robust = TRUE)
a <- summary(fit3)

#my.pval[kk] <- a$coefficients[,6]
coef[kk,]<- a$coefficients
}
colnames(coef) <- c("coef","se","lower.95" ,"upper.95","z", "p.value")
head(coef)
my.pval


fit.multi <- ah(Surv(cft+runif(length(cft)) ,del == 1) ~ . , ties = FALSE, data = n.data, robust = TRUE)
sumary.me <- summary(fit.multi)
write.table(my.pval, file='C:/Users/Dear User/Dropbox/AminiTapak collaboration/CoxBoost/pvalue195.txt', sep=',')
write.table(coef, file='E:/students/Salimi/univriate.regression.coef.txt', sep=',')
write.table(sumary.me$coefficients, file='E:/students/Salimi/multivriate.regression.coef.txt', sep=',')




multivariate.bet = sumary.me$coefficient[,1]
PI.sica = multivariate.bet%*% t(n.data[,-c(1,2)])


###########################################################
#Prognostic index                     #
###########################################################

ita <- drop(drop(PI.sica))
hist(ita)

###########################################################
#Risk groups for validation group                         #
###########################################################
rgroup <- ifelse(ita<median(ita),0,1)
table(rgroup )
summary(rgroup )
length(rgroup )



##########################################################
#Survival curves for risk groups                         #
##########################################################

library(survival)
surv_diff <- survdiff(Surv(cft, del) ~ rgroup )
surv_diff
model_fit <- survfit(Surv(S.time, S.status) ~  rgroup)


###################################################################
#           Figure in the manuscript 				      #  
###################################################################
plot(model_fit)
pdf(file = "E:/students/Salimi/risk groups.pdf", 
    # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
# Step 2: Create the plot with R code
plot(model_fit, lty = c("solid", "dashed"), 
col = c("black", "grey"), xlab = "Survival Time in Days", ylab = "Survival Probabilities")

legend("bottomleft", c("Low risk", "High risk"), lty = c("solid", "dashed"), col = c("black", "grey"))

# Step 3: Run dev.off() to create the file!
dev.off()






