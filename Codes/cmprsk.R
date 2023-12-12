
install.packages("cmprsk")

library(readxl)
library(cmprsk)
library(crrp)

mir_ok <- read_excel("E:/mir-ok.xlsx")

View(mir_ok)

dim(mir_ok)

cov.mic = mir_ok[  , 4:453]

(time.mic <- mir_ok[ , 2])

(status.mic <- mir_ok[ ,3])

time.mic <- as.vector(data.matrix(time.mic)) +runif(length(time.mic))

(status.mic <- as.vector(data.matrix(status.mic)))

cov.mic = data.frame(data.matrix(cov.mic))

rep.num = 10

(beta.lasso.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))
(beta.scad.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))
(beta.mcp.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))



data.mic1 = mir_ok 




for (i in 1: rep.num){
  
  
  index.m <- sample(1:nrow(data.mic1),round(0.7*nrow(data.mic1)))
  train.data.mic = data.mic1[index.m,]
  dim(train.data.mic)
  
  test.data.mic =  data.mic1[-index.m,]
  dim(test.data.mic)
  
  
  
  cov1.train.m =  train.data.mic[ , 4:453]
  library(crrp)
  
  cov1.test.m =  test.data.mic[ , 4:453]
  
  (time.train.m <- train.data.mic[ , 2])
  
  (status.train.m <-  train.data.mic[ ,3])
  
  time.train.m <- as.vector(data.matrix(time.train.m)) +runif(length(time.train.m))
  
  (status.train.m <- as.vector(data.matrix(status.train.m))    )
  
  
  (time.test.m <-  test.data.mic[ , 2])
  
  (status.test.m <-  test.data.mic[ ,3])
  
  time.test.m <- as.vector(data.matrix(time.test.m)) +runif(length(time.test.m))
  
  (status.test.m <- as.vector(data.matrix(status.test.m))    )
  
  
  ##lasso
  fit.m.1 <- crrp( time.train.m , status.train.m , cov1.train.m  , penalty="LASSO")
  
  a.m <-  which.min(fit.m.1$BIC)
  
  opt.lambda.m.1 = fit.m.1$lambda[a.m]
  
  fit.m.1 <- crrp( time.test.m , status.test.m , cov1.test.m ,lambda =opt.lambda.m.1 , penalty="LASSO")
  
  beta.lasso.mic[i, ] <- fit.m.1$beta[ , which.min(fit.m.1$BIC)]
  
  
  #scad
  
  fit.m.2 <- crrp( time.train.m , status.train.m , cov1.train.m  , penalty="SCAD")
  
  b.m <-  which.min(fit.m.2$BIC)
  
  opt.lambda.m.2 = fit.m.2$lambda[b.m]
  
  fit.m.2 <- crrp( time.test.m , status.test.m , cov1.test.m ,lambda =opt.lambda.m.2 , penalty="SCAD")
 
  
  beta.scad.mic[i,] <- fit.m.2$beta[, which.min(fit.m.2$BIC)]
  
  
  #mcp
  
  fit.m.3 <- crrp( time.train.m , status.train.m , cov1.train.m  , penalty="MCP")
  
  c.m <-  which.min(fit.m.3$BIC)
  
  opt.lambda.m.3 = fit.m.3$lambda[c.m]
  
  fit.m.3 <- crrp( time.test.m , status.test.m , cov1.test.m ,lambda =opt.lambda.m.3 , penalty="MCP")
  
  beta.scad.mic[i,] <- fit.m.3$beta[, which.min(fit.m.3$BIC)]
  
  
  print(i) }







(freq.lasso.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))
(freq.scad.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))
(freq.mcp.mic <- matrix(0,ncol = ncol(cov.mic), nrow = rep.num))





for(g in 1:rep.num){
  
  freq.lasso.mic[g,] <- ifelse(beta.lasso.mic[g,] !=0,1,0)
  
  freq.scad.mic[g,] <- ifelse(beta.scad.mic[g,] !=0,1,0)
  
  freq.mcp.mic[g,] <- ifelse(beta.mcp.mic[g,] !=0,1,0)
}


colSums(freq.lasso.mic)
colSums(freq.scad.mic)
colSums(freq.mcp.mic)

http://webinar.skums.ac.ir/b/hea-c2k-r42


