library(rpart)
set.seed(888, kind = NULL)
n<-500
d = 10
coe = 10

x<-matrix(0,n,d)
for (ddd in 1:d){
  x[,ddd]<-rnorm(n)
}

y<-rep(0,n)
for (i in 1:n){
  y[i]<--1+2*(runif(1)> exp( coe*sum(x[i,1:6])*(1-x[i,1]+x[i,2]-x[i,3]+x[i,4]-x[i,5]+
                                                  x[i,6]) )/
                ( 1+exp( coe*sum(x[i,1:6])*(1-x[i,1]+x[i,2]-x[i,3]+x[i,4]-x[i,5]+x[i,6]))))
}

probs<-rep(0,n)
for (i in 1:n){
  probs[i]<- 1 - 1 * exp( coe*sum(x[i,1:6])*(1-x[i,1]+x[i,2]-x[i,3]+x[i,4]-x[i,5]+
                                               x[i,6]) )/
    ( 1+exp( coe*sum(x[i,1:6])*(1-x[i,1]+x[i,2]-x[i,3]+x[i,4]-x[i,5]+x[i,6])))
}


ntest = 2500
xtest<-matrix(0,ntest,d)
for (ddd in 1:d){
  xtest[,ddd]<-rnorm(ntest)
}

ytest<-rep(ntest)
for (i in 1:ntest){
  ytest[i]<--1+2*(runif(1)> exp( coe*sum(xtest[i,1:6])*(1-xtest[i,1]+xtest[i,2]-xtest[i,3]+xtest[i,4]-xtest[i,5]+
                                                          xtest[i,6]) )/
                    ( 1+exp( coe*sum(xtest[i,1:6])*(1-xtest[i,1]+xtest[i,2]-xtest[i,3]+xtest[i,4]-xtest[i,5]+xtest[i,6]))))
}


ptest<-rep(0,ntest)
for (i in 1:ntest){
  ptest[i]<- 1 - 1 * exp( coe*sum(xtest[i,1:6])*(1-xtest[i,1]+xtest[i,2]-xtest[i,3]+xtest[i,4]-xtest[i,5]+
                                                   xtest[i,6]) )/
    ( 1+exp( coe*sum(xtest[i,1:6])*(1-xtest[i,1]+xtest[i,2]-xtest[i,3]+xtest[i,4]-xtest[i,5]+xtest[i,6])))
}

jb = JOUSboost(x, y, type = "over", num_iter = 250)
object = jb
newdata = xtest
phats = predict.JOUSboost(object, newdata)
sum(abs(phats - ptest)^2/2500)
head(ptest)
head(phats)
head(preds)
sum(phats)