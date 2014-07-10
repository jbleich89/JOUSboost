library(randomForest)
library(bartMachine)
library(JOUSboost)
library(ada)
set_bart_machine_memory(2000)

set.seed(10)
n = 200 
x1 = rbinom(n, 1, .5)
x0 = matrix(runif(n*200), nrow = n)
xmat = data.frame(cbind(x1,x0))
ps = ifelse(x1 == 1, .75, .25)
for(i in 1 : n){
  if(x1[i] == 1){
    y[i] = rbinom(1, 1,.75)
  }else{y[i] = rbinom(1, 1,.25)}
}

ytrain = y[1 : (n/2)]
ytest = y[(n/2 + 1) : n]
xtrain = xmat[1 : (n/2), ]
xtest = xmat[(n/2 + 1) : n, ]
ytrjb = ifelse(ytrain == 1, 1, -1)
ytejb = ifelse(ytrain == 1, 1, -1)
ps_test = ps[(n/2 + 1) : n]

jb = JOUSboost(X = xtrain, y = ytrjb, type = "under", num_iter = 200)
phats_jb = predict(jb, xtest)
hist(phats_jb)
summary(phats_jb)

rf = randomForest( y = as.factor(ytrain), x = xtrain, ntree = 1000)
rf$confusion


phats = predict(rf, newdata = xtest, type = "prob")
summary(phats[,2])
hist(phats[,2])

adab = ada(x = xtrain, y = ytrjb, iter = 50, nu = 1)
predict(adab, xtest)

bart = bartMachine(X = xtrain, y = as.factor(ytrain), run_in_sample = F)
bpreds = predict(bart, xtest)
hist(bpreds)

sum((ps_test - phats_jb)^2)/(.5*n)
sum((ps_test - phats[,2])^2)/(.5*n)
sum((ps_test - bpreds)^2)/(.5*n)



cbind(bpreds, phats[,2], ps)