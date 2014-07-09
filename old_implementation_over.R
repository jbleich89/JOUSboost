
###########for laptop with old R do predict.rpart instead of predict and type="vector" instead of type="prob"


library(LogitBoost, lib.loc="temp/LogitBoost")
library(randomForest, lib.loc="temp")
library(class)

dep<-3

logitboost<-function (xlearn, ylearn, xtest, mfinal, presel = 0, estimate = 0, 
                      verbose = FALSE) 
{
  if (nlevels(as.factor(ylearn)) == 2) {
    if (presel > 0) {
      s <- apply(xlearn, 2, score, ylearn)
      quality <- apply(rbind(s, -s + (sum(ylearn == 0) * 
                                        sum(ylearn == 1))), 2, max)
      genes <- rev(order(quality))[1:presel]
      xlearn <- xlearn[, genes]
      xtest <- xtest[, genes, drop = FALSE]
    }
    if (estimate > 0) {
      if (verbose) {
        print("Stopping Parameter Estimation")
      }
      likeli <- numeric(mfinal)
      probabs <- crossval(xlearn, ylearn, estimate, mfinal)$probs
      for (k in 1:mfinal) {
        a <- pmax(log(probabs[, k]), -1e+36)
        b <- pmax(log(1 - probabs[, k]), -1e+36)
        likeli[k] <- (ylearn %*% a) + ((1 - ylearn) %*% 
                                         b)
      }
    }
    learn <- dim(xlearn)[1]
    print(learn)
    test <- dim(xtest)[1]
    Flearn <- numeric(learn)
    Ftest <- numeric(test)
    flearn <- numeric(learn)
    ftest <- numeric(test)
    z <- numeric(learn)
    w <- numeric(learn)
    plearn <- rep(1/2, learn)
    ptest <- matrix(0, test, mfinal)
    if (verbose) {
      print("Boosting Iterations")
    }
    for (m in 1:mfinal) {
      w <- pmax(plearn * (1 - plearn), 1e-24)
      z <- (ylearn - plearn)/w
      cntrl <- rpart.control(maxdepth = dep, minsplit = 0,
                             maxcompete = 0, maxsurrogate = 0, cp = 0, 
                             xval = 0)
      xx <- xlearn
      fit <- rpart(z ~ xx, weights = w/mean(w), control = cntrl)
      #print(fit)
      flearn <- predict(fit)
      xx <- xtest
      ftest <- predict(fit, newdata = data.frame(xx))
      Flearn <- Flearn + (1/2) * flearn
      Ftest <- Ftest + (1/2) * ftest
      plearn <- 1/(1 + exp((-2) * Flearn))
      ptest[, m] <- 1/(1 + exp((-2) * Ftest))
    }
    output <- list(probs = ptest)
    if (estimate > 0) {
      output <- list(probs = ptest, loglikeli = matrix(likeli, 
                                                       nr = test, nc = mfinal, byrow = TRUE))
    }
  }
  if (nlevels(as.factor(ylearn)) > 2) {
    K <- nlevels(as.factor(ylearn))
    likeli <- array(0, c(dim(xtest)[1], mfinal, K))
    ptest <- array(0, c(dim(xtest)[1], mfinal, K))
    for (k in 0:(K - 1)) {
      yyl <- as.numeric(ylearn == k)
      if (presel > 0) {
        s <- apply(xlearn, 2, score, yyl)
        quality <- apply(rbind(s, -s + (sum(yyl == 0) * 
                                          sum(yyl == 1))), 2, max)
        genes <- rev(order(quality))[1:presel]
        xxl <- xlearn[, genes]
        xxt <- xtest[, genes, drop = FALSE]
      }
      if (estimate > 0) {
        if (verbose) {
          print("Stopping Parameter Estimation")
        }
        probabs <- crossval(xxl, yyl, estimate, mfinal)$probs
        for (i in 1:mfinal) {
          a <- pmax(log(probabs[, i]), -1e+36)
          b <- pmax(log(1 - probabs[, i]), -1e+36)
          for (q in 1:dim(xtest)[1]) {
            likeli[q, i, (k + 1)] <- (yyl %*% a) + ((1 - 
                                                       yyl) %*% b)
          }
        }
      }
      learn <- dim(xxl)[1]
      test <- dim(xxt)[1]
      Flearn <- numeric(learn)
      Ftest <- numeric(test)
      flearn <- numeric(learn)
      ftest <- numeric(test)
      z <- numeric(learn)
      w <- numeric(learn)
      plearn <- rep(1/2, learn)
      if (verbose) {
        print("Boosting Iterations")
      }
      for (m in 1:mfinal) {
        w <- pmax(plearn * (1 - plearn), 1e-24)
        z <- (yyl - plearn)/w
        cntrl <- rpart.control(maxdepth = 1, minsplit = learn - 
                                 1, xval = 0, maxcompete = 0, cp = 0, maxsurrogate = 0, 
        )
        xx <- xxl
        fit <- rpart(z ~ xx, weights = w/mean(w), control = cntrl)
        flearn <- predict(fit)
        xx <- xxt
        ftest <- predict(fit, newdata = data.frame(xx))
        Flearn <- Flearn + (1/2) * flearn
        Ftest <- Ftest + (1/2) * ftest
        plearn <- 1/(1 + exp((-2) * Flearn))
        ptest[, m, (k + 1)] <- 1/(1 + exp((-2) * Ftest))
      }
    }
    output <- list(probs = ptest)
    if (estimate > 0) {
      output <- list(probs = ptest, loglikeli = likeli)
    }
  }
  output
}








#set max iterations

iterations<-2

stop<-iterations

error<-rep(0,iterations)
undererror<-rep(0,iterations)
linkerror<-rep(0,iterations)
lblinkerror<-rep(0,iterations)

exploss<-rep(0,iterations)
underexploss<-rep(0,iterations)
linkexploss<-rep(0,iterations)
lblinkexploss<-rep(0,iterations)

logerror<-rep(0,iterations)
underlogerror<-rep(0,iterations)
loglinkerror<-rep(0,iterations)
lbloglinkerror<-rep(0,iterations)

misclass<-rep(0,iterations)
lbmisclass<-rep(0,iterations)

f90forplot<-rep(0,iterations)
f70forplot<-rep(0,iterations)
f50forplot<-rep(0,iterations)

lbf90forplot<-rep(0,iterations)
lbf70forplot<-rep(0,iterations)
lbf50forplot<-rep(0,iterations)


iter<-seq(1,iterations)



set.seed(1234, kind = NULL)


#####hold out data is largex and largey

number<-50
m<-number^2
d<-1000

largex<-matrix(0,m,d)
for (k in 1:d) {
  largex[,k]<-rnorm(m)
}




largey<-rep(0,number^2)
realp<-rep(0,number^2)




largey<-rep(0,(number^2))
for (i in 1:(number^2)){
  largey[i]<--1+2*(runif(1)>.03+.94*(sum(largex[i,1:d])>0))
}



realp<-rep(0,(number^2))
for (i in 1:(number^2)){
  realp[i]<-1-(.03+.94*(sum(largex[i,1:d])>0))
}




#####training data is x and y #############


n<-5000


realx<-matrix(0,n,d)
for (ddd in 1:d){
  realx[,ddd]<-rnorm(n)
}


x<-realx 



realy<-rep(0,n)
for (i in 1:n){
  realy[i]<--1+2*(runif(1)>.03+.94*(sum(x[i,1:d])>0))
}




y<-realy

############################################



#fit<-knn(x,largex,y,k=1,prob=T)
#fit<-as.numeric(fit)*2-3
#print("nn1 error")
#nnerror<-(1-sum((fit==largey)/length(largey)))
#print(nnerror)




noise<-2



ss<-sum(realy==1) ##number yes
ff<-sum(realy==-1)  ##number no

xplus1<-realx[realy==1,] ## matrix for 1s
xminus1<-realx[realy==-1,] ##matrix for -1s

yplus1<-realy[realy==1] ##the 1s
yminus1<-realy[realy==-1] ##the -1s

cutplus1<-round(ss/9,0) ##?? --divide into equal pieces
cutminus1<-round(ff/9,0)  ##?? -- divide into equal pieces


underx<-realx
underx10<-rbind(xplus1[1:(1*cutplus1),],xminus1[1:(ff),]) ##1 parts 1, 9 part - 1
underx20<-rbind(xplus1[1:(2*cutplus1),],xminus1[1:(8*cutminus1),]) ##2 vs 8
underx30<-rbind(xplus1[1:(3*cutplus1),],xminus1[1:(7*cutminus1),])
underx40<-rbind(xplus1[1:(4*cutplus1),],xminus1[1:(6*cutminus1),])
underx60<-rbind(xplus1[1:(6*cutplus1),],xminus1[1:(4*cutminus1),])
underx70<-rbind(xplus1[1:(7*cutplus1),],xminus1[1:(3*cutminus1),])
underx80<-rbind(xplus1[1:(8*cutplus1),],xminus1[1:(2*cutminus1),])
underx90<-rbind(xplus1[1:(ss),],xminus1[1:(1*cutminus1),])

undery<-realy
undery10<-c(yplus1[1:(1*cutplus1)],yminus1[1:(ff)]) ##1 parts 1, 9 part - 1
undery20<-c(yplus1[1:(2*cutplus1)],yminus1[1:(8*cutminus1)]) ##2 vs 8
undery30<-c(yplus1[1:(3*cutplus1)],yminus1[1:(7*cutminus1)])
undery40<-c(yplus1[1:(4*cutplus1)],yminus1[1:(6*cutminus1)])
undery60<-c(yplus1[1:(6*cutplus1)],yminus1[1:(4*cutminus1)])
undery70<-c(yplus1[1:(7*cutplus1)],yminus1[1:(3*cutminus1)])
undery80<-c(yplus1[1:(8*cutplus1)],yminus1[1:(2*cutminus1)])
undery90<-c(yplus1[1:(ss)],yminus1[1:(1*cutminus1)])




undern<-length(underx[,1])
undern10<-length(underx10[,1])
undern20<-length(underx20[,1])
undern30<-length(underx30[,1])
undern40<-length(underx40[,1])
undern60<-length(underx60[,1])
undern70<-length(underx70[,1])
undern80<-length(underx80[,1])
undern90<-length(underx90[,1])


underf<-rep(0,undern)
underf10<-rep(0,undern10)
underf20<-rep(0,undern20)
underf30<-rep(0,undern30)
underf40<-rep(0,undern40)
underf60<-rep(0,undern60)
underf70<-rep(0,undern70)
underf80<-rep(0,undern80)
underf90<-rep(0,undern90)

underlargef<-0
underlargef10<-0
underlargef20<-0
underlargef30<-0
underlargef40<-0
underlargef60<-0
underlargef70<-0
underlargef80<-0
underlargef90<-0

###############################





bayespred<-rep(0,m)
for (i in 1:(m)){
  bayespred[i]<--1+2*(sum(largex[i,1:d])<0)
}

print("bayes error")
bayeserror<-1-sum((as.numeric(bayespred)==largey))/length(bayespred)
print(bayeserror)


################################################################




#######################################################


jitters<-matrix(0,(8*ss),d) ##need to only jitter n-1 of the repeats
jitterf<-matrix(0,(8*ff),d) ##need to only jitter up to n-1 of the repeats


for (ghg in 1:d){
  jitters[,ghg]<-sqrt(var(realx[,ghg]))*noise*(runif(8*ss)-.5) ##beteween - noise*var and noise*var
}


for (ghg in 1:d){
  jitterf[,ghg]<-sqrt(var(realx[,ghg]))*noise*(runif(8*ff)-.5) 
}



runx<-realx
runy<-realy



runx10<-realx
runy10<-realy
for(k in 1:8){
  runx10<-rbind(runx10,(realx[realy==-1,]))
  runy10<-c(runy10,realy[realy==-1]     )   
}
for (ghg in 1:d){
  runx10[,ghg]<-runx10[,ghg]+c(rep(0,n),jitterf[1:(8*ff),ghg])
}




runx20<-realx
runy20<-realy
for(k in 1:7){
  runx20<-rbind(runx20,(realx[realy==-1,]))
  runy20<-c(runy20,realy[realy==-1]     )   
}
for(k in 1:1){
  runx20<-rbind(runx20,(realx[realy==1,]))
  runy20<-c(runy20,realy[realy==1]      )  
}
for (ghg in 1:d){
  runx20[,ghg]<-runx20[,ghg]+c(rep(0,n),jitterf[1:(7*ff),ghg],jitters[1:(1*ss),ghg])
}


runx30<-realx
runy30<-realy
for(k in 1:6){
  runx30<-rbind(runx30,(realx[realy==-1,]))
  runy30<-c(runy30,realy[realy==-1]     )   
}
for(k in 1:2){
  runx30<-rbind(runx30,(realx[realy==1,]))
  runy30<-c(runy30,realy[realy==1]      )  
}
for (ghg in 1:d){
  runx30[,ghg]<-runx30[,ghg]+c(rep(0,n),jitterf[1:(6*ff),ghg],jitters[1:(2*ss),ghg])
}



runx40<-realx
runy40<-realy
for(k in 1:5){
  runx40<-rbind(runx40,(realx[realy==-1,]))
  runy40<-c(runy40,realy[realy==-1]     )   
}
for(k in 1:3){
  runx40<-rbind(runx40,(realx[realy==1,]))
  runy40<-c(runy40,realy[realy==1]      )  
}

for (ghg in 1:d){
  runx40[,ghg]<-runx40[,ghg]+c(rep(0,n),jitterf[1:(5*ff),ghg],jitters[1:(3*ss),ghg])
}


runx60<-realx
runy60<-realy
for(k in 1:3){
  runx60<-rbind(runx60,(realx[realy==-1,]))
  runy60<-c(runy60,realy[realy==-1]     )   
}
for(k in 1:5){
  runx60<-rbind(runx60,(realx[realy==1,]))
  runy60<-c(runy60,realy[realy==1]      )  
}
for (ghg in 1:d){
  runx60[,ghg]<-runx60[,ghg]+c(rep(0,n),jitterf[1:(3*ff),ghg],jitters[1:(5*ss),ghg])
}


runx70<-realx
runy70<-realy
for(k in 1:2){
  runx70<-rbind(runx70,(realx[realy==-1,]))
  runy70<-c(runy70,realy[realy==-1]     )   
}
for(k in 1:6){
  runx70<-rbind(runx70,(realx[realy==1,]))
  runy70<-c(runy70,realy[realy==1]      )  
}
for (ghg in 1:d){
  runx70[,ghg]<-runx70[,ghg]+c(rep(0,n),jitterf[1:(2*ff),ghg],jitters[1:(6*ss),ghg])
}



runx80<-realx
runy80<-realy
for(k in 1:1){
  runx80<-rbind(runx80,(realx[realy==-1,]))
  runy80<-c(runy80,realy[realy==-1]     )   
}
for(k in 1:7){
  runx80<-rbind(runx80,(realx[realy==1,]))
  runy80<-c(runy80,realy[realy==1]      )  
}
for (ghg in 1:d){
  runx80[,ghg]<-runx80[,ghg]+c(rep(0,n),jitterf[1:(1*ff),ghg],jitters[1:(7*ss),ghg])
}


runx90<-realx
runy90<-realy
for(k in 1:8){
  runx90<-rbind(runx90,(realx[realy==1,]))
  runy90<-c(runy90,realy[realy==1]     )   
}
for (ghg in 1:d){
  runx90[,ghg]<-runx90[,ghg]+c(rep(0,n),jitters[1:(8*ss),ghg])
}






x<-runx
x10<-runx10
x20<-runx20
x30<-runx30
x40<-runx40
x60<-runx60
x70<-runx70
x80<-runx80
x90<-runx90

y<-runy
y10<-runy10
y20<-runy20
y30<-runy30
y40<-runy40
y60<-runy60
y70<-runy70
y80<-runy80
y90<-runy90






n<-length(x[,1])
n10<-length(x10[,1])
n20<-length(x20[,1])
n30<-length(x30[,1])
n40<-length(x40[,1])
n60<-length(x60[,1])
n70<-length(x70[,1])
n80<-length(x80[,1])
n90<-length(x90[,1])


library(rpart)


f<-rep(0,n)
f10<-rep(0,n10)
f20<-rep(0,n20)
f30<-rep(0,n30)
f40<-rep(0,n40)
f60<-rep(0,n60)
f70<-rep(0,n70)
f80<-rep(0,n80)
f90<-rep(0,n90)

largef<-0
largef10<-0
largef20<-0
largef30<-0
largef40<-0
largef60<-0
largef70<-0
largef80<-0
largef90<-0

i<-1
while(i<=iterations){
  
  
  
  
  w<-exp(-y*f)
  w10<-exp(-y10*f10)
  w20<-exp(-y20*f20)
  w30<-exp(-y30*f30)
  w40<-exp(-y40*f40)
  w60<-exp(-y60*f60)
  w70<-exp(-y70*f70)
  w80<-exp(-y80*f80)
  w90<-exp(-y90*f90)
  
  
  w<-w/max(w)
  w10<-w10/max(w10)
  w20<-w20/max(w20)
  w30<-w30/max(w30)
  w40<-w40/max(w40)
  w60<-w60/max(w60)
  w70<-w70/max(w70)
  w80<-w80/max(w80)
  w90<-w90/max(w90)
  
  w<-w*1.000001
  w10<-w10*1.000001
  w20<-w20*1.000001
  w30<-w30*1.000001
  w40<-w40*1.000001
  w60<-w60*1.000001
  w70<-w70*1.000001
  w80<-w80*1.000001
  w90<-w90*1.000001
  
  
  
  
  underw<-exp(-undery*underf)
  underw10<-exp(-undery10*underf10)
  underw20<-exp(-undery20*underf20)
  underw30<-exp(-undery30*underf30)
  underw40<-exp(-undery40*underf40)
  underw60<-exp(-undery60*underf60)
  underw70<-exp(-undery70*underf70)
  underw80<-exp(-undery80*underf80)
  underw90<-exp(-undery90*underf90)
  
  
  underw<-underw/max(underw)
  underw10<-underw10/max(underw10)
  underw20<-underw20/max(underw20)
  underw30<-underw30/max(underw30)
  underw40<-underw40/max(underw40)
  underw60<-underw60/max(underw60)
  underw70<-underw70/max(underw70)
  underw80<-underw80/max(underw80)
  underw90<-underw90/max(underw90)
  
  underw<-underw*1.000001
  underw10<-underw10*1.000001
  underw20<-underw20*1.000001
  underw30<-underw30*1.000001
  underw40<-underw40*1.000001
  underw60<-underw60*1.000001
  underw70<-underw70*1.000001
  underw80<-underw80*1.000001
  underw90<-underw90*1.000001
  
  
  treee<-rpart(y~.,data.frame(x),w,method="class",
               control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee10<-rpart(y10~.,data.frame(x10),w10,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee20<-rpart(y20~.,data.frame(x20),w20,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee30<-rpart(y30~.,data.frame(x30),w30,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee40<-rpart(y40~.,data.frame(x40),w40,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee60<-rpart(y60~.,data.frame(x60),w60,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee70<-rpart(y70~.,data.frame(x70),w70,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee80<-rpart(y80~.,data.frame(x80),w80,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  treee90<-rpart(y90~.,data.frame(x90),w90,method="class",
                 control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  
  
  
  undertreee<-rpart(undery~.,data.frame(underx),underw,method="class",
                    control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee10<-rpart(undery10~.,data.frame(underx10),underw10,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee20<-rpart(undery20~.,data.frame(underx20),underw20,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee30<-rpart(undery30~.,data.frame(underx30),underw30,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee40<-rpart(undery40~.,data.frame(underx40),underw40,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee60<-rpart(undery60~.,data.frame(underx60),underw60,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee70<-rpart(undery70~.,data.frame(underx70),underw70,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee80<-rpart(undery80~.,data.frame(underx80),underw80,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  undertreee90<-rpart(undery90~.,data.frame(underx90),underw90,method="class",
                      control=rpart.control(minsplit=0,minbucket=1,cp=-1,maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  
  
  gmin<-(-1)+2*(predict(treee,data.frame(x),type="prob")[,1]>
                  predict(treee,data.frame(x),type="prob")[,2])
  gmin10<-(-1)+2*(predict(treee10,data.frame(x10),type="prob")[,1]>
                    predict(treee10,data.frame(x10),type="prob")[,2])
  gmin20<-(-1)+2*(predict(treee20,data.frame(x20),type="prob")[,1]>
                    predict(treee20,data.frame(x20),type="prob")[,2])
  gmin30<-(-1)+2*(predict(treee30,data.frame(x30),type="prob")[,1]>
                    predict(treee30,data.frame(x30),type="prob")[,2])
  gmin40<-(-1)+2*(predict(treee40,data.frame(x40),type="prob")[,1]>
                    predict(treee40,data.frame(x40),type="prob")[,2])
  gmin60<-(-1)+2*(predict(treee60,data.frame(x60),type="prob")[,1]>
                    predict(treee60,data.frame(x60),type="prob")[,2])
  gmin70<-(-1)+2*(predict(treee70,data.frame(x70),type="prob")[,1]>
                    predict(treee70,data.frame(x70),type="prob")[,2])
  gmin80<-(-1)+2*(predict(treee80,data.frame(x80),type="prob")[,1]>
                    predict(treee80,data.frame(x80),type="prob")[,2])
  gmin90<-(-1)+2*(predict(treee90,data.frame(x90),type="prob")[,1]>
                    predict(treee90,data.frame(x90),type="prob")[,2])
  
  
  
  
  
  undergmin<-(-1)+2*(predict(undertreee,data.frame(underx),type="prob")[,1]>
                       predict(undertreee,data.frame(underx),type="prob")[,2])
  undergmin10<-(-1)+2*(predict(undertreee10,data.frame(underx10),type="prob")[,1]>
                         predict(undertreee10,data.frame(underx10),type="prob")[,2])
  undergmin20<-(-1)+2*(predict(undertreee20,data.frame(underx20),type="prob")[,1]>
                         predict(undertreee20,data.frame(underx20),type="prob")[,2])
  undergmin30<-(-1)+2*(predict(undertreee30,data.frame(underx30),type="prob")[,1]>
                         predict(undertreee30,data.frame(underx30),type="prob")[,2])
  undergmin40<-(-1)+2*(predict(undertreee40,data.frame(underx40),type="prob")[,1]>
                         predict(undertreee40,data.frame(underx40),type="prob")[,2])
  undergmin60<-(-1)+2*(predict(undertreee60,data.frame(underx60),type="prob")[,1]>
                         predict(undertreee60,data.frame(underx60),type="prob")[,2])
  undergmin70<-(-1)+2*(predict(undertreee70,data.frame(underx70),type="prob")[,1]>
                         predict(undertreee70,data.frame(underx70),type="prob")[,2])
  undergmin80<-(-1)+2*(predict(undertreee80,data.frame(underx80),type="prob")[,1]>
                         predict(undertreee80,data.frame(underx80),type="prob")[,2])
  undergmin90<-(-1)+2*(predict(undertreee90,data.frame(underx90),type="prob")[,1]>
                         predict(undertreee90,data.frame(underx90),type="prob")[,2])
  
  
  
  largegmin<-(-1)+2*(predict(treee,data.frame(largex),type="prob")[,1]>
                       predict(treee,data.frame(largex),type="prob")[,2])
  largegmin10<-(-1)+2*(predict(treee10,data.frame(largex),type="prob")[,1]>
                         predict(treee10,data.frame(largex),type="prob")[,2])
  largegmin20<-(-1)+2*(predict(treee20,data.frame(largex),type="prob")[,1]>
                         predict(treee20,data.frame(largex),type="prob")[,2])
  largegmin30<-(-1)+2*(predict(treee30,data.frame(largex),type="prob")[,1]>
                         predict(treee30,data.frame(largex),type="prob")[,2])
  largegmin40<-(-1)+2*(predict(treee40,data.frame(largex),type="prob")[,1]>
                         predict(treee40,data.frame(largex),type="prob")[,2])
  largegmin60<-(-1)+2*(predict(treee60,data.frame(largex),type="prob")[,1]>
                         predict(treee60,data.frame(largex),type="prob")[,2])
  largegmin70<-(-1)+2*(predict(treee70,data.frame(largex),type="prob")[,1]>
                         predict(treee70,data.frame(largex),type="prob")[,2])
  largegmin80<-(-1)+2*(predict(treee80,data.frame(largex),type="prob")[,1]>
                         predict(treee80,data.frame(largex),type="prob")[,2])
  largegmin90<-(-1)+2*(predict(treee90,data.frame(largex),type="prob")[,1]>
                         predict(treee90,data.frame(largex),type="prob")[,2])
  
  
  
  
  underlargegmin<-(-1)+2*(predict(undertreee,data.frame(largex),type="prob")[,1]>
                            predict(undertreee,data.frame(largex),type="prob")[,2])
  underlargegmin10<-(-1)+2*(predict(undertreee10,data.frame(largex),type="prob")[,1]>
                              predict(undertreee10,data.frame(largex),type="prob")[,2])
  underlargegmin20<-(-1)+2*(predict(undertreee20,data.frame(largex),type="prob")[,1]>
                              predict(undertreee20,data.frame(largex),type="prob")[,2])
  underlargegmin30<-(-1)+2*(predict(undertreee30,data.frame(largex),type="prob")[,1]>
                              predict(undertreee30,data.frame(largex),type="prob")[,2])
  underlargegmin40<-(-1)+2*(predict(undertreee40,data.frame(largex),type="prob")[,1]>
                              predict(undertreee40,data.frame(largex),type="prob")[,2])
  underlargegmin60<-(-1)+2*(predict(undertreee60,data.frame(largex),type="prob")[,1]>
                              predict(undertreee60,data.frame(largex),type="prob")[,2])
  underlargegmin70<-(-1)+2*(predict(undertreee70,data.frame(largex),type="prob")[,1]>
                              predict(undertreee70,data.frame(largex),type="prob")[,2])
  underlargegmin80<-(-1)+2*(predict(undertreee80,data.frame(largex),type="prob")[,1]>
                              predict(undertreee80,data.frame(largex),type="prob")[,2])
  underlargegmin90<-(-1)+2*(predict(undertreee90,data.frame(largex),type="prob")[,1]>
                              predict(undertreee90,data.frame(largex),type="prob")[,2])
  
  
  
  
  e<-sum(exp(-y*f)*(1-(y==gmin)))
  e10<-sum(exp(-y10*f10)*(1-(y10==gmin10)))
  e20<-sum(exp(-y20*f20)*(1-(y20==gmin20)))
  e30<-sum(exp(-y30*f30)*(1-(y30==gmin30)))
  e40<-sum(exp(-y40*f40)*(1-(y40==gmin40)))
  e60<-sum(exp(-y60*f60)*(1-(y60==gmin60)))
  e70<-sum(exp(-y70*f70)*(1-(y70==gmin70)))
  e80<-sum(exp(-y80*f80)*(1-(y80==gmin80)))
  e90<-sum(exp(-y90*f90)*(1-(y90==gmin90)))
  
  
  
  
  undere<-sum(exp(-undery*underf)*(1-(undery==undergmin)))
  undere10<-sum(exp(-undery10*underf10)*(1-(undery10==undergmin10)))
  undere20<-sum(exp(-undery20*underf20)*(1-(undery20==undergmin20)))
  undere30<-sum(exp(-undery30*underf30)*(1-(undery30==undergmin30)))
  undere40<-sum(exp(-undery40*underf40)*(1-(undery40==undergmin40)))
  undere60<-sum(exp(-undery60*underf60)*(1-(undery60==undergmin60)))
  undere70<-sum(exp(-undery70*underf70)*(1-(undery70==undergmin70)))
  undere80<-sum(exp(-undery80*underf80)*(1-(undery80==undergmin80)))
  undere90<-sum(exp(-undery90*underf90)*(1-(undery90==undergmin90)))
  
  
  
  
  t<-sum(exp(-y*f))
  t10<-sum(exp(-y10*f10))
  t20<-sum(exp(-y20*f20))
  t30<-sum(exp(-y30*f30))
  t40<-sum(exp(-y40*f40))
  t60<-sum(exp(-y60*f60))
  t70<-sum(exp(-y70*f70))
  t80<-sum(exp(-y80*f80))
  t90<-sum(exp(-y90*f90))
  
  
  
  
  undert<-sum(exp(-undery*f))
  undert10<-sum(exp(-undery10*underf10))
  undert20<-sum(exp(-undery20*underf20))
  undert30<-sum(exp(-undery30*underf30))
  undert40<-sum(exp(-undery40*underf40))
  undert60<-sum(exp(-undery60*underf60))
  undert70<-sum(exp(-undery70*underf70))
  undert80<-sum(exp(-undery80*underf80))
  undert90<-sum(exp(-undery90*underf90))
  
  
  
  
  alpha<-.5*log ( (t-e) / e )
  alpha10<-.5*log ( (t10-e10) / e10 )
  alpha20<-.5*log ( (t20-e20) / e20 )
  alpha30<-.5*log ( (t30-e30) / e30 )
  alpha40<-.5*log ( (t40-e40) / e40 )
  alpha60<-.5*log ( (t60-e60) / e60 )
  alpha70<-.5*log ( (t70-e70) / e70 )
  alpha80<-.5*log ( (t80-e80) / e80 )
  alpha90<-.5*log ( (t90-e90) / e90 )
  
  
  
  underalpha<-.5*log ( (undert-undere) / undere )
  underalpha10<-.5*log ( (undert10-undere10) / undere10 )
  underalpha20<-.5*log ( (undert20-undere20) / undere20 )
  underalpha30<-.5*log ( (undert30-undere30) / undere30 )
  underalpha40<-.5*log ( (undert40-undere40) / undere40 )
  underalpha60<-.5*log ( (undert60-undere60) / undere60 )
  underalpha70<-.5*log ( (undert70-undere70) / undere70 )
  underalpha80<-.5*log ( (undert80-undere80) / undere80 )
  underalpha90<-.5*log ( (undert90-undere90) / undere90 )
  
  
  
  
  largef<-largef+alpha*largegmin
  largef10<-largef10+alpha10*largegmin10
  largef20<-largef20+alpha20*largegmin20
  largef30<-largef30+alpha30*largegmin30
  largef40<-largef40+alpha40*largegmin40
  largef60<-largef60+alpha60*largegmin60
  largef70<-largef70+alpha70*largegmin70
  largef80<-largef80+alpha80*largegmin80
  largef90<-largef90+alpha90*largegmin90
  
  
  
  
  
  underlargef<-underlargef+underalpha*underlargegmin
  underlargef10<-underlargef10+underalpha10*underlargegmin10
  underlargef20<-underlargef20+underalpha20*underlargegmin20
  underlargef30<-underlargef30+underalpha30*underlargegmin30
  underlargef40<-underlargef40+underalpha40*underlargegmin40
  underlargef60<-underlargef60+underalpha60*underlargegmin60
  underlargef70<-underlargef70+underalpha70*underlargegmin70
  underlargef80<-underlargef80+underalpha80*underlargegmin80
  underlargef90<-underlargef90+underalpha90*underlargegmin90
  
  
  
  f<-f+alpha*gmin
  f10<-f10+alpha10*gmin10
  f20<-f20+alpha20*gmin20
  f30<-f30+alpha30*gmin30
  f40<-f40+alpha40*gmin40
  f60<-f60+alpha60*gmin60
  f70<-f70+alpha70*gmin70
  f80<-f80+alpha80*gmin80
  f90<-f90+alpha90*gmin90
  
  
  
  
  underf<-underf+underalpha*undergmin
  underf10<-underf10+underalpha10*undergmin10
  underf20<-underf20+underalpha20*undergmin20
  underf30<-underf30+underalpha30*undergmin30
  underf40<-underf40+underalpha40*undergmin40
  underf60<-underf60+underalpha60*undergmin60
  underf70<-underf70+underalpha70*undergmin70
  underf80<-underf80+underalpha80*undergmin80
  underf90<-underf90+underalpha90*undergmin90
  
  
  phat<-(.45)+.10*(largef>0)
  
  phat<-phat+.10*(largef>0)*(largef40>0)
  phat<-phat+.10*(largef>0)*(largef40>0)*(largef30>0)
  phat<-phat+.10*(largef>0)*(largef40>0)*(largef30>0)*(largef20>0)
  phat<-phat+.10*(largef>0)*(largef40>0)*(largef30>0)*(largef20>0)*(largef10>0)
  
  phat<-phat-.10*(largef<0)*(largef60<0)
  phat<-phat-.10*(largef<0)*(largef60<0)*(largef70<0)
  phat<-phat-.10*(largef<0)*(largef60<0)*(largef70<0)*(largef80<0)
  phat<-phat-.10*(largef<0)*(largef60<0)*(largef70<0)*(largef80<0)*(largef90<0)
  
  
  
  
  underphat<-(.45)+.10*(underlargef>0)
  
  underphat<-underphat+.10*(underlargef>0)*(underlargef40>0)
  underphat<-underphat+.10*(underlargef>0)*(underlargef40>0)*(underlargef30>0)
  underphat<-underphat+.10*(underlargef>0)*(underlargef40>0)*(underlargef30>0)*(underlargef20>0)
  underphat<-underphat+.10*(underlargef>0)*(underlargef40>0)*(underlargef30>0)*
    (underlargef20>0)*(underlargef10>0)
  
  underphat<-underphat-.10*(underlargef<0)*(underlargef60<0)
  underphat<-underphat-.10*(underlargef<0)*(underlargef60<0)*(underlargef70<0)
  underphat<-underphat-.10*(underlargef<0)*(underlargef60<0)*(underlargef70<0)*(underlargef80<0)
  underphat<-underphat-.10*(underlargef<0)*(underlargef60<0)*(underlargef70<0)*
    (underlargef80<0)*(underlargef90<0)
  
  
  linkphat<-plogis(2*largef)
  
  
  
  
  
  error[i]<-sum( abs(phat-realp)^2 /length(largey)  )
  undererror[i]<-sum( abs(underphat-realp)^2 /length(largey)  )
  linkerror[i]<-sum( abs(linkphat-realp)^2 /length(largey)  )
  
  
  
  f50forplot[i]<-median(largef[realp<.60&realp>.50])
  f70forplot[i]<-median(largef[realp<.80&realp>.70])
  f90forplot[i]<-median(largef[realp>.90])
  
  
  
  linkphat[linkphat<.05]<-.05
  linkphat[linkphat>.95]<-.95
  
  
  
  logerror[i]<-1/length(realp)*sum( -1*realp*log(phat)-(1-realp)*log(1-phat)  )
  underlogerror[i]<-1/length(realp)*sum( -1*realp*log(underphat)-(1-realp)*log(1-underphat)  )
  loglinkerror[i]<-1/length(realp)*sum( -1*realp*log(linkphat)-(1-realp)*log(1-linkphat)  )
  
  
  exploss[i]<-1/length(realp)*
    sum( 1*realp*sqrt((1-phat)/phat)+(1-realp)*sqrt(phat/(1-phat)))
  underexploss[i]<-1/length(realp)*
    sum( 1*realp*sqrt((1-underphat)/underphat)+(1-realp)*sqrt(underphat/(1-underphat)))
  linkexploss[i]<-1/length(realp)*
    sum( 1*realp*sqrt((1-linkphat)/linkphat)+(1-realp)*sqrt(linkphat/(1-linkphat)))
  
  
  
  te<-sum(1*f*y<0)/n
  te90<-sum(1*f90*y90<0)/length(y90)
  te10<-sum(1*f10*y10<0)/length(y10)
  misclass[i]<-sum(1*largef*largey<0)/length(largef)
  
  
  print(c(i,te,error[i],linkerror[i],logerror[i],loglinkerror[i],misclass[i],te90,te10))
  
  
  
  i<-i+1
  
  
}




y<-(y+1)/2


model<-logitboost(x, y, largex, (iterations), presel = 0, estimate = 0,verbose = TRUE)
alllblinkphat<-model$prob[,1:(iterations)]

summarize(model,(largey+1)/2, mout=(iterations))

y<-2*y-1


for(i in 1:iterations){
  lblinkphat<-alllblinkphat[,i]
  lblargef<-qlogis(lblinkphat)
  
  lbf50forplot[i]<-median(lblargef[realp<.60&realp>.50])
  lbf70forplot[i]<-median(lblargef[realp<.80&realp>.70])
  lbf90forplot[i]<-median(lblargef[realp>.90])
  
  lblinkerror[i]<-sum( abs(lblinkphat-realp)^2 /length(largey)  )
  
  lblinkphat[lblinkphat<.05]<-.05
  lblinkphat[lblinkphat>.95]<-.95
  
  lbloglinkerror[i]<-1/length(realp)*
    sum( -1*realp*log(lblinkphat)-(1-realp)*log(1-lblinkphat)  )
  
  lblinkexploss[i]<-1/length(realp)*
    sum( 1*realp*sqrt((1-lblinkphat)/lblinkphat)+(1-realp)*sqrt(lblinkphat/(1-lblinkphat)))
  
  lbmisclass[i]<-sum(1*lblargef*largey<0)/length(largef)
  
}




print("ours")
print(error[1:100])
print(error[(stop)])
print(min(error[1:(stop)]))

print("link")
print(linkerror[1:100])
print(linkerror[(stop)])
print(min(linkerror[1:(stop)]))



#set zero iter

proportion<-sum(y==1)/n
zeroprob<-rep(proportion,m)

if (proportion>=.5) {
  zeroclass<-rep(1,m)
}

if (proportion<.5) {
  zeroclass<-rep(-1,m)
}


misclass<-c(sum(1*zeroclass*largey<0)/m,misclass)
lbmisclass<-c(sum(1*zeroclass*largey<0)/m,lbmisclass)


error<-c( sum((realp-zeroprob)^2)/m ,error)      
undererror<-c( sum((realp-zeroprob)^2)/m ,undererror)
linkerror<-c( sum((realp-zeroprob)^2)/m ,linkerror)
lblinkerror<-c( sum((realp-zeroprob)^2)/m ,lblinkerror)


logerror<-c(1/m*sum( -1*realp*log(zeroprob)-(1-realp)*log(1-zeroprob)),logerror)
underlogerror<-c(1/m*sum( -1*realp*log(zeroprob)-(1-realp)*log(1-zeroprob)),underlogerror)
loglinkerror<-c(1/m*sum( -1*realp*log(zeroprob)-(1-realp)*log(1-zeroprob)),loglinkerror)
lbloglinkerror<-c(1/m*sum( -1*realp*log(zeroprob)-(1-realp)*log(1-zeroprob)),lbloglinkerror)

exploss<-c(1/m*
             sum( 1*realp*sqrt((1-zeroprob)/zeroprob)+(1-realp)*sqrt(zeroprob/(1-zeroprob))),exploss)   
underexploss<-c(1/m*
                  sum( 1*realp*sqrt((1-zeroprob)/zeroprob)+(1-realp)*sqrt(zeroprob/(1-zeroprob))),underexploss)  
linkexploss<-c(1/m*
                 sum( 1*realp*sqrt((1-zeroprob)/zeroprob)+(1-realp)*sqrt(zeroprob/(1-zeroprob))),linkexploss)  
lblinkexploss<-c(1/m*
                   sum( 1*realp*sqrt((1-zeroprob)/zeroprob)+(1-realp)*sqrt(zeroprob/(1-zeroprob))),lblinkexploss)  
iter<-c(0,iter)




#######################################################


postscript("1000dsqauredloss.ps")
par(pty="s",cex=1.5,cex.axis=1.8,cex.lab=2.5)
plot(c(iter,iter,iter,iter),c(error,linkerror,lblinkerror,undererror),type="n",
     xlab="Iterations",ylab="Squared Loss")
points(0,error[1],pch=19,cex=3)
lines(iter,error,lwd=3,col="red")
lines(iter,linkerror,lwd=3)
lines(iter,lblinkerror,lwd=3,col="blue")
lines(iter,undererror,lwd=3,col="darkgreen")
dev.off()



postscript("1000dmisclass.ps")
par(pty="s",cex=1.5,cex.axis=1.8,cex.lab=2.5)
plot(c(iter,iter),c(misclass,lbmisclass),
     type="n",
     xlab="Iterations",ylab="Misclassification Error")
points(0,misclass[1],pch=19,cex=3)
lines(iter,misclass,lwd=3)
lines(iter,lbmisclass,lwd=3,col="blue")
dev.off()


postscript("1000dlogloss.ps")
par(pty="s",cex=1.5,cex.axis=1.8,cex.lab=2.5)
plot(c(iter,iter,iter,iter),c(logerror,loglinkerror,lbloglinkerror,underlogerror),
     type="n",
     xlab="Iterations",ylab="Log Loss")
points(0,logerror[1],pch=19,cex=3)
lines(iter,logerror,lwd=3,col="red")
lines(iter,loglinkerror,lwd=3)
lines(iter,lbloglinkerror,lwd=3,col="blue")
lines(iter,underlogerror,lwd=3,col="darkgreen")
dev.off()



postscript("1000dexploss.ps")
par(pty="s",cex=1.5,cex.axis=1.8,cex.lab=2.5)
plot(c(iter,iter,iter,iter),c(exploss,linkexploss,lblinkexploss,underexploss),
     type="n",
     xlab="Iterations",ylab="Exponential Loss")
points(0,exploss[1],pch=19,cex=3)
lines(iter,exploss,lwd=3,col="red")
lines(iter,linkexploss,lwd=3)
lines(iter,lblinkexploss,lwd=3,col="blue")
lines(iter,underexploss,lwd=3,col="darkgreen")
dev.off()

iter<-iter[-c(1)]


postscript("1000df50plot.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,f50forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,f50forplot)
dev.off()


postscript("1000df70plot.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,f70forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,f70forplot)
dev.off()

postscript("1000df90plot.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,f90forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,f90forplot)
dev.off()


postscript("1000df50plotlb.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,lbf50forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,lbf50forplot)
dev.off()



postscript("1000df70plotlb.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,lbf70forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,lbf70forplot)
dev.off()


postscript("1000df90plotlb.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(iter,lbf90forplot,type="n",xlab="Iterations",ylab=expression(F[m]))
lines(iter,lbf90forplot)
dev.off()


postscript("1000df50andf90plot.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(c(iter,iter),c(f90forplot,f50forplot),
     type="n",ylab=expression(F[m]),xlab="Iterations")
lines(iter,f90forplot)
lines(iter,f50forplot)
dev.off()


postscript("1000df50andf90plotlb.ps")
par(pty="s",cex.axis=2,cex.lab=2)
plot(c(iter,iter),c(lbf90forplot,lbf50forplot),
     type="n",ylab=expression(F[m]),xlab="Iterations")
lines(iter,lbf90forplot)
lines(iter,lbf50forplot)
dev.off()



postscript("1000drealphist.ps")
par(pty="s",cex.axis=2,cex.lab=2)
hist.default(realp,main="",xlab ="", col="lightblue",breaks=10)
dev.off()


postscript("1000dlinkphist.ps")
par(pty="s",cex.axis=2,cex.lab=2)
hist.default(linkphat,main="",xlab ="", col="lightblue",breaks=10)
dev.off()


postscript("1000dlinkphistlb.ps")
par(pty="s",cex.axis=2,cex.lab=2)
hist.default(lblinkphat,main="",xlab ="", col="lightblue",breaks=10)
dev.off()

postscript("1000dourphist.ps")
par(pty="s",cex.axis=2,cex.lab=2)
hist.default(phat,main="",xlab ="", col="lightblue",breaks=10)
dev.off()


postscript("1000dourphistunder.ps")
par(pty="s",cex.axis=2,cex.lab=2)
hist.default(underphat,main="",xlab ="", col="lightblue",breaks=10)
dev.off()






########comparison######################
i<-i+1

print("##########")
print("Over Sampling")
print("Misclassification Error")
print(misclass[i])
print("Squared Error")
print(error[i])
print("Log Loss")
print(logerror[i])
print("Exponential Loss")
print(exploss[i])


print("##########")
print("Under Sampling")
print("Misclassification Error")
print(misclass[i])
print("Squared Error")
print(undererror[i])
print("Log Loss")
print(underlogerror[i])
print("Exponential Loss")
print(underexploss[i])


print("##########")
print("AdaBoost Link")
print("Misclassification Error")
print(misclass[i])
print("Squared Error")
print(linkerror[i])
print("Log Loss")
print(loglinkerror[i])
print("Exponential Loss")
print(linkexploss[i])


print("##########")
print("LogitBoost Link")
print("Misclassification Error")
print(lbmisclass[i])
print("Squared Error")
print(lblinkerror[i])
print("Log Loss")
print(lbloglinkerror[i])
print("Exponential Loss")
print(lblinkexploss[i])

print("##########")
print("CART")
treee<-rpart(y~.,data.frame(x),method="class")

strawphat<-predict(treee,data.frame(largex))[,2]
strawclass<-(-1)+2*(strawphat>=.5)

print("Misclassification Error")
print(sum(1*strawclass*largey<0)/m)

print("Squared Error")
print(sum((realp-strawphat)^2)/m)


strawphat[strawphat<.05]<-.05
strawphat[strawphat>.95]<-.95


print("Log Loss")
print(1/m*sum( -1*realp*log(strawphat)-(1-realp)*log(1-strawphat)))

print("Exponential Loss")
print(1/m*sum( 1*realp*sqrt((1-strawphat)/strawphat)+(1-realp)*sqrt(strawphat/(1-strawphat))))





print("Random Forest")
print("##########")

y<-as.factor(y)

rf<-randomForest(x,y)
strawphat<-predict(rf,largex,type="prob")[,2]
strawclass<-(-1)+2*(strawphat>=.5)


y<-(-3)+2*as.numeric(y)


print("Misclassification Error")
print(sum(1*strawclass*largey<0)/m)

print("Squared Error")
print(sum((realp-strawphat)^2)/m)


strawphat[strawphat<.05]<-.05
strawphat[strawphat>.95]<-.95


print("Log Loss")
print(1/m*sum( -1*realp*log(strawphat)-(1-realp)*log(1-strawphat)))

print("Exponential Loss")
print(1/m*sum( 1*realp*sqrt((1-strawphat)/strawphat)+(1-realp)*sqrt(strawphat/(1-strawphat))))



print("K Nearest Neighbor")
print("##########")

for(aaaa in 1:30){
  fit<-knn(x,largex,y,k=aaaa,prob=T)
  ppp<-attributes(fit)$prob
  ppp[as.numeric(fit)==1]<-1-ppp[as.numeric(fit)==1]
  strawphat<-ppp
  strawclass<-(-1)+2*(strawphat>=.5)
  
  print("........")
  print(aaaa)
  
  print("Misclassification Error")
  print(sum(1*strawclass*largey<0)/m)
  
  print("Squared Error")
  print(sum((realp-strawphat)^2)/m)
  
  strawphat[strawphat<.05]<-.05
  strawphat[strawphat>.95]<-.95
  
  print("Log Loss")
  print(1/m*sum( -1*realp*log(strawphat)-(1-realp)*log(1-strawphat)))
  
  print("Exponential Loss")
  print(1/m*sum( 1*realp*sqrt((1-strawphat)/strawphat)+(1-realp)*sqrt(strawphat/(1-strawphat))))
  
  
}

