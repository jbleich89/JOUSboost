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



jb = jousboost(x, y, type = "over", num_iter = 250)

jousboost = function(X, y, type = "under", num_iter = 100, delta = 10, nu = 1, tree_depth = 3){
  ##check if y is +-1 
  if(!all(y %in% c(-1,1))) stop("All y must be either -1 or +1")
  ##Need to worry about existence of median adaboost... 
  ##for now don't allow grids that don't have it
  if((50 %% (100 /delta)) != 0) stop("There will be no median adaboost for prediction. To be Fixed Later.")
  
  ##basic info
  n = length(y)
  p = ncol(x)
  ix_pos = y == 1
  ix_neg = y == -1
  num_pos = sum(ix_pos)
  num_neg = length(y) - num_pos
  
  xpos = X[ix_pos, ]
  xneg = X[ix_neg, ]
  ypos = y[ix_pos]
  yneg = y[ix_neg]
  
  num_cuts = delta - 1
  cut_pos = round(num_pos / num_cuts, 0)
  cut_neg = round(num_neg / num_cuts, 0)
  
  ##set up info for over/undersampling
  data_sets = list()
  median_cut = median(1 : num_cuts)
  
  if(type == "over"){
    ##set up jittering
    jitter_pos = sapply(1 : p, function(s) sqrt(var(x[,s]))*nu*runif((num_cuts - 1) * num_pos, -1, 1)) ##only need n - of repeats to be jittered
    jitter_neg = sapply(1 : p, function(s) sqrt(var(x[,s]))*nu*runif((num_cuts - 1) * num_neg, -1, 1)) ##only need n - of repeats to be jittered
    
    for(cut in 1 : num_cuts){
      ##start with 1 replicate -- 10% is 9 -1's and 1 1's
      pos_replicates = cut - 1
      neg_replicates = num_cuts - cut
      
      x_temp = x
      y_temp = y
      if(neg_replicates > 0){
        ##augment some negative replicates
        for(j in 1 : neg_replicates){
          x_temp = rbind(x_temp, xneg)
          y_temp = c(y_temp, yneg)
        }
      }
      if(pos_replicates > 0){
        ##augment some positive replicates
        for(j in 1 : pos_replicates){
          x_temp = rbind(x_temp, xpos)
          y_temp = c(y_temp, ypos)
        }
      }
      ##add the noise
      if(cut == 1){
        x_jit = sapply(1 : p, function(s) x_temp[,s] + c(rep(0,n), jitter_neg[,s]))
      }else if(cut == num_cuts){
        x_jit = sapply(1 : p, function(s) x_temp[,s] + c(rep(0,n), jitter_pos[,s]))
      }else{
        x_jit = sapply(1 : p, function(s) x_temp[,s] + c(rep(0,n), 
        jitter_neg[1 : (neg_replicates * num_neg) , s], jitter_pos[1 : (pos_replicates * num_pos) , s]))  
      }
      
      data_sets[[cut]] = list(x = x_jit, y = y_temp)
    } 
  }else{ ##undersample
    for(cut in 1 : num_cuts){ ## 1 corresponds to all of the negative data so its the 90th quantile cut
      if(cut == median_cut){
        x_temp = X
        y_temp = y
      }else{
        pos_iter = cut * cut_pos
        neg_iter = (delta - cut) * cut_neg
        if(cut == 1){
          neg_iter = num_neg
        }else if(cut == num_cuts){
          pos_iter = num_pos
        } 
        x_temp = rbind(xpos[1 : pos_iter, ], xneg[1 : neg_iter, ])
        y_temp = c(ypos[1 : pos_iter], yneg[1: neg_iter])
      }
      data_sets[[cut]] = list(x  = x_temp, y = y_temp)
    }
  }
  
  #do the boosting on each data set here
  out_list = list()
    
  ##Run the grid
  for(i in 1 : num_cuts){
    cat(paste("Cut Point:", i, "Iter: "))
    out_list[[i]] = list()

    y = data_sets[[i]]$y
    x = data_sets[[i]]$x
       
    #out_list[[i]] = adaboost_jb(x = x, y = y, num_iter = num_iter)
    out_list[[i]] = adaboost_dm(x = x, y = y, num_iter = num_iter)

    cat("\n")
  }
  class(out_list) = "JOUSboost"
  out_list[["num_iter"]] = num_iter
  out_list[["type"]] = "undersampled"
  out_list[["delta"]] = delta
  out_list[["nu"]] = nu
  out_list[["tree_depth"]] = tree_depth
  
  out_list 
}

adaboost_jb = function(x, y, num_iter){
  f = numeric(length(y))
  ada_list = list()
  weights = rep(1, length(y)) / length(y)
  for(iter in 1 : num_iter){
    
    ##call rpart
    tree_mod = rpart(y ~ ., data = data.frame(x), weights = weights, method = "class", y = FALSE,
                     control = rpart.control(minsplit = 0, minbucket = 1, cp = -1, 
                                             maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, 
                                             xval = 0, maxdepth = tree_depth))
    
    tree_preds = predict(tree_mod, data.frame(x), type = "prob")
    pred_classes = -1 + 2 * (tree_preds[ ,1] > tree_preds[ ,2])
    
    ##next set of lines is the boosting update
    err = sum(weights * (pred_classes != y))
    if( (1-err) == 1 | err == 1 ){
      err = (1 - err) * 0.0001 + err * .9999
    }
    
    alpha = .5 * log((1 - err)/ err)
    f = f + alpha * pred_classes
    weights = weights * exp(- alpha *pred_classes * y)
    weights = weights/sum(weights)
    
    
    ada_list[[iter]] = list(tree_mod = tree_mod, f = f, alpha = alpha)  ##storage for forecasting later
    if(iter %% 10 == 0) cat(paste(iter,""))
  }
  ada_list
}

object = jb
newdata = xtest
phats = predict.JOUSboost(object, newdata)
sum(abs(phats - ptest)^2/2500)
head(ptest)
head(phats)
head(preds)
sum(phats)

predict.JOUSboost = function(object, newdata){
  ##first compute the boosted output -- then compute the probability 
  num_cuts = object$delta - 1
  preds = matrix(NA, nrow = nrow(newdata), ncol = num_cuts) ##matrix of full set of predicted values
  for(i in 1 : num_cuts){
    ##go over iterations
    f = numeric(nrow(newdata))
    
    for(iter in 1:object$num_iter){
      tree_preds = predict(object[[i]][[iter]]$tree_mod, data.frame(newdata), type = "prob") ##predict
      pred_classes = -1 + 2 * (tree_preds[ ,1] > tree_preds[ ,2]) ##convert to -1/1
      f = f + object[[i]][[iter]]$alpha * pred_classes ##do update based on training
    }
    
    preds[ ,i] = f
  }
  ##post-process results 
  phats = apply(preds, 1, get_gridded_prob, delta = object$delta)
  phats 
}


##private function for dealing with grid 
get_gridded_prob = function(obs_vec, delta){
  quants = seq(1 - 1/delta, 1/delta, by = -1/delta)
  median_col = median(1 : (delta - 1))
  if(obs_vec[median_col] < 0){
    ix_to_check = (median_col + 1) : length(obs_vec)
    temp = which(obs_vec[ix_to_check] > 0)[1] ##first occurence less than 0 
    if(length(temp == 1) & !is.na(temp)){
     phat =  .5 - abs(median_col - ix_to_check[temp]) * delta/100 + 1/(2 * delta)
    }else{
      phat = 1/(2 *delta)
    }
  }else{
    ix_to_check = 1 : (median_col - 1)
    temp = which(obs_vec[ix_to_check] < 0)[ length(which(obs_vec[ix_to_check] < 0))] ##last occurence greater than 0 
    if(length(temp) == 1){
     phat = .5 + abs(median_col - temp) * delta/100 - 1/(2 * delta)
    }else{
    phat =  1 - 1/(2 * delta)
    }
   
  }
  phat
}


####

adaboost_dm = function(x, y, num_iter){
  f = numeric(length(y))
  ada_list = list()
  for(iter in 1 : num_iter){
    weights = exp(-y * f)
    weights = weights / max(weights)
    weights = weights * 1.000001 ##stability?
    
    ##call rpart
    tree_mod = rpart(y ~ ., data = data.frame(x), weights = weights, method = "class", y = FALSE,
                     control = rpart.control(minsplit = 0, minbucket = 1, cp = -1, 
                                             maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, 
                                             xval = 0, maxdepth = tree_depth))
    
    tree_preds = predict(tree_mod, data.frame(x), type = "prob")
    pred_classes = -1 + 2 * (tree_preds[ ,1] > tree_preds[ ,2])
    
    ##next set of lines is the boosting update
    temp = exp(-y * f)
    e = sum(temp * (1 - (y == pred_classes)))
    t = sum(temp)
    alpha = .5 * log( (t - e) / e)
    f = f + alpha * pred_classes
    
    ada_list[[iter]] = list(tree_mod = tree_mod, f = f, alpha = alpha)  ##storage for forecasting later
    if(iter %% 10 == 0) cat(paste(iter,""))
  }
  ada_list
}

