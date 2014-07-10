
JOUSboost = function(X, y, 
                     type = "under",
                     num_iter = 100, 
                     delta = 10, 
                     nu = 1, 
                     tree_depth = 3, 
                     verbose = TRUE){
  ##check if y is +-1 
  if(!all(y %in% c(-1,1))) stop("All y must be either -1 or +1")
  ##Need to worry about existence of median adaboost... 
  ##for now don't allow grids that don't have it
  if((50 %% (100 /delta)) != 0) stop("There will be no median adaboost for prediction. To be Fixed Later.")
  
  ##basic info
  n = length(y)
  p = ncol(X)
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
    jitter_pos = sapply(1 : p, function(s) sqrt(var(X[,s]))*nu*runif((num_cuts - 1) * num_pos, -1, 1)) ##only need n - of repeats to be jittered
    jitter_neg = sapply(1 : p, function(s) sqrt(var(X[,s]))*nu*runif((num_cuts - 1) * num_neg, -1, 1)) ##only need n - of repeats to be jittered
    
    for(cut in 1 : num_cuts){
      ##start with 1 replicate -- 10% is 9 -1's and 1 1's
      if(cut == median_cut){
          x_jit = X
          y_temp = y
      }else{
        pos_replicates = cut - 1
        neg_replicates = num_cuts - cut
        
        x_temp = X
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
  boost_list = list()
  rpart_formula = NA  
  ##Run the grid
  for(i in 1 : num_cuts){
   if(verbose) cat(paste("Cut Point:", i, "Iter: "))
    boost_list[[i]] = list()

    y = data_sets[[i]]$y
    x = data_sets[[i]]$x
       
    boost_list[[i]] = adaboost_v2(x = x, y = y, num_iter = num_iter, 
                                  tree_depth = tree_depth, verbose = verbose)
   
   ##Kill all the formulas because they waste a lot of memory. 
   if(i == 1) rpart_formula = boost_list[[i]][[1]]$tree_mod$terms
   for(j in 1 : num_iter){
     boost_list[[i]][[j]]$tree_mod$terms = NULL
   }

    if(verbose) cat("\n")
  }
  class(boost_list) = "JOUSboost"
  boost_list[["num_iter"]] = num_iter
  boost_list[["type"]] = "undersampled"
  boost_list[["delta"]] = delta
  boost_list[["nu"]] = nu
  boost_list[["tree_depth"]] = tree_depth
  boost_list[["rpart_formula"]] = rpart_formula
  boost_list[["verbose"]] = verbose
  boost_list 
}


JOUSboost2 = function(X, y, 
                     type = "under",
                     num_iter = 100, 
                     delta = 10, 
                     nu = 1,
                     jitter_factor = 1,
                     bag_frac = 1,
                     tree_depth = 3, 
                     verbose = TRUE){
  ##check if y is +-1 
  if(!all(y %in% c(-1,1))) stop("All y must be either -1 or +1")
  ##Need to worry about existence of median adaboost... 
  ##for now don't allow grids that don't have it
  if((50 %% (100 /delta)) != 0) stop("There will be no median adaboost for prediction. To be Fixed Later.")
  
  ##basic info
  n = length(y)
  p = ncol(X)
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
    jitter_pos = sapply(1 : p, function(s) sqrt(var(X[,s]))*nu*runif((num_cuts - 1) * num_pos, -1, 1)) ##only need n - of repeats to be jittered
    jitter_neg = sapply(1 : p, function(s) sqrt(var(X[,s]))*nu*runif((num_cuts - 1) * num_neg, -1, 1)) ##only need n - of repeats to be jittered
    
    for(cut in 1 : num_cuts){
      ##start with 1 replicate -- 10% is 9 -1's and 1 1's
      if(cut == median_cut){
        x_jit = X
        y_temp = y
      }else{
        pos_replicates = cut - 1
        neg_replicates = num_cuts - cut
        
        x_temp = X
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
  boost_list = list()
  boost_list[["adaboost_models"]] = list()
  rpart_formula = NA  
  control = rpart.control(minsplit = 0, minbucket = 1, cp = -1, 
                          maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, 
                          xval = 0, maxdepth = tree_depth)
  
  ##Run the grid
  for(i in 1 : num_cuts){
    if(verbose) cat(paste("Cut Point:", i, "Iter: "))
    boost_list[[i]] = list()
    
    y = data_sets[[i]]$y
    x = data_sets[[i]]$x
    
    ##need to add rpart control stuff here... 
    boost_list[["adaboost_models"]][[i]] = ada(x = x, y = y, bag.frac = bag_frac,iter = num_iter,
                                               nu = nu, verbose = T, control = control)
    
    ##Kill all the formulas because they waste a lot of memory. 
    if(i == 1) rpart_formula = boost_list[["adaboost_models"]][[i]]$model$trees[[1]]$terms
    for(j in 1 : num_iter){
      boost_list[["adaboost_models"]][[i]]$model$trees[[i]]$terms = NULL
    }
    
    if(verbose) cat("\n")
  }
  class(boost_list) = "JOUSboost"
  boost_list[["num_iter"]] = num_iter
  boost_list[["type"]] = "undersampled"
  boost_list[["delta"]] = delta
  boost_list[["nu"]] = nu
  boost_list[["tree_depth"]] = tree_depth
  boost_list[["rpart_formula"]] = rpart_formula
  boost_list[["verbose"]] = verbose
  boost_list 
}
