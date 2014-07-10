##Coded Adaboost Implementations which can be substituted into the JOUSboost function
##For internal use only for now.
adaboost_v1 = function(x, y, num_iter, tree_depth, verbose){
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
    
    ada_list[[iter]] = list(tree_mod = tree_mod, alpha = alpha)  ##storage for forecasting later
    if(verbose & iter %% 10 == 0) cat(paste(iter,""))
  }
  ada_list
}


adaboost_v2 = function(x, y, num_iter, tree_depth, verbose){
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
    
    
    ada_list[[iter]] = list(tree_mod = tree_mod, alpha = alpha)  ##storage for forecasting later
    if(verbose & iter %% 10 == 0) cat(paste(iter,""))
  }
  ada_list
}