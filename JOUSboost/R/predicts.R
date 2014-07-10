predict.JOUSboost = function(object, newdata, ...){
  ##first compute the boosted output -- then compute the probability 
  num_cuts = object$delta - 1
  preds = matrix(NA, nrow = nrow(newdata), ncol = num_cuts) ##matrix of full set of predicted values
  for(i in 1 : num_cuts){
    ##go over iterations
    f = numeric(nrow(newdata))
    
    for(iter in 1:object$num_iter){
      object[[i]][[iter]]$tree_mod$terms = object$rpart_formula
      tree_preds = predict(object[[i]][[iter]]$tree_mod, data.frame(newdata), type = "prob") ##predict
      pred_classes = -1 + 2 * (tree_preds[ ,1] > tree_preds[ ,2]) ##convert to -1/1
      f = f + object[[i]][[iter]]$alpha * pred_classes ##do update based on training
      object[[i]][[iter]]$tree_mod$terms = NULL
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