#' NEAREST NEIGHBORS
#' 
#' @param x an nxm matrix where each row is a reference time series. There are n series of length m.
#' @param obs a 1Xm matrix that represents a series of length m which we want to classify by comparing with the reference series.
#' @param k number of nearest neighbors
#' @param FUN distance measure.
#' @param p some extra parameter of the FUN function
nearest_neighbors = function(x, obs, k, FUN, p=NULL){
  # Check if the number of observations is equal
  if (ncol(x) != length(obs)){stop('Series must have the same length')}
  
  # Calculate the distance, considering p by Minkowski
  dist = apply(x, 1, FUN, obs, p=p)
  
  # Find nearest neighbor
  idx = order(dist)[1:k]
  distances = dist[idx]
  if (length(unique(distances)) != k){
    warning(
      paste("Several variables with equal distance. k was used:", k)
    )
  }
  
  # Create a nearest_neighbors object class
  ret = list(
    neighbor_ind = idx,
    distances = distances,
    obs = obs,
    k = k,
    FUN = match.fun(FUN),
    p = p,
    x = x
  )
  class(ret) <- "nearest_neighbors"
  return(ret)
}


#' KNN Prediction Function
#'
#' @param x a vector of class labels
#' @return the predicted class label
#' @export
knn_prediction = function(x){
  pred = names(which.max(table(x)))
  return(pred)
}


#' KNN prediction for several series
#'
#' @param x_fit a data frame with reference time series (by row) and first column as class label
#' @param x_pred a data frame with time series to classify (by row) and first column as class label
#' @param k number of nearest neighbors
#' @param func distance measure
#' @param weighted_pred whether to use weighted prediction or not
#' @param p extra parameter for the distance measure function
#' @return a vector of predicted class labels
#' @export
knn = function(x_fit, x_pred, k,
               func = "dist", weighted_pred = FALSE, p = NULL){
  # We initialize the predictions
  predictions = character(nrow(x_pred))
  
  # For each observation, we get the prediction
  for (i in seq_len(nrow(x_pred))){
    neighbors = nearest_neighbors(x_fit[,-1],
                                  x_pred[i,-1], k, FUN = func, p = p)
    
    if (weighted_pred){
      pred = knn_prediction(x_fit[neighbors$neighbor_ind, 1], neighbors$distances)
    }else{
      pred = knn_prediction(x_fit[neighbors$neighbor_ind, 1])
    }
    
    # If there is more than 1 predicted class, make a prediction with k + 1
    if(length(pred) > 1){
      pred = knn(x_fit, x_pred[i,], k=k+1,
                 func = func, weighted_pred = weighted_pred, p=p)
    }
    predictions[i] = pred
  }
  
  return(predictions)
}
