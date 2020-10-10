group_lasso = function(X, Z, Y, activeSet, betahat, numLevels, lambda, family, tol, maxIter, verbose){

  #get size of each group, number of groups in each category, and total number of groups
  groupSizes = get_group_sizes(activeSet, numLevels)
  numGroups = sapply(activeSet, function(x) if (is.null(x)) 0 else nrow(x))
  totalGroups = sum(numGroups)

  #if active set is empty, just return estimate of the intercept
  if (totalGroups == 0){
    res = Y - mean(Y)
    betahat = ifelse(family=="gaussian",
                     mean(Y),
                     -log(1/mean(Y)-1))
    objValue = ifelse(family=="gaussian",
                      sum(res^2)/(2*length(Y)),
                      -mean(Y)*betahat[[1]]+log(1/(1-mean(Y))))
    return(list(betahat=betahat, activeSet=activeSet, res=res, objValue=objValue))
  }

  n = length(Y)
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)

  #fit and get new betahat, res, objValue
  fit = .Call("R_gl_solver", X, Z, Y, n, betahat[1], betahat[-1], numLevels, numGroups, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, lambda, tol, 0.1, maxIter, ifelse(family=="gaussian", 0, 1), ifelse(verbose, 1, 0))
  res = fit$res
  objValue = fit$objValue

  #get the nonzero parts of betahat and update activeSet
  idx = .Call("R_retrieve_beta", fit$coefficients, groupSizes, totalGroups, integer(totalGroups), integer(length(fit$coefficients)))
  beta = c(fit$mu, fit$coefficients[idx$betaIdx != 0])
  range = c(0, cumsum(numGroups))
  activeSet = lapply(1:5, function(i){
    if (numGroups[i] > 0){
      index = which(idx$idx[(range[i]+1):range[i+1]] != 0)
      if (length(index) > 0) matrix(activeSet[[i]][index, ], nrow=length(index))
      else NULL
    }
    else NULL
  })
  names(activeSet) = c("cat", "cont", "catcat", "contcont", "catcont")

  #output
  list(betahat=beta, activeSet=activeSet, res=res, objValue=objValue)
}




