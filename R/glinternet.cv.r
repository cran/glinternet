glinternet.cv = function(X, Y, numLevels, nFolds=10, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, interactionCandidates=NULL, screenLimit=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1){

                                        #get call and family
  thisCall = match.call()
  family = match.arg(family)

                                        #make sure inputs are valid
  n = length(Y)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")

    fullfitted = glinternet(X, Y, numLevels, lambda, nLambda, lambdaMinRatio, interactionCandidates,screenLimit, family=family, tol=tol, maxIter=maxIter, verbose=verbose, numCores=numCores)
    if(verbose)cat("\n Done fit on all data\n")

    lambda=fullfitted$lambda
    nlambda=length(lambda)
                                        #create the folds
  folds = sample(rep(1:nFolds, ceiling(n/nFolds)), n, replace=FALSE)

                                        #perform cv
  compute_loss = function(y, yhat, family){
    if (family == "gaussian") return (sum((y-yhat)^2)/(2*length(y)))
    yhat = sapply(yhat, function(x) min(max(1e-15, x), 1-1e-15))
    -(t(y)%*%log(yhat) + t(1-y)%*%log(1-yhat))/length(y)
  }
  loss = matrix(0, nFolds, nLambda)

    X=as.matrix(X)
    for (fold in 1:nFolds){
        testIndex= folds == fold
        trainIndex = !testIndex
        fitted = glinternet(X[trainIndex,,drop=FALSE], Y[trainIndex], numLevels, lambda, nLambda, lambdaMinRatio, interactionCandidates,screenLimit, numToFind=NULL, family, tol, maxIter, verbose, numCores)
        YtestHat = predict(fitted, X[testIndex,,drop=FALSE], "response")
        loss[fold, ] = apply(YtestHat, 2, function(yhat) compute_loss(Y[testIndex], yhat, family))
        if(verbose)cat("\n Done fold",fold,"\n")
  }

                                        #compute cv errors and get minimum
  cv = apply(loss, 2, mean)
  cvStd = apply(loss, 2, sd)
  bestIndex1Std = which(cv <= min(cv)+cvStd[which.min(cv)])
  bestIndex = which.min(cv)
  if (length(bestIndex1Std) == nLambda) lambdaHat1Std = lambdaHat = lambda[bestIndex]
  else {
    lambdaHat1Std = lambda[bestIndex1Std[1]]
    lambdaHat = lambda[bestIndex]
  }

### Return fit on full dataset with chosen lambda
 output = list(call=thisCall, glinternetFit=fullfitted, fitted=fullfitted$fitted[, bestIndex], activeSet=fullfitted$activeSet[bestIndex], betahat=fullfitted$betahat[bestIndex], lambda=lambda, lambdaHat=lambdaHat, lambdaHat1Std=lambdaHat1Std, cvErr=cv, cvErrStd=cvStd, family=family, numLevels=numLevels, nFolds=nFolds)
  class(output) = "glinternet.cv"
  return (output)
}
