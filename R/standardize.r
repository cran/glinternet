standardize = function(x){
  if (length(unique(x)) == 1) return(rep(0, length(x)))
  result = x - mean(x)
  return(result / rep(mynorm(result), length(result)))
}
