#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <R.h>
#include <Rinternals.h>

static const double eps = 0.0;

void x_times_beta(int *restrict x, double *restrict z, double *restrict beta, int *nRows, int *nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict result){
  int n = *nRows;
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int i, p, nLevels, allzero, offset = 0;
  int *restrict xOffsetPtr;
  double *restrict zOffsetPtr;
  double factor;
  /* categorical */
  if (pCat > 0){
    factor = sqrt(n);
    for (p=0; p<pCat; p++){
      nLevels = numLevels[catIndices[p]-1];
      /* check if beta for this variable is all zero. If so, move on */
      allzero = 1;
      for (i=0; i<nLevels; i++){
	if (fabs(beta[offset + i]) > eps){
	  allzero = 0;
	  break;
	}
      }
      if (allzero){
	offset += nLevels;
	continue;
      }
      xOffsetPtr = x + (catIndices[p]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i]] / factor;
      }
      offset += nLevels;
    }
  }
  /* continuous */
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      /* check if beta is zero */
      if (fabs(beta[offset]) < eps){
	++offset;
	continue;
      }
      zOffsetPtr = z + (contIndices[p]-1)*n;
      for (i=0; i<n; i++){
	result[i] += zOffsetPtr[i] * beta[offset];
      }
      ++offset;
    }
  }
  /* categorical-categorical */
  if (pCatCat > 0){
    int len;
    int *restrict yOffsetPtr;
    factor = sqrt(n);
    for (p=0; p<pCatCat; p+=2){
      nLevels = numLevels[catcatIndices[p]-1];
      len = nLevels * numLevels[catcatIndices[p+1]-1];
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<len; i++){
	if (fabs(beta[offset + i]) > eps){
	  allzero = 0;
	  break;
	}
      }
      if (allzero){
	offset += len;
	continue;
      }
      xOffsetPtr = x + (catcatIndices[p]-1)*n;
      yOffsetPtr = x + (catcatIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i] + nLevels*yOffsetPtr[i]] / factor;
      }
      offset += len;
    }
  }
  /* continuous-continuous */
  if (pContCont > 0){
    double *restrict wOffsetPtr;
    factor = sqrt(3);
    double *restrict product = malloc(n * sizeof *product);
    double mean, norm;
    for (p=0; p<pContCont; p+=2){
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<3; i++){
	if (fabs(beta[offset + i]) > eps){
	  allzero = 0;
	  break;
	}
      }
      if (allzero){
	offset += 3;
	continue;
      }
      wOffsetPtr = z + (contcontIndices[p]-1)*n;
      zOffsetPtr = z + (contcontIndices[p+1]-1)*n;
      mean = norm = 0.0;
      for (i=0; i<n; i++){
	result[i] += (wOffsetPtr[i]*beta[offset] + zOffsetPtr[i]*beta[offset+1]) / factor;
	product[i] = wOffsetPtr[i] * zOffsetPtr[i];
	mean += product[i];
	norm += product[i]*product[i];
      }
      if (norm > 0){
	mean /= n;
	norm = sqrt(3 * (norm-n*pow(mean, 2)));
	for (i=0; i<n; i++){
	  result[i] += (product[i]-mean) * beta[offset+2] / norm;
	}
      }
      offset += 3;
    }
    free(product);
  }
  /* categorical-continuous */
  if (pCatCont > 0){
    factor = sqrt(2*n);
    double factorZ = sqrt(2);
    for (p=0; p<pCatCont; p+=2){
      nLevels = numLevels[catcontIndices[p]-1];
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<2*nLevels; i++){
	if (fabs(beta[offset + i]) > eps){
	  allzero = 0;
	  break;
	}
      }
      if (allzero){
	offset += 2*nLevels;
	continue;
      }
      xOffsetPtr = x + (catcontIndices[p]-1)*n;
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i]] / factor;
	result[i] += zOffsetPtr[i] * beta[offset + nLevels + xOffsetPtr[i]] / factorZ;
      }
      offset += 2*nLevels;
    }
  }
}

double compute_loglik(const double *restrict y, const double *restrict linear, const double *restrict intercept, const int *restrict nRows, const int *restrict family){
  double result = 0.0, mu = *intercept;
  int i, n = *nRows;
  if (*family == 0){
    for (i=0; i<n; i++){
      result += pow(y[i]-mu-linear[i], 2);
    }
    result /= (2*n);
  }
  else {
    for (i=0; i<n; i++){
      result += -y[i]*(mu+linear[i]) + log(1+exp(mu+linear[i]));
    }
    result /= n;
  }
  return result;
}

void compute_objective(const double *restrict y, const double *restrict res, const double *restrict linear, const double *restrict intercept, const double *restrict beta, const int *restrict nRows, const int *restrict numGroups, const int *restrict groupSizes, const double *restrict lambda, double *restrict objValue, const int *restrict family){
  int i, j, size, n = *nRows, numgroups = *numGroups, offset = 0;
  double loglik = 0.0, penalty = 0.0, mu = *intercept, temp;
  if (*family == 0){
    for (i=0; i<n; i++) loglik += res[i]*res[i];
    loglik /= (2*n);
  }
  else {
    for (i=0; i<n; i++) loglik += -y[i]*(mu+linear[i]) + log(1+exp(mu+linear[i]));
    loglik /= n;
  }
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    temp = 0.0;
    for (j=0; j<size; j++){
      temp += pow(beta[offset+j], 2);
    }
    penalty += sqrt(temp);
    offset += size;
  }
  *objValue = loglik + penalty*(*lambda);
}

void compute_gradient(int *restrict x, double *restrict z, double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict gradient){
  int n = *nRows;
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int i, p, offset = 0;
  int *restrict xOffsetPtr;
  double *restrict zOffsetPtr;
  double factor;
  /* categorical */
  if (pCat > 0){
    factor = sqrt(n);
    for (p=0; p<pCat; p++){
      xOffsetPtr = x + (catIndices[p]-1)*n;
      for (i=0; i<n; i++){
	gradient[offset + xOffsetPtr[i]] += r[i];
      }
      offset += numLevels[catIndices[p]-1];
    }
    for (i=0; i<offset; i++){
      gradient[i] /= factor;
    }
  }
  /* continuous */
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      zOffsetPtr = z + (contIndices[p]-1)*n;
      for (i=0; i<n; i++){
	gradient[offset] += zOffsetPtr[i] * r[i];
      }
      ++offset;
    }
  }
  /* categorical-categorical */
  if (pCatCat > 0){
    factor = sqrt(n);
    int nLevels, start = offset;
    int *restrict yOffsetPtr;
    for (p=0; p<pCatCat; p+=2){
      xOffsetPtr = x + (catcatIndices[p]-1)*n;
      yOffsetPtr = x + (catcatIndices[p+1]-1)*n;
      nLevels = numLevels[catcatIndices[p]-1];
      for (i=0; i<n; i++){
	gradient[offset + xOffsetPtr[i] + nLevels*yOffsetPtr[i]] += r[i];
      }
      offset += nLevels * numLevels[catcatIndices[p+1]-1];
    }
    for (i=start; i<offset; i++){
      gradient[i] /= factor;
    }
  }
  /* continuous-continuous */
  if (pContCont > 0){
    factor = sqrt(3);
    double *restrict wOffsetPtr;
    double *restrict product = malloc(n * sizeof *product);
    double mean, norm;
    for (p=0; p<pContCont; p+=2){
      wOffsetPtr = z + (contcontIndices[p]-1)*n;
      zOffsetPtr = z + (contcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	gradient[offset] += wOffsetPtr[i] * r[i];
	gradient[offset + 1] += zOffsetPtr[i] * r[i];
      }
      gradient[offset] /= factor;
      gradient[offset + 1] /= factor;
      mean = norm = 0.0;
      for (i=0; i<n; i++){
	product[i] = wOffsetPtr[i] * zOffsetPtr[i];
	mean += product[i];
	norm += product[i]*product[i];
      }
      if (norm > 0){
	mean /= n;
	norm = sqrt(3 * (norm-n*pow(mean, 2)));
	for (i=0; i<n; i++){
	  gradient[offset + 2] += (product[i]-mean) * r[i];
	}
	  gradient[offset + 2] /= norm;
      }
      offset += 3;
    }
    free(product);
  }
  /* categorical-continuous */
  if (pCatCont > 0){
    factor = sqrt(2*n);
    double factorZ = sqrt(2);
    int nLevels;
    for (p=0; p<pCatCont; p+=2){
      nLevels = numLevels[catcontIndices[p]-1];
      xOffsetPtr = x + (catcontIndices[p]-1)*n;
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	gradient[offset + xOffsetPtr[i]] += r[i];
	gradient[offset + nLevels + xOffsetPtr[i]] += zOffsetPtr[i] *r[i];
      }
      for (i=offset; i<offset+nLevels; i++){
	gradient[i] /= factor;
      }
      for (i=offset+nLevels; i<offset+2*nLevels; i++){
	gradient[i] /= factorZ;
      }
      offset += 2*nLevels;
    }
  }
  /* normalize by n */
  for (i=0; i<offset; i++){
    gradient[i] /= -n;
  }
}

void compute_group_info(const int *nVars, const int *restrict numLevels, const int *restrict catIndices, const int *restrict contIndices, const int *restrict catcatIndices, const int *restrict contcontIndices, const int *restrict catcontIndices, int *length, int *groupSizes){
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int p, counter = 0, len = 0;
  if (pCat > 0){
    for (p=0; p<pCat; p++){
      groupSizes[counter] = numLevels[catIndices[p]-1];
      len += groupSizes[counter++];
    }
  }
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      groupSizes[counter++] = 1;
      ++len;
    }
  }
  if (pCatCat > 0){
    for (p=0; p<pCatCat; p+=2){
      groupSizes[counter] = numLevels[catcatIndices[p]-1] * numLevels[catcatIndices[p+1]-1];
      len += groupSizes[counter++];
    }
  }
  if (pContCont > 0){
    for (p=0; p<pContCont; p+=2){
      groupSizes[counter++] = 3;
      len += 3;
    }
  }
  if (pCatCont > 0){
    for (p=0; p<pCatCont; p+=2){
      groupSizes[counter] = 2 * numLevels[catcontIndices[p]-1];
      len += groupSizes[counter++];
    }
  }
  *length = len;
}

void compute_update(const double *restrict beta, double *restrict betaUpdated, const double *restrict gradient, const int *restrict groupSizes, const int *restrict numGroups, const double *restrict stepsize, const double *restrict lambda){
  int i, j, size, offset = 0, numgroups = *numGroups;
  double step = *stepsize, factor = step * (*lambda);
  double norm;
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    norm = 0.0;
    for (j=0; j<size; j++){
      betaUpdated[offset+j] = beta[offset+j] - step*gradient[offset+j];
      norm += betaUpdated[offset+j]*betaUpdated[offset+j];
    }
    norm = sqrt(norm);
    for (j=0; j<size; j++){
      betaUpdated[offset+j] *= fmax(0.0, 1.0 - factor/norm);
    }
    offset += size;
  }
}

void optimize_step(int *restrict x, double *restrict z, const double *restrict y, const double *restrict residual, double *restrict linear, int *restrict nRows, int *restrict numGroups, int *restrict groupSizes, int *restrict gradientLength, const double *restrict intercept, double *restrict beta, double *restrict betaUpdated, const double *restrict gradient, double *restrict stepsize, const double *restrict lambda, const double *restrict alpha, int *restrict nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, const int *restrict family){
  int i, n = *nRows, len = *gradientLength;
  double step = *stepsize;
  double loglik, loglikUpdated;
  loglik = compute_loglik(y, linear, intercept, nRows, family);
  double *restrict delta = malloc(len * sizeof *delta);
  double gradientTimesDelta, deltaTimesDelta;
  double factor = *alpha;
  double mu = 0.0;
  while (1){
    gradientTimesDelta = 0.0;
    deltaTimesDelta = 0.0;
    compute_update(beta, betaUpdated, gradient, groupSizes, numGroups, &step, lambda);
    for (i=0; i<len; i++){
      delta[i] = betaUpdated[i] - beta[i];
      gradientTimesDelta += gradient[i] * delta[i];
      deltaTimesDelta += delta[i]*delta[i];
    }
    memset(linear, 0, n * sizeof *linear);
    if (*family == 0){/* gaussian case */
      x_times_beta(x, z, delta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
      loglikUpdated = compute_loglik(residual, linear, &mu, nRows, family); /* residual already had intercept subtracted from it */
    }
    else {/* binomial case */
      x_times_beta(x, z, betaUpdated, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
      loglikUpdated = compute_loglik(y, linear, intercept, nRows, family);
    }
    if (loglikUpdated <= loglik + gradientTimesDelta + deltaTimesDelta/(2*step)) break;
    step *= factor;
  }
  *stepsize = step;
  free(delta);
}

int check_convergence(const double *restrict beta, const double *restrict gradient, const int *restrict groupSizes, const int *restrict numGroups, const double *restrict lambda, const double *restrict tol){
  int i, j, size, allzero, offset = 0, numgroups = *numGroups;
  double norm, factor = *lambda, tolerance = *tol;
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    /* check if beta is zero and compute norm of the gradient */
    allzero = 1;
    for (j=0; j<size; j++){
      if (fabs(beta[offset+j]) > eps){
	allzero = 0;
	break;
      }
    }
    norm = 0.0;
    for (j=0; j<size; j++){
      norm += gradient[offset+j]*gradient[offset+j];
    }
    norm = sqrt(norm);
    /* if nonzero, check for kkt condition */
    if (!allzero && fabs(norm-factor)/factor>tolerance) return 0;
    /* if zero, make sure gradient norm is < lambda */
    if (allzero && norm>factor) return 0;
    offset += size;
  }
  return 1;
} 

double update_theta(const double *restrict beta, const double *restrict intermediate, const double *restrict intermediateOld, const int gradientLength, const double theta){
  int i;
  double value = 0.0;
  for (i=0; i<gradientLength; i++){
    value += (beta[i]-intermediate[i]) * (intermediate[i]-intermediateOld[i]);
  }
  return value > 0.0 ? 1.0 : theta;
}

void update_intercept(const double *restrict y, const int *restrict nRows, const double *restrict linear, double *restrict intercept, double *restrict residual, const int *restrict family){
  int i, n = *nRows;
  double residualMean = 0.0, mu = *intercept;
  if (*family == 0){
    for (i=0; i<n; i++){
      residual[i] = y[i] - mu - linear[i];
      residualMean += residual[i];
    }
    residualMean /= n;
    *intercept += residualMean;
    for (i=0; i<n; i++){
      residual[i] -= residualMean;
    }
  }
  else {
    const double xmax = -log(DBL_EPSILON);
    const double xmin = log(DBL_MIN);
    double *restrict temp = malloc(n * sizeof *temp);
    double *restrict exponent = malloc(n * sizeof *exponent);
    double f = 0.0, fPrime, sumY = 0.0, expMu;
    double sum;
    expMu = exp(-mu);
    for (i=0; i<n; i++){
      exponent[i] = exp(-linear[i]);
      temp[i] = expMu * exponent[i];
      sumY += y[i];
      sum = mu + linear[i];
      f += y[i] - (sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i])));
    }
    int iter = 0;
    while (iter<1000 && fabs(f)>1e-2){
      fPrime = 0.0;
      for (i=0; i<n; i++){
	sum = mu + linear[i];
	fPrime -= (sum > xmax || sum < xmin) ? 0.0 : temp[i]/pow(1+temp[i], 2);
      }
      mu -= f/fPrime;
      expMu = exp(-mu);
      f = sumY;
      for (i=0; i<n; i++){
	temp[i] = expMu * exponent[i];
	sum = mu + linear[i];
	f -= sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i]));
      }
      ++iter;
    }
    *intercept = mu;
    for (i=0; i<n; i++){
      sum = mu + linear[i];
      residual[i] = y[i] - (sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i])));
    }
    free(temp);
    free(exponent);
  }
}      

double compute_stepsize(const double *restrict gradient, const double *restrict gradientOld, const double *restrict beta, const double *restrict betaOld, const int gradientLength){
  int i;
  double normBeta = 0.0, normGradient = 0.0;
  for (i=0; i<gradientLength; i++){
    normBeta += (beta[i]-betaOld[i])*(beta[i]-betaOld[i]);
    normGradient += (gradient[i]-gradientOld[i])*(gradient[i]-gradientOld[i]);
  }
  return sqrt(normBeta/normGradient);
}

void gl_solver(int *restrict x, double *restrict z, double *restrict y, int *restrict nRows, double *restrict intercept, double *restrict beta, double *restrict residual, double *restrict linear, int *restrict numLevels, int *restrict nVars, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict lambda, double *restrict tol, double *restrict alpha, int *restrict maxIter, int *restrict convergedFlag, double *restrict objValue, double *restrict steps, int *restrict family){
  /* initialize required variables */
  int i, gradientLength, iter, converged, n = *nRows;
  int numGroups = nVars[0] + nVars[1] + nVars[2] + nVars[3] + nVars[4];
  int *restrict groupSizes = malloc(numGroups * sizeof *groupSizes);
  double theta, thetaOld, momentum, stepsize;
  /* get grouping information and initialize gradient */
  compute_group_info(nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, &gradientLength, groupSizes);
  double *restrict gradient = malloc(gradientLength * sizeof *gradient);
  double *restrict intermediate = malloc(gradientLength * sizeof *intermediate);
  double *restrict intermediateOld = malloc(gradientLength * sizeof *intermediateOld);
  memcpy(intermediate, beta, gradientLength * sizeof *beta);
  double *restrict betaOld = malloc(gradientLength * sizeof *betaOld);
  double *restrict gradientOld = malloc(gradientLength * sizeof *gradientOld);
  /* compute residual from initialized beta */
  x_times_beta(x, z, beta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
  update_intercept(y, nRows, linear, intercept, residual, family);
  /* start accelerated FISTA */
  iter = 0;
  theta = 1.0;
  *convergedFlag = 0;
  while (iter < *maxIter) {
    /* compute gradient */
    memcpy(gradientOld, gradient, gradientLength * sizeof *gradient);
    memset(gradient, 0, gradientLength * sizeof *gradient);
    compute_gradient(x, z, residual, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, gradient);
    /* check convergence */
    converged = check_convergence(beta, gradient, groupSizes, &numGroups, lambda, tol);
    if (converged){
      *convergedFlag = 1;
      break;
    }
    /* compute intermediate update and stepsize */
    memcpy(intermediateOld, intermediate, gradientLength * sizeof *intermediate);
    stepsize = iter > 0 ? compute_stepsize(gradient, gradientOld, beta, betaOld, gradientLength) : 1.0;
    optimize_step(x, z, y, residual, linear, nRows, &numGroups, groupSizes, &gradientLength, intercept, beta, intermediate, gradient, &stepsize, lambda, alpha, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, family);
    /* check if restart required  and update theta*/
    thetaOld = update_theta(beta, intermediate, intermediateOld, gradientLength, theta);
    /* update momentum */
    theta = (1+sqrt(1+4*pow(thetaOld, 2))) / 2;
    momentum = (thetaOld-1) / theta;
    /* update beta */
    memcpy(betaOld, beta, gradientLength * sizeof *beta);
    for (i=0; i<gradientLength; i++){
      beta[i] = intermediate[i] + momentum*(intermediate[i]-intermediateOld[i]);
    }
    /* update residual and mu */
    memset(linear, 0, n * sizeof *linear);
    x_times_beta(x, z, beta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
    update_intercept(y, nRows, linear, intercept, residual, family);
    /* update iteration count */
    //compute_objective(y, residual, linear, intercept, beta, nRows, &numGroups, groupSizes, lambda, objValue+iter, family);
    steps[iter] = stepsize;
    ++iter;
  }
  compute_objective(y, residual, linear, intercept, beta, nRows, &numGroups, groupSizes, lambda, objValue, family);
  free(groupSizes);
  free(gradient);
  free(intermediate);
  free(intermediateOld);
  free(betaOld);
  free(gradientOld);
}

SEXP R_gl_solver(SEXP R_x, SEXP R_z, SEXP R_y, SEXP R_nRows, SEXP R_intercept, SEXP R_beta, SEXP R_residual, SEXP R_linear, SEXP R_numLevels, SEXP R_nVars, SEXP R_catIndices, SEXP R_contIndices, SEXP R_catcatIndices, SEXP R_contcontIndices, SEXP R_catcontIndices, SEXP R_lambda, SEXP R_tol, SEXP R_alpha, SEXP R_maxIter, SEXP R_convergedFlag, SEXP R_objValue, SEXP R_steps, SEXP R_family){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_z = coerceVector(R_z, REALSXP));
  PROTECT(R_y = coerceVector(R_y, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_intercept = coerceVector(R_intercept, REALSXP));
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  PROTECT(R_residual = coerceVector(R_residual, REALSXP));
  PROTECT(R_linear = coerceVector(R_linear, REALSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_catIndices = coerceVector(R_catIndices, INTSXP));
  PROTECT(R_contIndices = coerceVector(R_contIndices, INTSXP));
  PROTECT(R_catcatIndices = coerceVector(R_catcatIndices, INTSXP));
  PROTECT(R_contcontIndices = coerceVector(R_contcontIndices, INTSXP));
  PROTECT(R_catcontIndices = coerceVector(R_catcontIndices, INTSXP));
  PROTECT(R_lambda = coerceVector(R_lambda, REALSXP));
  PROTECT(R_tol = coerceVector(R_tol, REALSXP));
  PROTECT(R_alpha = coerceVector(R_alpha, REALSXP));
  PROTECT(R_maxIter = coerceVector(R_maxIter, INTSXP));
  PROTECT(R_convergedFlag = coerceVector(R_convergedFlag, INTSXP));
  PROTECT(R_objValue = coerceVector(R_objValue, REALSXP));
  PROTECT(R_steps = coerceVector(R_steps, REALSXP));
  PROTECT(R_family = coerceVector(R_family, INTSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict z = REAL(R_z);
  double *restrict y = REAL(R_y);
  int *restrict nRows = INTEGER(R_nRows);
  double *restrict intercept = REAL(R_intercept);
  double *restrict beta = REAL(R_beta);
  double *restrict residual = REAL(R_residual);
  double *restrict linear = REAL(R_linear);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict catIndices = INTEGER(R_catIndices);
  int *restrict contIndices = INTEGER(R_contIndices);
  int *restrict catcatIndices = INTEGER(R_catcatIndices);
  int *restrict contcontIndices = INTEGER(R_contcontIndices);
  int *restrict catcontIndices = INTEGER(R_catcontIndices);
  double *restrict lambda = REAL(R_lambda);
  double *restrict tol = REAL(R_tol);
  double *restrict alpha = REAL(R_alpha);
  int *restrict maxIter = INTEGER(R_maxIter);
  int *restrict convergedFlag = INTEGER(R_convergedFlag);
  double *restrict objValue = REAL(R_objValue);
  double *restrict steps = REAL(R_steps);
  int *restrict family = INTEGER(R_family);
  gl_solver(x, z, y, nRows, intercept, beta, residual, linear, numLevels, nVars, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, lambda, tol, alpha, maxIter, convergedFlag, objValue, steps, family);
  UNPROTECT(23);
  SEXP result = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(result, 0, R_intercept);
  SET_VECTOR_ELT(result, 1, R_beta);
  SET_VECTOR_ELT(result, 2, R_residual);
  SET_VECTOR_ELT(result, 3, R_objValue);
  const char *names[4] = {"mu", "coefficients", "res", "objValue"};
  SEXP sNames = PROTECT(allocVector(STRSXP, 4));
  int i;
  for (i=0; i<4; i++) SET_STRING_ELT(sNames, i, mkChar(names[i]));
  setAttrib(result, R_NamesSymbol, sNames);
  UNPROTECT(2);
  return result;
}
