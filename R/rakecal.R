
#' Raking calibration of design weights
#'
#' Design weights (also known as sampling weights, survey weights, or inverse probability weights) are calibrated using the raking method. The source code is similar to calib() of 'sampling' package.
#' @param A a matrix of auxiliary statistics for weight calibration; this should have the size of phase 2 (subsample of a cohort)
#' @param w0 a vector of design weights; this should have the size of phase 2 (subsample of a cohort)
#' @param total a vector of cohort totals of A
#' @param eta logical; if TRUE, a vector of estimates of weight calibration parameters is attributed
#' @param max_iter maximum number of iterations
#' @param EPS epsilon of machine controlling the precision
#' @keywords ahazcal
#' @return a vector of calibrated weights of phase 2 samples with calibration parameter estimates (eta)
#' @export 
#' @seealso \code{\link{sample_data}}, \code{\link{aux.ahaz}}
#' @examples #### STEP1: Predict Missing Covariates ####
#' fit_Z2 = glm(Z2 ~ X1 + X2 + Z1 + U, family = binomial(link = "logit"), 
#'              data = sample_data, subset = (ind.ph2 == 1), weights = 1/incl.prob)
#' 
#' sample_data.tilde = sample_data
#' sample_data.tilde$Z2 = predict(fit_Z2, type = "response", newdata = sample_data)
#' 
#' 
#' #### STEP2: Obtain Auxiliary Statistics (Influence Functions) ####
#' aux = with(sample_data.tilde, aux.ahaz(outcome = ind.fail, etime = eventime, tau = 8, X = cbind(X1, X2), Z = cbind(Z1, Z2)))
#' 
#' aux1 = cbind(1, aux[,c("gamma1","gamma2")]) # Auxiliary 1: Influence functions of gamma estimates 
#' aux2 = cbind(1, aux) # Auxiliary 2: Influence functions of Beta and gamma estimates
#' 
#' 
#' #### STEP3: Calibrate Weights ####
#' cal1 = with(sample_data.tilde, rakecal(A = aux1[ind.ph2,], w0 = 1/incl.prob[ind.ph2], 
#'                                        total = colSums(aux1), eta = T, EPS = 1e-12))
#' cal2 = with(sample_data.tilde, rakecal(A = aux2[ind.ph2,], w0 = 1/incl.prob[ind.ph2], 
#'                                        total = colSums(aux2), eta = T, EPS = 1e-12))
#' 
rakecal = function(A, w0, total, eta = T, max_iter = 500, EPS = 1e-11){
  
  lambda = as.matrix(rep(0, ncol(A)))
  w1 = as.vector(w0 * exp(A %*% lambda))
  
  for (l in 1:max_iter) {
    phi = crossprod(A, w1) - total
    phiprim = crossprod(A * w1, A)
    lambda1 = lambda - ginv(phiprim, tol = .Machine$double.eps) %*% phi
    w1 = as.vector(w0 * exp(A %*% lambda1))
    
    if (any(is.na(w1)) | any(is.infinite(w1))) {
      warning("No convergence")
      g = NULL
      break
    }
    tr = crossprod(A, w1)
    
    if (max(abs(tr - total)) < EPS) break else lambda = lambda1
  }
  
  if (l == max_iter){
    # warning("No convergence")
    g = NULL
  }
  else g = w1
  
  if(eta == T) attributes(g) <- list(eta = lambda1)
  
  return(g)
} 
