
#' Estimating A Covariate-specific Pure Risk and Variance Estimate
#'
#' A pure risk is estimated from additive hazards model fit for a given vector of covariates. Its variance estimate (standard error) is also computed based on influence functions of model coefficients. 
#' @param fit an object of ahaz.wcal()
#' @param xzval a n x (p+q) matrix of X- and Z-levels for covariate-specific pure risk estimation
#' @param tau pure risk projection time
#' @keywords ahazcal
#' @return Pure risk estimate and standard error for given covariate level
#' @export 
#' @seealso \code{\link{sample_data}}, \code{\link{aux.ahaz}}, \code{\link{rakecal}}, \code{\link{ahaz.wcal}}
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
#' 
#' #### STEP4: Fit Additive Hazards Model ####
#' # with design weights; no calibration 
#' fit_ahaz_w = with(sample_data, ahaz.wcal(outcome = ind.fail, etime = eventime, tau = 8, 
#'                                  X = cbind(X1, X2), Z = cbind(Z1, Z2),
#'                                  ind.ph2 = ind.ph2, incl.prob = incl.prob, m = 1, nrisk = nrisk,
#'                                  A = NULL, eta = NULL))
#'                                  
#' # with calibrated weights with aux1
#' fit_ahaz_wcal_g = with(sample_data, ahaz.wcal(outcome = ind.fail, etime = eventime, tau = 8, 
#'                                       X = cbind(X1, X2), Z = cbind(Z1, Z2),
#'                                       ind.ph2 = ind.ph2, incl.prob = incl.prob, m = 1, nrisk = nrisk,
#'                                       A = aux1, eta = attr(cal1, "eta")))
#'                                       
#' # with calibrated weights with aux2
#' fit_ahaz_wcal_Bg = with(sample_data, ahaz.wcal(outcome = ind.fail, etime = eventime, tau = 8, 
#'                                        X = cbind(X1, X2), Z = cbind(Z1, Z2),
#'                                        ind.ph2 = ind.ph2, incl.prob = incl.prob, m = 1, nrisk = nrisk,
#'                                        A = aux2, eta = attr(cal2, "eta")))
#'                                        
#' #### Result: Estimates of Beta and gamma Coefficients ####
#' cbind(Estimate = fit_ahaz_w$est, StandardError = fit_ahaz_w$se)
#' cbind(Estimate = fit_ahaz_wcal_g$est, StandardError = fit_ahaz_wcal_g$se)
#' cbind(Estimate = fit_ahaz_wcal_Bg$est, StandardError = fit_ahaz_wcal_Bg$se)
#' 
#' 
#' #### Result: Estimates of Covariate-specific Pure Risks ####
#' XZ = (rbind(r1=c(0,0,0,0), r2=c(0,0,1,1), r3=c(0,1,0,1), r4=c(1,1,0,0), r5=c(1,0,1,0), r6=c(1,1,1,1)))
#' colnames(XZ) = c(paste("X", 1:2, sep=""), paste("Z", 1:2, sep=""))
#' 
#' t(apply(XZ, 1, prisk, fit = fit_ahaz_w, tau = 8))
#' t(apply(XZ, 1, prisk, fit = fit_ahaz_wcal_g, tau = 8))
#' t(apply(XZ, 1, prisk, fit = fit_ahaz_wcal_Bg, tau = 8))
#' 
prisk = function(fit, xzval, tau){
  p = sum(substr(names(fit$est), 1, 4) == "Beta") - 1
  q = sum(substr(names(fit$est), 1, 4) == "gamm")
  
  xval = xzval[1:p]; zval = xzval[(p+1):(p+q)]
  
  est = fit$est
  rhat = as.numeric(1-exp(-(est[1] + sum(est[(1+1):(1+p)]*xval) + tau * sum(est[(1+p+1):(1+p+q)]*zval))))
  
  dr_dest = c(1, xval, zval*tau)*(rhat-1)
  infl_prisk = fit$infl %*% dr_dest
  infl_omg_prisk = fit$infl.omg %*% dr_dest
  
  psig = crossprod(infl_prisk*fit$sqrtprob)
  pomg = crossprod(infl_omg_prisk[fit$ind.omega,], fit$jcov.byp) %*% infl_omg_prisk[fit$ind.omega,]
  
  rse = sqdg(psig+pomg)
  return(c(xzval, RiskEstimate = rhat, StandardError = rse))
}
