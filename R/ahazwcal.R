#' Fit Additive Hazards Model from Nested Case-control Samples (NCC)
#'
#' This function fits the additive hazards model with time-varying coefficients (beta) and time-invariant coefficients (gamma) from nested case-control study with or without weight calibration.
#' Weight calibration parameters and auxiliary statistics are required to fit the model with weight calibration. Otherwise, it will provide the result without weight calibration.
#' @param outcome binary outcome status (1 = case; 0 = non-case)
#' @param etime time-to-event (either case or censored, whichever comes first)
#' @param tau projection time for pure risk or beta coefficient
#' @param ind.ph2 a vector of binary status whether included in NCC (phase 2) or not
#' @param m number of controls per case
#' @param X a matrix of model covariates with time-varying effects (Beta); one or more columns may include NA if missing by NCC
#' @param Z a matrix of model covariates with time-invariant effects (gamma); one or more columns may include NA if missing by NCC
#' @param incl.prob probability of being included in NCC; e.g., this can be computed from KMprob() of 'multipleNCC'; if NULL, it is internally computed based on outcome, etime, and m
#' @param nrisk number of at-risk subjects at each 'etime'; if NULL, it is internally computed based on outcome, etime, and m
#' @param A a matrix of auxiliary statistics; e.g., cbind(1, AUX) where AUX is an object of aux.ahaz(); if A = NULL, it returns estimation results with design weights i.e., without weight calibration
#' @param eta a vector of estimates of weight calibration parameters obtained from the attributes of a rakecal() object; if A is NULL, eta must also be NULL
#' @keywords ahazcal
#' @return A list of beta-gamma estimates and their variance estimates with influence functions
#' @export 
#' @seealso \code{\link{sample_data}}, \code{\link{aux.ahaz}}, \code{\link{rakecal}}
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
ahaz.wcal = function(outcome, etime, tau, ind.ph2, m, X, Z, incl.prob = NULL, nrisk = NULL, A = NULL, eta = NULL){
  
  otime = sort(etime); otime1 = c(0, otime); dtime = diff(otime1)
  
  Yt = t(mapply(function(x){as.numeric(otime <= x)}, x = etime)); 
  N = nrow(Yt); nt = ncol(Yt)
  
  Xt = array(NA, dim = c(N, nt, ncol(X)))
  Zt = array(NA, dim = c(N, nt, ncol(Z)))
  for(j in 1:ncol(X)){Xt[,,j] = sweep(Yt, 1, X[,j], "*")}
  for(j in 1:ncol(Z)){Zt[,,j] = sweep(Yt, 1, Z[,j], "*")}
  
  evnt = apply(sweep(Yt, 1, outcome, "*"), 1, function(x){ifelse(all(x==0), 0, max(which(x==1)))})
  
  dNt = matrix(0, N, nt); dNt[cbind(1:N, evnt)] = 1; dNt = sweep(dNt, 1, outcome, "*")
  
  ind.ph2 = (ind.ph2 == 1); n = sum(ind.ph2)
  
  if(is.null(nrisk)) {
    nrisk = mapply(function(i){ sum(etime >= etime[i]) }, 1:N)
  }
  if(is.null(incl.prob)) { 
    incl.prob = mapply(function(i){1-prod(1-m/(nrisk[which((etime < etime[i]) & (outcome == 1))]-1))}, 1:N)
  }
  incl.prob[outcome == 1] = 1
  
  wgt = 1/incl.prob[ind.ph2]
  if(!is.null(A)) wgt = ((1/incl.prob) * exp(A %*% eta))[ind.ph2]
  
  XX = simplify2array(mapply(function(k){wcrossprod(cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,]), wgt)}, k=1:nt, SIMPLIFY = F))
  ZZ = simplify2array(mapply(function(k){wcrossprod(Zt[ind.ph2,k,], wgt)}, k=1:nt, SIMPLIFY = F))
  XZ = simplify2array(mapply(function(k){wcrossprod(cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,]), wgt, Zt[ind.ph2,k,])}, k = 1:nt, SIMPLIFY = F))
  XZdt = sweep(XZ, 3, dtime, "*")
  ZHZ = simplify2array(mapply(function(k){ZZ[,,k] - (crossprod(XZ[,,k], ginv(XX[,,k])) %*% XZ[,,k])}, k=1:nt, SIMPLIFY = F))
  ZdN = simplify2array(mapply(function(k){wcrossprod(Zt[ind.ph2,k,], wgt, dNt[ind.ph2,k,drop=F])}, k=1:nt, SIMPLIFY = F))
  XdN = simplify2array(mapply(function(k){wcrossprod(cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,]), wgt, dNt[ind.ph2,k,drop=F])}, k=1:nt, SIMPLIFY = F))
  ZHdN = simplify2array(mapply(function(k){ZdN[,,k] - (crossprod(XZ[,,k], ginv(XX[,,k])) %*% XdN[,,k])}, k=1:nt, SIMPLIFY = F))
  
  intZHZ = apply(sweep(ZHZ, 3, dtime, "*"), c(1,2), sum)
  intZHdN = apply(ZHdN, c(1,2), sum)
  
  #### ESTIMATES ####
  gam.hat = solve(intZHZ, intZHdN)
  
  dBt.hat = t(mapply(function(k){ginv(XX[,,k]) %*% (XdN[,,k] - XZdt[,,k] %*% gam.hat)}, k=1:nt, SIMPLIFY = T))
  Bt.hat = cbind(time = otime, apply(dBt.hat, 2, cumsum))
  
  #### MARTINGALE RESIDUALS ####
  XdBt = mapply(function(k){cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,]) %*% dBt.hat[k,]}, k=1:nt, SIMPLIFY = T) 
  Zgamdt = mapply(function(k){Zt[ind.ph2,k,] %*% gam.hat * dtime[k]}, k = 1:nt)
  
  dMt = sweep(dNt[ind.ph2,] - XdBt - Zgamdt, 1, wgt, "*")
  
  #### INFLUENCE FUNCTIONS ####
  ZHt = simplify2array(mapply(function(k){crossprod(tcrossprod(ginv(XX[,,k]), cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,])), XZ[,,k])},
                              k = 1:nt, SIMPLIFY = F))
  intZHdMi = apply(simplify2array(mapply(function(k){sweep(Zt[ind.ph2,k,]-ZHt[,,k], 1, dMt[,k], "*")}, k = 1:nt, SIMPLIFY = F)), c(1,2), sum)
  
  if.gam = if.gam.copy = matrix(0, nrow = N, ncol = ncol(Z))
  if0.dBt = if.dBt = array(0, dim=c(N, 1+ncol(X), nt))
  if.gam[ind.ph2,] = if.gam.copy[ind.ph2,] = intZHdMi %*% ginv(intZHZ)
  
  XdM = simplify2array(mapply(function(k){sweep(cbind(Yt[ind.ph2,k], Xt[ind.ph2,k,]), 1, dMt[,k], "*")}, k=1:nt, SIMPLIFY = F))
  res.dBt.gam = simplify2array(mapply(function(k){tcrossprod(if.gam, XZdt[,,k])}, k = 1:nt, SIMPLIFY = F))
  if0.dBt[ind.ph2,,] = simplify2array(mapply(function(k){XdM[,,k] %*% ginv(XX[,,k])}, k=1:nt, SIMPLIFY = F))
  if1.dBt.gam = simplify2array(mapply(function(k){res.dBt.gam[,,k] %*% ginv(XX[,,k])}, k = 1:nt, SIMPLIFY = F))
  if.dBt[ind.ph2,,] = simplify2array(mapply(function(k){(XdM[,,k]-res.dBt.gam[ind.ph2,,k]) %*% ginv(XX[,,k])}, k=1:nt, SIMPLIFY = F))
  
  id.tau = max(which(otime<tau))
  
  #### VARIANCE FOR WEIGHT-CALIBRATED ESTIMATES ####
  if.eta.omg = matrix(0, nrow = N, 2)
  ZA = matrix(0, ncol(Z), 2); XA = array(0, dim=c(1+ncol(X), 2, nt))
  
  if(!is.null(A)){
    AA = wcrossprod(A[ind.ph2,], wgt)
    if.eta.num = -A 
    if.eta.num[ind.ph2,] = if.eta.num[ind.ph2,] + sweep(A[ind.ph2,], 1, wgt,"*")
    if.eta = if.eta.num %*% ginv(AA)
    
    ZA = crossprod(intZHdMi, A[ind.ph2,])
    XA = simplify2array(mapply(function(k){crossprod(XdM[,,k], A[ind.ph2,])}, k = 1:nt, SIMPLIFY = F))
    
    if.gam = if.gam.copy - if.eta %*% crossprod(ZA, ginv(intZHZ))
    if1.dBt.gam = simplify2array(mapply(function(k){tcrossprod(if.gam, ginv(XX[,,k]) %*% XZdt[,,k])}, k = 1:nt, SIMPLIFY = F))
    if.dBt = if0.dBt - if1.dBt.gam - simplify2array(mapply(function(k){if.eta %*% crossprod(XA[,,k], ginv(XX[,,k]))}, k = 1:nt, SIMPLIFY = F))
    
    if.eta.num.omg = matrix(0, nrow = N, ncol = ncol(A))
    if.eta.num.omg[ind.ph2,] = sweep(A[ind.ph2,], 1, wgt,"*")
    if.eta.omg = if.eta.num.omg %*% ginv(AA)
  }
  
  if.gam.omg = if.gam.copy - if.eta.omg %*% crossprod(ZA, ginv(intZHZ))
  
  if1.dBt.gam.omg = simplify2array(mapply(function(k){tcrossprod(if.gam.omg, ginv(XX[,,k]) %*% XZdt[,,k])}, k = 1:nt, SIMPLIFY = F))
  if.dBt.omg = if0.dBt - if1.dBt.gam.omg - simplify2array(mapply(function(k){if.eta.omg %*% crossprod(XA[,,k], ginv(XX[,,k]))}, k = 1:nt, SIMPLIFY = F))
  
  if.Bt.tau.omg = apply(if.dBt.omg[,,1:id.tau], c(1,2), sum)
  infl.omg = cbind(if.Bt.tau.omg, if.gam.omg)
  
  
  Bt.hat.tau = Bt.hat[id.tau,]
  if0.Bt.tau = apply(if0.dBt[,,1:id.tau], c(1,2), sum)
  if.Bt.tau = apply(if.dBt[,,1:id.tau], c(1,2), sum)
  
  est = c(Bt.hat.tau[-1], gam.hat)
  infl = cbind(if.Bt.tau, if.gam)
  
  # VARIANCE 
  sqrt.prob = ifelse(ind.ph2, sqrt(incl.prob), 1)
  sig = crossprod(infl*sqrt.prob)
  omega = matrix(0, nrow(sig), ncol(sig)); ind.omega = jcov.byp = NULL
  if(!is.null(m)){
    VP = jointVP.ncc(outcome = outcome[ind.ph2], etime = etime[ind.ph2], tau = tau, m = m, 
                     nrisk = nrisk[ind.ph2], incl.prob = incl.prob[ind.ph2]); 
    jcov.byp = with(VP, V/P)
    ind.omega = (ind.ph2) & (outcome == 0)
    omega = crossprod(infl.omg[ind.omega,], jcov.byp) %*% infl.omg[ind.omega,]
  }
  
  se = sqdg((sig+omega)*N/(N-1))[1:(1+ncol(X)+ncol(Z))]
  
  names(est) = names(se) = colnames(infl) = c(paste("Beta", 0:ncol(X), sep=""), paste("gamma", 1:ncol(Z), sep=""))
  
  return(list(est = est, se = se,
              infl = infl, infl.omg = infl.omg,
              sqrt.prob = sqrt.prob, ind.omega = ind.omega, jcov.byp = jcov.byp))
}

