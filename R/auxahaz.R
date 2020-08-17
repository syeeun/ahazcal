#' Auxiliary Statistics of Pure Risk Estimations from Additive Hazards Model
#'
#' This function provides influence functions of additive hazards model coefficients: beta are time-varying coefficients of X; and gamma are time-invariant coefficients of Z. 
#' The resulting values with observed/predicted covariates in a full cohort can be used as auxiliary statistics of pure risk estimates for weight calibration.
#' All arguments for this function except tau must be complete with the same length of vectors.
#' @param outcome binary outcome status (1 = case; 0 = non-case)
#' @param etime time to event (either case or censored, whichever comes first)
#' @param tau risk projection time
#' @param X a matrix of observed or predicted model covariates with time-varying effects (Beta)
#' @param Z a matrix of observed or predicted model covariates with time-invariant effects (gamma)
#' @keywords ahazcal
#' @return a matrix of auxiliary statistics of pure risk estimates (cumulated influence functions of beta and influence functions of gamma) to calibrate design weights
#' @export 
#' @seealso \code{\link{sample_data}}
#' @examples #### STEP1: Predict Missing Covariates ####
#' fit_Z2 = glm(Z2 ~ X1 + X2 + Z1 + U, family = binomial(link = "logit"), 
#'              data = sample_data, subset = (ind.ph2 == 1), weights = 1/incl.prob)
#' 
#' sample_data.tilde = sample_data
#' sample_data.tilde$Z2 = predict(fit_Z2, type = "response", newdata = sample_data)
#' 
#' 
#' #### STEP2: Obtain Auxiliary Statistics (Influence Functions) ####
#' aux = with(sample_data.tilde, aux.ahaz(outcome = ind.fail, etime = eventime, tau = 8, 
#'                                        X = cbind(X1, X2), Z = cbind(Z1, Z2)))
#' 
#' aux1 = cbind(1, aux[,c("gamma1","gamma2")]) # Auxiliary 1: Influence functions of gamma estimates 
#' aux2 = cbind(1, aux) # Auxiliary 2: Influence functions of Beta and gamma estimates
#' 
aux.ahaz = function(outcome, etime, tau, X, Z){
  
  otime = sort(etime); otime1 = c(0, otime); dtime = diff(otime1)
  
  Yt = t(mapply(function(x){as.numeric(otime <= x)}, x = etime)); 
  n = nrow(Yt); nt = ncol(Yt)
  
  Xt = array(NA, dim = c(n, nt, ncol(X)))
  Zt = array(NA, dim = c(n, nt, ncol(Z)))
  for(j in 1:ncol(X)){Xt[,,j] = sweep(Yt, 1, X[,j], "*")}
  for(j in 1:ncol(Z)){Zt[,,j] = sweep(Yt, 1, Z[,j], "*")}
  
  evnt = apply(sweep(Yt, 1, outcome, "*"), 1, function(x){ifelse(all(x==0), 0, max(which(x==1)))})
  
  dNt = matrix(0, n, nt); dNt[cbind(1:n, evnt)] = 1; dNt = sweep(dNt, 1, outcome, "*")
  
  XX = simplify2array(mapply(function(k){crossprod(cbind(Yt[,k], Xt[,k,]))}, k=1:nt, SIMPLIFY = F))
  ZZ = simplify2array(mapply(function(k){crossprod(Zt[,k,])}, k=1:nt, SIMPLIFY = F))
  XZ = simplify2array(mapply(function(k){crossprod(cbind(Yt[,k], Xt[,k,]), Zt[,k,])}, k = 1:nt, SIMPLIFY = F))
  XZdt = sweep(XZ, 3, dtime, "*")
  ZHZ = simplify2array(mapply(function(k){ZZ[,,k] - (crossprod(XZ[,,k], ginv(XX[,,k])) %*% XZ[,,k])}, k=1:nt, SIMPLIFY = F))
  ZdN = simplify2array(mapply(function(k){crossprod(Zt[,k,], dNt[,k,drop=F])}, k=1:nt, SIMPLIFY = F))
  XdN = simplify2array(mapply(function(k){crossprod(cbind(Yt[,k], Xt[,k,]), dNt[,k,drop=F])}, k=1:nt, SIMPLIFY = F))
  ZHdN = simplify2array(mapply(function(k){ZdN[,,k] - (crossprod(XZ[,,k], ginv(XX[,,k])) %*% XdN[,,k])}, k=1:nt, SIMPLIFY = F))
  
  intZHZ = apply(sweep(ZHZ, 3, dtime, "*"), c(1,2), sum)
  intZHdN = apply(ZHdN, c(1,2), sum)
  
  #### ESTIMATES ####
  gam.hat = solve(intZHZ, intZHdN)
  
  dBt.hat = t(mapply(function(k){ginv(XX[,,k]) %*% (XdN[,,k] - XZdt[,,k] %*% gam.hat)}, k=1:nt, SIMPLIFY = T))
  Bt.hat = cbind(time = otime, apply(dBt.hat, 2, cumsum))
  
  #### MARTINGALE RESIDUALS ####
  XdBt = mapply(function(k){cbind(Yt[,k], Xt[,k,]) %*% dBt.hat[k,]}, k=1:nt, SIMPLIFY = T) 
  Zgamdt = mapply(function(k){Zt[,k,] %*% gam.hat * dtime[k]}, k = 1:nt)
  
  dMt = dNt - XdBt - Zgamdt
  
  #### INFLUENCE FUNCTIONS ####
  ZHt = simplify2array(mapply(function(k){crossprod(tcrossprod(ginv(XX[,,k]), cbind(Yt[,k], Xt[,k,])), XZ[,,k])},
                              k = 1:nt, SIMPLIFY = F))
  intZHdMi = apply(simplify2array(mapply(function(k){sweep(Zt[,k,]-ZHt[,,k], 1, dMt[,k], "*")}, k = 1:nt, SIMPLIFY = F)), c(1,2), sum)
  
  if.gam = if.gam.copy = intZHdMi %*% ginv(intZHZ)
  
  XdM = simplify2array(mapply(function(k){sweep(cbind(Yt[,k], Xt[,k,]), 1, dMt[,k], "*")}, k=1:nt, SIMPLIFY = F))
  res.dBt.gam = simplify2array(mapply(function(k){tcrossprod(if.gam, XZdt[,,k])}, k = 1:nt, SIMPLIFY = F))
  if0.dBt = simplify2array(mapply(function(k){XdM[,,k] %*% ginv(XX[,,k])}, k=1:nt, SIMPLIFY = F))
  if1.dBt.gam = simplify2array(mapply(function(k){res.dBt.gam[,,k] %*% ginv(XX[,,k])}, k = 1:nt, SIMPLIFY = F))
  if.dBt = simplify2array(mapply(function(k){(XdM[,,k]-res.dBt.gam[,,k]) %*% ginv(XX[,,k])}, k=1:nt, SIMPLIFY = F))
  
  id.tau = max(which(otime<tau))
  
  aux = cbind(apply((XdM-res.dBt.gam)[,,1:id.tau], c(1,2), sum), intZHdMi)
  colnames(aux) = c(paste("Beta", 0:ncol(X), sep=""), paste("gamma", 1:ncol(Z), sep=""))
  return(aux)
}

