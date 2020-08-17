
#' Joint Probabilities of Nested Case-control Samples
#'
#' Internal use to compute joint probabilities of nested case-control samples. The source code is similar to 'multipleNCC' package.
#' @param outcome binary outcome status (1 = case; 0 = non-case)
#' @param etime time to event (either case or censored, whichever comes first)
#' @param tau risk projection time
#' @param m number of controls per case
#' @param incl.prob probability of being included in NCC; e.g., this can be computed from KMprob() of 'multipleNCC'; if NULL, it is internally computed based on outcome, etime, and m
#' @param nrisk = number of at-risk subjects at each 'etime'; if NULL, it is internally computed based on outcome, etime, and m
#' @keywords ahazcal
#' @export 
#' 
jointVP.ncc = function(outcome, etime, tau, m, incl.prob = NULL, nrisk = NULL){
  
  n = length(outcome)
  control.num = sum(outcome==0)
  case.times = etime[outcome == 1]
  control.times = etime[outcome == 0]
  
  if(is.null(nrisk)) {
    nrisk = mapply(function(i){ sum(etime >= etime[i]) }, 1:n)
  }
  if(is.null(incl.prob)) { 
    incl.prob = mapply(function(i){1-prod(1-m/(nrisk[which((etime < etime[i]) & (outcome == 1))]-1))}, 1:n)
  }
  incl.prob[outcome == 1] = 1
  
  nrisko1 = nrisk[outcome==1]
  H = ((1-2*m/(nrisko1-1)+m*(m-1)/((nrisko1-1)*(nrisko1-2)))/(1-m/(nrisko1-1))^2)
  
  Rho = matrix(0, ncol = control.num, nrow = control.num)
  for(i in (1:control.num-1)) {
    for(j in ((i+1):control.num))  {
      
      ind = 1:length(case.times)
      ind = ind[which(case.times<control.times[i] & case.times<control.times[j])]
      
      Rho[i,j] = prod(H[ind])-1
      
    }
  }
  
  Rhony = Rho + t(Rho)
  p0 = incl.prob[outcome == 0]
  
  Vij = Rhony * tcrossprod(1-p0) + diag((1-p0)*(p0)) 
  Pij = Rhony * tcrossprod(1-p0) + tcrossprod(p0); 
  diag(Pij) = p0 
  
  return(list(V = Vij, P = Pij))
}



#' Extract Standard Deviation from A Variance-covariance Matrix
#'
#' Internal use to compute the square roots of diagonals of a given variance-covariance matrix.
#' @param x A n x n variance-covariance matrix
#' @keywords ahazcal
#' @export 
#' 
sqdg = function(x){ sqrt(diag(x)) }

#' Weighted Matrix Crossproduct
#'
#' Internal use to compute the weighted crossproduct of one or two matrices. (x'wy)
#' @param x n x p matrix
#' @param y n x q matrix; if NULL, y = x
#' @param w n vector of weights
#' @keywords ahazcal
#' @export 
#'
wcrossprod = function(x, w, y = NULL){
  if(is.null(y)) y = x
  if(is.null(w)) w = rep(1, nrow(x))
  wy = sweep(y, 1, w, "*")
  return(crossprod(x, wy))
}