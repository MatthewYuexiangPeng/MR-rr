#' Title
#'
#' @param Y 
#' @param X 
#' @param r 
#' @param sqrt_Gamma 
#' @param sqrt_Gamma_inv 
#'
#' @return
#' @export
#'
#' @examples
RRR = function(Y, X, r, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {
  
  XtY = crossprod(X, Y)
  XtX_inv = solve(crossprod(X))
  
  if (is.null(sqrt_Gamma)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  }
  
  V = sqrt_Gamma %*% t(XtY) %*% XtX_inv %*% XtY %*% sqrt_Gamma / nrow(Y)
  V = eigen(V)$vec[, 1:r, drop = FALSE]
  
  if (is.null(sqrt_Gamma_inv)) sqrt_Gamma_inv = solve(sqrt_Gamma)
  
  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(XtY) %*% XtX_inv)
  
  return(list(A = A, B = B, AB = A %*% B))
  
}


#' Title
#'
#' @param Y 
#' @param X 
#' @param r 
#' @param Sigma_X 
#' @param sqrt_Gamma 
#' @param sqrt_Gamma_inv 
#'
#' @return
#' @export
#'
#' @examples
debiased_RRR = function(Y, X, r, Sigma_X, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {
  
  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  
  if (is.null(sqrt_Gamma)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  }
  
  V = sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv %*% sigmaxy %*% sqrt_Gamma
  V = eigen(V)$vec[, 1:r, drop = FALSE]
  
  if (is.null(sqrt_Gamma_inv)) sqrt_Gamma_inv = solve(sqrt_Gamma)
  
  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv)
  
  return(list(A = A, B = B, AB = A %*% B))
  
}


#' Title
#'
#' @param Y 
#' @param X 
#' @param r 
#' @param Sigma_X 
#' @param ph 
#' @param sqrt_Gamma 
#' @param sqrt_Gamma_inv 
#'
#' @return
#' @export
#'
#' @examples
debiased_RRR_modified = function(Y, X, r, Sigma_X, ph = NULL, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {
  
  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  stable_Sigma_xx = debiased_Sigma_xx + ph * debiased_Sigma_xx_inv
  stable_Sigma_xx_inv = solve(stable_Sigma_xx)
  
  
  if (is.null(sqrt_Gamma)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  }
  
  V = sqrt_Gamma %*% t(sigmaxy) %*% stable_Sigma_xx_inv %*% sigmaxy %*% sqrt_Gamma
  V = eigen(V)$vec[, 1:r, drop = FALSE]
  
  if (is.null(sqrt_Gamma_inv)) sqrt_Gamma_inv = solve(sqrt_Gamma)
  
  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(sigmaxy) %*% stable_Sigma_xx_inv)
  
  return(list(A = A, B = B, AB = A %*% B))
  
}