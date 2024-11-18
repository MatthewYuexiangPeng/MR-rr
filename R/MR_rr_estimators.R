#' Naive MR-rr estimator
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param sqrt_Gamma a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the square root of the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param sqrt_Gamma_inv the inverse of sqrt_Gamma.
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the naive MR-rr estimator, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the naive MR-rr estimator. See more details in the manuscript.
#' @export
mr_rr_naive = function(Y, X, r, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {

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


#' MR-rr estimator
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param Sigma_X a px by px numeric matrix, the average conditional covariance matrix of the coefficients of regressing exposures on SNPs, which can be culculated by averaging the regression SEs across all SNPs. It corresponds to the \eqn{\Sigma_X} in the manuscript, see more details in the estimator section of the manuscript.
#' @param sqrt_Gamma a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the square root of the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param sqrt_Gamma_inv the inverse of sqrt_Gamma.
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the MR-rr estimator, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the MR-rr estimator. See more details in the manuscript.
#' @export
mr_rr = function(Y, X, r, Sigma_X, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {

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


#' MR-rr estimator with spectral regularization
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param Sigma_X a px by px numeric matrix, the average conditional covariance matrix of the coefficients of regressing exposures on SNPs, which can be culculated by averaging the regression SEs across all SNPs. It corresponds to the \eqn{\Sigma_X} in the manuscript, see more details in the estimator section of the manuscript.
#' @param regularization_rate a numeric value, the regularization rate used in the MR-rr estimator with spectral regularization. It corresponds to the \eqn{\phi} in the manuscript, see more details in the estimator section of the manuscript.
#' @param sqrt_Gamma a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the square root of the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param sqrt_Gamma_inv the inverse of sqrt_Gamma.
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the MR-rr estimator with spectral regularization, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the MR-rr estimator with spectral regularization. See more details in the manuscript.
#' @export
mr_rr_regularized = function(Y, X, r, Sigma_X, regularization_rate = 1e-13, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {

  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  stable_Sigma_xx = debiased_Sigma_xx + regularization_rate * debiased_Sigma_xx_inv
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
