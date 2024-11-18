# # new sim for sparse 10.30 ----
# soft_threshold <- function(x, lambda) {
#   sign(x) * pmax(abs(x) - lambda, 0)
# }
#
#
# parameters = get_parameters(me_weight=0.1, px = 24, r_RR = 5)
# pz = n = 100 #(pz)
# lambda = rep(1, parameters$px)
#
# # generate B sparse para
# py = parameters$py
# px = parameters$px
# var_Z = parameters$var_Z
# VX_tilde = parameters$VX_tilde
# Sigma_X = parameters$Sigma_X
# Sigma_Y = parameters$Sigma_Y
# SigmaXX = parameters$SigmaXX
# SigmaYX = parameters$SigmaYX
# SigmaXY = parameters$SigmaXY
# C = parameters$C
# r_RR = parameters$r_RR
# VY_tilde = parameters$VY_tilde
# W = parameters$weight.matrix
# W_sqrt = .sqrt_matrix(parameters$weight.matrix)
# true.A_sparse = parameters$true.A_sparse
# true.B_sparse = parameters$true.B_sparse
#
# # # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
# # gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
# # Gamma_j_star = gamma_j_star %*% t(C)
# #
# # # sample x_j_hat, y_j_hat
# # x_j_hat = matrix(0, n, px)
# # y_j_hat = matrix(0, n, py)
# # for (j in 1:n) {
# #   x_j_hat[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
# #   y_j_hat[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
# # }
#
# # compute sample cov matrix
# hat_Sigma_xhat_xhat = t(x_j_hat) %*% x_j_hat / n
# hat_Sigma_x_x = hat_Sigma_xhat_xhat - Sigma_X
#
# # check if PSD
# is_psd <- function(matrix) {
#   eigenvalues <- eigen(matrix)$values
#   all(eigenvalues >= 0)
# }
#
#
# is_psd(hat_Sigma_x_x)
#
# # start to implement algorithm
# # step 1: debiasing the objective function
# # chol decomposition
# tilde_gamma = t(chol(hat_Sigma_x_x)) * sqrt(n)
# # check if chol is correct
# all.equal(tilde_gamma %*% t(tilde_gamma) / n, hat_Sigma_x_x)
#
# tilde_Gamma = t(solve(tilde_gamma)%*%t(x_j_hat)%*%y_j_hat)
#
# # check the definition of tilde_Gamma
# all.equal(tilde_gamma %*% t(tilde_Gamma), t(x_j_hat)%*%y_j_hat)
#
# # step 2: solve the objective function
# # store iteration values of A and B
# A_hat_list = list()
# B_hat_list = list()
# # initial values of A_hat and B_hat using dRRR
# result <- mr_rr(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X)
# A_hat = result$A
# B_hat = result$B
# A_hat_list[[1]] = A_hat
# B_hat_list[[1]] = B_hat
# t(A_hat) %*% W %*% A_hat
#
#
# # given B minimize A:  Orthogonal Procrustes ----
# U = svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)$u
# V = svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)$v
#
# A_hat = solve(W_sqrt) %*% V %*% t(U)
# # regularization, make t(A_hat) %*% W %*% A_hat= I
#
# A_hat_list = c(A_hat_list, list(A_hat))
# #####
#
# # given A minimize B:  Cyclic Coordinate Descent ----
#
# # 计算目标函数的第一部分
# first_term <- (1 / n) * sum(sapply(1:px, function(j) {
#   W_j <- W[, j, drop = FALSE]       # 取矩阵 W 的第 j 列
#   gamma_tilde_j <- tilde_gamma[, j, drop = FALSE] # 取矩阵 tilde_gamma 的第 j 列
#   norm(t(A_hat) %*% W_j - B_hat %*% gamma_tilde_j, type = "2")^2
# }))
#
# # 计算目标函数的第二部分
# second_term <- sum(sapply(1:px, function(k) {
#   lambda[k] * sum(abs(B_hat[, k])) # 对矩阵 B_hat 的第 k 列计算 L1 范数
# }))
#
# # 最终目标函数的值
# objective_function_value <- first_term + second_term
#
#
#
# # coordinate decent
# # tilde_gamma = x_j_hat
# # tilde_Gamma = t(y_j_hat)
# # W = diag(1, py)
# max_change = 1000
# iter = 1
# # B_hat = matrix(0, r_RR, px)
# B_hat_list[[1]] = B_hat
# # while iter too large also stop
# while (iter < 1000) {
#   for (k in 1:px){
#     # for given k, compute all r_j, j=1,...,p_Y
#     updata_k = 0
#     for (j in 1:px){ # px
#       # compute the summation of the ith column vect of B_hat times the ji entry of tilde gamma
#       sum = 0
#       for (i in 1:px){
#         if (i != k){
#           sum = sum + B_hat[,i] * tilde_gamma[j,i]
#         }
#       }
#       r_j = t(A_hat) %*% W_sqrt %*% W_sqrt %*% tilde_Gamma[,j] - sum
#       r_j_times_tilde_gamma = r_j * tilde_gamma[j,k]
#       updata_k = updata_k + r_j_times_tilde_gamma
#     }
#     updata_k = updata_k/pz
#     # update B_hat
#     B_hat[,k] = soft_threshold(B_hat[,k] + updata_k, lambda[k])
#   }
#   mean_change = mean(abs(B_hat-B_hat_list[[iter]]))
#   B_hat_list[[iter+1]] = B_hat
#   iter = iter +1
#   # print iter and max_change
#   print(c(iter, mean_change))
# }
#
#
# # TODO: test extreme case
#
#
#
#
# # we should see the change of objective function value
#
#
#
# # standard coordinate decent (univariate y)----
# # Implementation of Cyclic Coordinate Descent for solving Lasso
#
# # Soft-thresholding operator
# soft_threshold <- function(x, lambda) {
#   return(sign(x) * pmax(abs(x) - lambda, 0))
# }
#
# # Generate some example data
# set.seed(123)
# X <- matrix(rnorm(100 * 10), 100, 10)  # 100 samples, 10 predictors
# y <- rnorm(100)  # response vector
# lambda <- 0.1  # regularization parameter
#
# # Initialize parameters
# beta <- rep(0, ncol(X))  # initial coefficients set to zero
# max_iter <- 1000  # maximum number of iterations
# tol <- 1e-6  # convergence tolerance
#
# # Coordinate Descent Algorithm
# for (iter in 1:max_iter) {
#   beta_old <- beta
#   for (j in 1:ncol(X)) {
#     # Compute partial residual without the contribution from the j-th predictor
#     r_j <- y - X[, -j] %*% beta[-j]
#     # Update the j-th coefficient using the soft-thresholding operator
#     beta[j] <- soft_threshold(mean(X[, j] * r_j), lambda)
#   }
#   # Check for convergence
#   if (max(abs(beta - beta_old)) < tol) {
#     cat("Algorithm converged after", iter, "iterations.\n")
#     break
#   }
# }
#
# # Print the final coefficients
# cat("Estimated coefficients:", beta, "\n")

# old ver: group lasso----

# may not apply, this is group lasso
# AT_W_TildeG = t(A_hat) %*% W %*% tilde_Gamma
# fit <- glmnet(t(tilde_gamma), t(AT_W_TildeG), family = "mgaussian", alpha = 1, intercept = FALSE)
# # output B hat as a matrix
# # empty
# B_hat = matrix(0, px, py)
# for (i in length(coef(fit, s = 0.1))){
#   if (i == 1){
#     B_hat = coef(fit, s = 0.1)$y[,1][-1]
#   } else {
#     B_hat = cbind(B_hat, coef(fit, s = 0.1)$beta[,i][-1])
#   }
# }
#
# beta_final <- lapply(fit$beta, function(mat) mat[, ncol(mat)])
# B_hat <- t(do.call(cbind, beta_final))
# B_hat_list = c(B_hat_list, B_hat)

#####
