# library(mr.divw)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap)
# library(readxl)
# library(glmnet)


source("R/MR_rr_simulation.R")
source("R/MR_rr_estimators.R")
data("lip_data")
data("lip_corr")
data("lip_samplesize")

set.seed(123)


# non-parametric bootstrap method to estimate variance of the Drr estimator ----
#' Title
#'
#' @param me_weight a character value of a positive number, recommend to be chosen between 0.05 and 1. The argument indicates the measurement error weight, which means the scaler we multiply to the measurement error of the coefficients of exposure to SNPs, which is the \eqn{\Sigma_X} in the manuscript. Smaller measurement error weight implies larger IV strength.
#' @param regularized can be TRUE or FALSE, indicating whether we use the MR-rr estiomator with the spectral regularization. See the manuscript for more details.
#' @param regularization_rate a small positive numerical value, indicating the regularization rate. It is only used when regularized = TRUE. The larger regularization_rate make the estimator having smaller variance, but the rate can not be too large to ensure consistency holds, see more details in the manuscript.
#'
#' @return The average coverage rate of the entry-wise 95% confidence interval produced by non-parametric bootstrap method.
#' @export
#'
#' @examples
#' nonpara_bootstrap(me_weight = 0.2, regularized = TRUE, regularization_rate = 1e-13)
nonpara_bootstrap <- function(me_weight = 0.2, regularized = TRUE, regularization_rate = 1e-13){
  #### (a) get the parameters
  parameters = .get_parameters(me_weight = me_weight, px = 24, r_RR = 5)
  C = parameters$C
  A = parameters$A
  B = parameters$B
  A_d = parameters$A_d
  B_d = parameters$B_d
  C_r = parameters$C_r
  C_r_vec = as.vector(C_r)
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
  Sigma_Y_hat = parameters$Sigma_Y

  #### (b) in each iteration, gen data
  iteration = 100
  # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
  coverage_rate = matrix(NA, px*py, iteration)

  for (loop in 1:iteration){
    #### (b1) gen data in the trial
    simulation_result = .simulation(parameters, regularized = regularized, regularization_rate = regularization_rate)
    A_hat = simulation_result[[1]]
    B_hat = simulation_result[[2]]
    AB_hat = simulation_result[[3]]
    A_d_hat = simulation_result[[4]]
    B_d_hat = simulation_result[[5]]
    AB_d_hat = simulation_result[[6]]
    x_j_hat = simulation_result[[7]]
    y_j_hat = simulation_result[[8]]

    #### (b2) estimate V_X, C

    #### (b3) bootstrap
    bootstrap_size = 100  # sample size to calculate CI in each bt loop
    pz = 1000
    C_drr_hat_entrywise_bt_list = list()
    for (i in 1:(px*py)){
      C_drr_hat_entrywise_bt_list[[i]] = rep(0, bootstrap_size)
    }
    #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
    for (bt in 1:bootstrap_size){
      # sample bt_sample_number rows from x_j_hat, y_j_hat
      index_bt = sample(1:pz, pz, replace = TRUE)
      x_j_hat_star = x_j_hat[index_bt,]
      y_j_hat_star = y_j_hat[index_bt,]

      # estimate C_drr_hat
      if (regularized == FALSE) {
        result_d_bt <- mr_rr(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = .sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat)
      } else {
        result_d_bt <- mr_rr_regularized(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = .sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat, regularization_rate=regularization_rate)
      }
      result_d_bt_vec = as.vector(result_d_bt$AB)
      for (j in 1:(px*py)){
        C_drr_hat_entrywise_bt_list[[j]][bt] = result_d_bt_vec[j]
      }
    }

    #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
    C_drr_hat_entrywise_bt_95_list = list()
    for (i in 1:(px*py)){
      C_drr_hat_entrywise_bt_95_list[[i]] = quantile(C_drr_hat_entrywise_bt_list[[i]], c(0.025, 0.975))
    }
    #### (b5) calculate the coverage rate for each entry of C_r
    index_in95 = rep(1, px*py)
    for (i in 1:(px*py)){
      if (C_r_vec[i] < C_drr_hat_entrywise_bt_95_list[[i]][1] | C_r_vec[i] > C_drr_hat_entrywise_bt_95_list[[i]][2]){
        index_in95[i] = 0
      }
    }
    coverage_rate[,loop] = index_in95
  }

  #### (c) calculate the average coverage rate
  entry_coverage_rate = rowMeans(coverage_rate)
  average_coverage_rate = mean(entry_coverage_rate)
  average_coverage_rate
}


#####
## result:
# weight=0.2, regularized=TRUE, regularization_rate=1e-13 , average_coverage_rate=0.9065833
# weight=0.5, regularized=TRUE, regularization_rate=1e-13 , average_coverage_rate=0.9202917


# regularized = TRUE
# # parametric bootstrap method to estimate variance of the Drr estimator ----
# #### (a) get the parameters
# parameters = .get_parameters(me_weight = 0.2, px = 24, r_RR = 5)
# C = parameters$C
# A = parameters$ATRUE
# B = parameters$B
# A_d = parameters$A_d
# B_d = parameters$B_d
# C_r = parameters$C_r
# py = parameters$py
# px = parameters$px
# Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
# Sigma_Y_hat = parameters$Sigma_Y
#
# #### (b) in each iteration, gen data
# iteration = 100
# # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
# coverage_rate = matrix(NA, px*py, iteration)
#
# for (loop in 1:iteration){
#   #### (b1) gen data in the trial
#   simulation_result = .simulation(parameters, regularized = regularized)
#   A_hat = simulation_result[[1]]
#   B_hat = simulation_result[[2]]
#   AB_hat = simulation_result[[3]]
#   A_d_hat = simulation_result[[4]]
#   B_d_hat = simulation_result[[5]]
#   AB_d_hat = simulation_result[[6]]
#   x_j_hat = simulation_result[[7]]
#   y_j_hat = simulation_result[[8]]
#
#   #### (b2) estimate V_X, C
#   V_X_tilde_star = t(x_j_hat) %*% x_j_hat / 1000 - Sigma_X_hat
#   C_star = AB_d_hat
#   C_r_vec = as.vector(C_r)
#
#   #### (b3) bootstrap
#   bootstrap_size = 100  # sample size to calculate CI in each bt loop
#   pz = 1000
#   C_drr_hat_entrywise_bt_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_list[[i]] = rep(0, bootstrap_size)
#   }
#   #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
#   for (i in 1:bootstrap_size){
#     # sample true effect x_j_star and y_j_star
#     x_j_star = mvrnorm(n = pz, mu = rep(0, px), Sigma = V_X_tilde_star, tol = 100)
#     y_j_star = x_j_star %*% t(C)
#
#     # sample x_j_hat, y_j_hat
#     x_j_hat_star = matrix(0, pz, px)
#     y_j_hat_star = matrix(0, pz, py)
#     for (j in 1:pz) {
#       x_j_hat_star[j,] = mvrnorm(n = 1, mu = x_j_star[j,], Sigma = Sigma_X_hat, tol = 100)
#       y_j_hat_star[j,] = mvrnorm(n = 1, mu = y_j_star[j,], Sigma = Sigma_Y_hat, tol = 100)
#     }
#
#     # estimate C_drr_hat
#     if (regularized == FALSE) {
#       result_d_bt <- mr_rr(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = .sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat)
#     } else {
#       result_d_bt <- mr_rr_regularized(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = .sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat, regularization_rate=regularization_rate)
#     }
#     result_d_bt_vec = as.vector(result_d_bt$AB)
#     for (j in 1:(px*py)){
#       C_drr_hat_entrywise_bt_list[[j]][i] = result_d_bt_vec[j]
#     }
#   }
#
#   #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
#   C_drr_hat_entrywise_bt_95_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_95_list[[i]] = quantile(C_drr_hat_entrywise_bt_list[[i]], c(0.05, 0.95))
#   }
#   #### (b5) calculate the coverage rate for each entry of C_r
#   C_r_vec = as.vector(C_r)
#   index_in95 = rep(1, px*py)
#   for (i in 1:(px*py)){
#     if (C_r_vec[i] < C_drr_hat_entrywise_bt_95_list[[i]][1] | C_r_vec[i] > C_drr_hat_entrywise_bt_95_list[[i]][2]){
#       index_in95[i] = 0
#     }
#   }
#   coverage_rate[,loop] = index_in95
# }
#
# #### (c) calculate the average coverage rate
# entry_coverage_rate = rowMeans(coverage_rate)
# average_coverage_rate = mean(entry_coverage_rate)
# average_coverage_rate
# #####
