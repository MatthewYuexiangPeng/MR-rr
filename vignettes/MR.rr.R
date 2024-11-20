## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(MR.rr)
set.seed(123)

sim_result = run_simulation(regularized = TRUE, regularization_rate = 1e-13)
bias_naive_MRrr = sim_result[[1]]
bias_MRrr_regularized = sim_result[[2]]
print("This is a demo of the output format of the run_simulation function:")
str(bias_naive_MRrr)
print("This is a summary of the entrywise bias of the naive MR-rr estimator:")
summary(unlist(bias_naive_MRrr))
print("This is a summary of the entrywise bias of the MR-rr estimator with regularization:")
summary(unlist(bias_MRrr_regularized))

