## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(Romeb)

## ----example-complete, eval = TRUE--------------------------------------------
set.seed(123)
Y <- matrix(rnorm(300*5), nrow = 300, ncol = 5)  # tiny complete data set
result_full <- Romeb("no missing", data = Y, Time = 5, seed = 123)
print(result_full)

## ----example-mcar, eval = FALSE-----------------------------------------------
# set.seed(456)
# Y <- matrix(rnorm(300 * 5), nrow = 300)
# miss <- runif(length(Y)) < 0.1       # 10% missing completely at random
# Y[miss] <- NA
# result_mcar <- Romeb("MCAR", data = Y, Time = 5, seed = 456)

## ----example-mnar-aux, eval = FALSE-------------------------------------------
# set.seed(789)
# X  <- matrix(rnorm(300 * 2), 300, 2)     # two auxiliaries
# Y  <- matrix(rnorm(300 * 5), 300, 5)
# Data <- cbind(X, Y)
# result_mnar <- Romeb("MNAR", data = Data, Time = 5, K = 2, seed = 789)

## ----traceplot, fig.cap = "Trace plot for the first chain (par[1])", eval = TRUE----
# Uses the tiny example result_full from above
coda::traceplot(result_full$samps_full[,'par[1]'])

