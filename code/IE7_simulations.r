library(fixest)
library(skimr)
library(sandwich)
library(lmtest)

rm(list = ls())
setwd("C:/Users/lucas/Dropbox/whales_own")

# Parameters
N <- 1000    # number of individuals
T <- 50     # number of time periods
set.seed(123)

# 1. Simulate data
mu <- rnorm(N, mean = 0, sd = 3)       # individual fixed effects
y_mat <- matrix(NA, nrow = N, ncol = T)
for (i in 1:N) {
  for (t in 1:T) {
    y_mat[i, t] <- mu[i] + rnorm(1)   # y_ct = α_i + u_ct
  }
}

# 2. Build long‐format data frame and create lag
df <- data.frame(
  id   = rep(1:N, each = T),
  time = rep(1:T, times = N),
  y    = as.vector(t(y_mat))
)
df <- df[order(df$id, df$time), ]

 
max_lag <- 5
for (i in 1:max_lag) {
  df[[paste0("lag", i, "_y")]] <- ave(df$y, df$id, FUN = function(x) {
    n <- length(x)
    c(rep(NA, i), x[1:(n - i)])
  })
}

df$avg_y <- ave(df$y, df$id, FUN = function(x) {
  out <- numeric(length(x))
  for(i in seq_along(x)) {
    out[i] <- if(i == 1) NA else mean(x[1:(i-1)])
  }
  out
})


# 3. Regressions
# (a) Without individual fixed effects
m1 <- feols(y ~ lag1_y, data = df)
summary(m1)

# (b) With individual fixed effects (via dummies)
#m2 <- lm(y ~ lag_y + factor(id), data = df)
#summary(m2)

m2 <- feols(y ~ 1 + lag1_y | id , data = df)
summary(m2)    


#it is weird that if you set all the captain FE to 0 you still get a negative
m3 <- feols(y ~ 1 + lag1_y + avg_y  | id , data = df)

m4 <- feols(y ~ 1 + lag1_y + lag2_y + avg_y  | id , data = df)
summary(m4)

m5 <- feols(y ~ 1 + lag1_y + lag2_y + lag3_y + avg_y  | id , data = df)
summary(m4)

 esttex(   m1, m2,   m3, m4, m5,  
    file   = "Writeup/Tables/simulations.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Simulations",
    digits = 3,
    notes = "" )