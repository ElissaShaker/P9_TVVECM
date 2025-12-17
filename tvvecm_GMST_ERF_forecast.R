# source("~/Desktop/P9/P9_TimeVaryingCointegration/Radiative_forcings.r")

p # lag order
m <- 5 # Chebyshev polynomials
k # number of series
r # cointegration rank

result_final <- tvcoint(y, p, m) # building the tvvecm model

beta_t <- betat[, (k * (m - 1) + 1):(k * m)] # time-varying beta
resu_hat  <- result_final$resu    # gives U_t, residuals
betau_hat <- result_final$betau   # parameter estimates for short-run eq.

Gamma_1 <- betau_hat[2:(k+1), ] # since we have p=2 we only need gamma1
mu <- betau_hat[1, ]    # intercept in Δy-equation (epsilon_t)

evect <- result_final$eigenvectors   # gives q
ec_vecs <- result_final$eigenvectors[1:k, 1:r]
S01_baseline <- result_final$S01[, 1:k]   # extract first 3 columns
alpha <- result_final$S00inv %*% S01_baseline %*% ec_vecs # time invariant alpha

# forecast
tvvecm_forecast <- function(y, beta_t, alpha, Gamma_1, mu, h = 5) {
  T_y <- nrow(y)
  k   <- ncol(y)
  
  alpha <- as.matrix(alpha)
  Gamma_1 <- as.matrix(Gamma_1)
  mu <- matrix(mu, nrow = k, ncol = 1)
  
  beta_last <- as.matrix(beta_t[nrow(beta_t), ])
  if (ncol(beta_last) != ncol(alpha)) {
    beta_last <- matrix(beta_last, nrow = k, ncol = ncol(alpha))
  }
  
  y_last  <- matrix(as.numeric(y[T_y, ]), ncol = 1)
  dy_last <- matrix(as.numeric(y[T_y, ] - y[T_y-1, ]), ncol = 1)
  
  y_fore <- matrix(NA, nrow = h, ncol = k)   # h forecast steps only
  colnames(y_fore) <- colnames(y)
  
  for (i in 1:h) {
    ec_t <- t(beta_last) %*% y_last
    dy_new <- mu + alpha %*% ec_t + Gamma_1 %*% dy_last
    y_new <- y_last + dy_new
    
    y_fore[i, ] <- as.numeric(y_new)
    
    y_last <- y_new
    dy_last <- dy_new
  }
  
  return(y_fore)
}

y_pred <- tvvecm_forecast(y, beta_t, alpha, Gamma_1, mu, h = 5) 

#############################
# confidence interval
tvvecm_forecast_ci <- function(y, y_fore, resu_hat, level = 0.95) {
  
  T_y <- nrow(y)
  h   <- nrow(y_fore)
  k   <- ncol(y)
  
  Sigma_u <- cov(resu_hat)
  z <- qnorm(1 - (1 - level) / 2)
  
  lower <- upper <- matrix(NA, nrow = h, ncol = k)
  
  for (i in 1:h) {
    sd_h <- sqrt(diag(i * Sigma_u))
    lower[i, ] <- y_fore[i, ] - z * sd_h
    upper[i, ] <- y_fore[i, ] + z * sd_h
  }
  
  colnames(lower) <- colnames(y)
  colnames(upper) <- colnames(y)
  
  list(lower = lower, upper = upper)
}

ci <- tvvecm_forecast_ci(y, y_pred, resu_hat)

#############################
# plot
plot_tvvecm_forecast <- function(y, y_fore, ci, varname, vartext, years = merged_data$year) {
  
  T_y <- nrow(y)
  h   <- nrow(y_fore)
  
  # Observed years
  time_obs <- years
  
  # Forecast years: continue the sequence
  time_fc <- seq(from = max(years) + 1, by = 1, length.out = h)
  
  # Extract series
  y_obs <- y[[varname]]
  y_fc    <- y_fore[, varname]
  y_lower <- ci$lower[, varname]
  y_upper <- ci$upper[, varname]
  
  plot(
    time_obs, y_obs,
    type = "l",
    lwd  = 2,
    xlim = c(min(time_obs), max(time_fc)),
    ylim = range(c(y_obs, y_lower, y_upper)),
    xlab = "Year",
    ylab = vartext,
    main = paste("TV-VECM", h , "year forecast for", vartext)
  )
  
  lines(time_fc, y_fc, lwd = 2, lty = 2)
  lines(time_fc, y_lower, lty = 3)
  lines(time_fc, y_upper, lty = 3)
  
  abline(v = max(time_obs), lty = 3)
  
  legend(
    "topleft",
    legend = c("Observed", "Forecast", "95% CI"),
    lty    = c(1, 2, 3),
    lwd    = c(2, 2, 1),
    bty    = "n"
  )
}

plot_tvvecm_forecast(y, y_pred, ci, "GMST", "GMST")



### model validation
merged_data_lille <- merged_data[merged_data$year<=2019, ]
merged_data_test <- merged_data[merged_data$year > 2019, ]
y1_l <- merged_data_lille$GMST          # response
y2_l <- merged_data_lille$ERF           # forcing

# Bind into a matrix
y_l <- as.data.table(cbind(y1_l, y2_l))

y_test <- as.data.table(cbind(
  y1 = merged_data_test$GMST,
  y2 = merged_data_test$ERF
))

# Training data
y_l <- as.data.table(cbind(
  y1 = merged_data_lille$GMST,
  y2 = merged_data_lille$ERF
))

# Test data (realized values)
y_test <- as.data.table(cbind(
  y1 = merged_data_test$GMST,
  y2 = merged_data_test$ERF
))

ci_l <- tvvecm_forecast_ci(y_l, y_pred_l, resu_hat)
y_pred_l <- tvvecm_forecast(y_l, beta_t, alpha, Gamma_1, mu, h = 5)


plot_tvvecm_forecast_oos <- function(y_train, y_test, y_fore, ci, varname, 
                                     vartext, years_train, years_test){
  
  h <- nrow(y_fore)
  
  # Time indices
  time_obs <- years_train
  time_fc  <- years_test
  
  # Series
  y_obs   <- y_train[[varname]]
  y_real  <- y_test[[varname]]
  y_fc    <- y_fore[, varname]
  y_lower <- ci$lower[, varname]
  y_upper <- ci$upper[, varname]
  
  plot(
    time_obs, y_obs,
    type = "l",
    lwd  = 2,
    xlim = c(min(time_obs), max(time_fc)),
    ylim = range(c(y_obs, y_real, y_lower, y_upper)),
    xlab = "Year",
    ylab = vartext,
    main = paste("TV-VECM", h , "year forecast for", vartext)
    
  )
  
  # Forecast and CI (forecast stays dashed)
  lines(time_fc, y_fc, lwd = 2, lty = 2)
  lines(time_fc, y_lower, lty = 3)
  lines(time_fc, y_upper, lty = 3)
  
  # Realized values (solid line)
  lines(time_fc, y_real, lwd = 2, lty = 1)
  
  # Forecast origin
  abline(v = max(time_obs), lty = 3)
  
  legend(
    "topleft",
    legend = c("Observed", "Forecast", "95% CI"),
    lty    = c(1, 2, 3),
    lwd    = c(2, 2, 1),
    bty    = "n"
  )
}

plot_tvvecm_forecast_oos(
  y_train = y_l,
  y_test  = y_test,
  y_fore  = y_pred_l,
  ci      = ci_l,
  varname = "y1",
  vartext = "GMST",
  years_train = merged_data_lille$year,
  years_test  = merged_data_test$year
)
