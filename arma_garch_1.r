####################################
# ARMA-GARCH Workflow
# Done for the 20/03 forecast
# Goal is to predict the spread between BTP and BUND 
# with an ARMA-GARCH model, done via an "automated" function
####################################

#### Packaged and random stuff ####

rm(list = ls()) # clean enviroment
gc()            # same as above
start_time <- proc.time() # to see for how much the code run

library(tseries)
library(quantmod)
library(forecast)
library(rugarch)
library(FinTS)
library(ggplot2)
library(dplyr)
library(httpgd)
library(R.utils)
library(xts)
library(readxl)

# to see plot on vscode: open the link at the start of the output
# if on Rstudio eliminate the line (creates problems??)
hgd() 

#### Load and prepare data ####
#data  <- read.csv("data/spread_investing.csv")
#data  <- xts(data[, -1], order.by = as.Date(data$date))
#head(data)

data <- read_excel("data/spread_refinitiv.xlsx")
data <- data[-1,]
data$`IT10DE10=TWEB (BID)` <- as.numeric(data$`IT10DE10=TWEB (BID)`) * -1

names(data) <- c("date", "spread")
data  <- xts(data[, -1], order.by = as.Date(data$date))
head(data)

spread <- na.omit(data$spread)
ds  <- na.omit(diff(spread))

# Split data into train/test
n_test  <- round(length(ds) * 0.2)
n_train <- length(ds) - n_test

train <- ds[1:n_train]
test  <- ds[(n_train + 1):length(ds)]

p_spread <- ggplot(data, aes(x = index(data), y = spread)) +
            geom_line(linewidth = 0.8) +
            labs(title = "10 Years BTP-BUND Weekly Spread", x = "", y = "Basis points") +
            scale_x_date(date_breaks = "1 year", date_label = "%Y") +
            theme_minimal()
            
print(p_spread)

cat(sprintf("Training observations: %d\n", n_train))
cat(sprintf("Test observations:     %d\n", n_test))
cat(sprintf("Train period: %s -> %s\n", index(train)[1], index(train)[length(train)]))
cat(sprintf("Test period:  %s -> %s\n", index(test)[1],  index(test)[length(test)]))

#### Volatility analysis ####
par(mfrow = c(4, 1), mar = c(3, 4, 3, 2))

# ds
plot(index(ds), coredata(ds), type = "l", lwd = 0.7, 
     main = "Differenced spread (ds)", xlab = "", ylab = "ds")
abline(h = 0, col = "black", lwd = 0.5)

# ACF
Acf(train, lag.max = 30, main = "ACF of ds")

# PACF 
Pacf(train, lag.max = 30, main = "PACF of ds")

# ACF of squared ds (volatility clustering)
Acf(train^2, lag.max = 30, main = "ACF of Squared ds")

par(mfrow = c(1, 1))

#### Model 1: ARMA ####
# Check if ARMA model is needed
cat("\n============================================================\n")
cat("\nCheck if mean modelling is needed:\n")

# Stationarity test
adf_test  <- adf.test(train)
cat("\nADF Test on ds p-value:\n", adf_test$p.value, "\n")
if(adf_test$p.value < 0.05) cat("Stationarity detected\n") else cat("Non-Stationarity detected\n")

# Autocorrelation test
lb_ret <- Box.test(train, lag = 10, type = "Ljung-Box")
cat("\nLB Test on ds p-value:\n", lb_ret$p.value, "\n")
if(lb_ret$p.value < 0.05) cat("Autocorrelation detected, need mean modelling\n") else cat("ds are White Noise\n")

# White noise check
lb_test <- Box.test(train^2, lag = 10, type = "Ljung-Box")
cat(sprintf("\nLjung-Box on Squared ds at lag 10:\n %.4f\n", lb_test$p.value))
if(lb_test$p.value < 0.05) cat("Volatility clustering confirmed\n") else cat("No volatility clustering detected\n")

# ARCH test
arch_test <- ArchTest(train, lags = 5)
cat(sprintf("\nARCH Test (5 lags) LM Test - statistic: %.4f, p-value: %.4f\n", arch_test$statistic, arch_test$p.value))
if(arch_test$p.value < 0.05) cat("ARCH effects confirmed. GARCH modelling needed\n") else cat("No ARCH effects confirmed\n")

# Best ARIMA model
arima_model <- auto.arima(train, 
                          start.p = 1, start.q = 1, 
                          max.p = 5, max.q = 5, 
                          seasonal = FALSE, 
                          stepwise = FALSE, trace = FALSE, method = "ML")

cat("\n============================================================\n")
cat("\nBest ARIMA model:\n")
print(summary(arima_model))
cat("\n")

# ARIMA assumptions
cat("\n============================================================\n")

arma_assumptions <- function(model) {
  cat("\nARIMA model assumptions:\n")
  
  # Check if outside the unit root
  ar_roots <- polyroot(c(1, -model$model$phi))
  ma_roots <- polyroot(c(1, model$model$theta))
  
  is_stationary <- all(abs(ar_roots) > 1)
  is_invertible <- all(abs(ma_roots) > 1)
  
  cat("Stationary:", is_stationary, "\n")
  cat("Invertible:", is_invertible, "\n")
  
  # Check for common roots
  common_roots <- FALSE
  if(length(ar_roots) > 0 && length(ma_roots) > 0) {
    for(ar in ar_roots) {
      for(ma in ma_roots) {
        if(isTRUE(all.equal(ar, ma, tolerance = 1e-3))) common_roots <- TRUE
      }
    }
  }
  
  if(common_roots) {
    cat("Common roots detected, assumption violation\n")
  } else {
    cat("No common roots detected, assumption validated\n")
  }
}

print(arma_assumptions(arima_model))

#### Residual's check ####
cat("\n============================================================\n")
cat("\nARMA residuals check:\n")
# ARIMA residuals
arima_res <- residuals(arima_model)

arma_residual <- function(res) {

  lb_arima <- Box.test(res, lag = 10, type = "Ljung-Box")
  cat("\nLB test on ARIMA residuals (10 lags) p-value: ", round(lb_arima$p.value, 4), "\n")
  if(lb_arima$p.value >= 0.05) cat("ARIMA model cleared the mean\n") else cat("ARIMA model DID NOT clear the mean\n")
  
  lb_arima_2 <- Box.test(res^2, lag = 10, type = "Ljung-Box")
  cat("\nLB test on ARIMA squared residuals (10 lags) p-value: ", round(lb_arima_2$p.value, 4), "\n")
  if(lb_arima_2$p.value < 0.05) cat("ARIMA model has residual variance to explain. GARCH model needed\n") else cat("ARIMA model residual variance is OK\n")
  
  # ARCH test
  arch_test_res <- ArchTest(res, lags = 5)
  cat(sprintf("\nARCH Test on ARIMA (5 lags) residuals:\n LM Test - statistic: %.4f, p-value: %.4f\n", 
              arch_test_res$statistic, arch_test_res$p.value))
  if(arch_test_res$p.value < 0.05) cat("ARCH effects confirmed. GARCH modelling needed\n") else cat("No ARCH effects confirmed\n")
}

print(arma_residual(arima_res))

#### Model 2: ARMA-GARCH ####
cat("\n============================================================\n")
cat("\nFind best ARMA-GARCH model\n")

# Create a function that simulates 1 iteration of a for cycle, then instead of a for cycle use parrallel (mclapply)
# It will run the function over n times using n CPU simultanealy
# Can change the critera or modify the parameters, if too much combinations the code will run for ages

best_garch <- function(y,
                       ar_lags  = 0:2,
                       ma_lags  = 0:2,
                       garch_p  = 1:2,
                       garch_q  = 1:2,
                       models   = c("sGARCH", "eGARCH"),
                       dists    = c("norm", "std", "sstd"),
                       criteria = "aic", 
                       timeout  = 20) {

  library(parallel) # To use all CPU-1 and run more iterations simultanealy

  # Create every combinations possible for the parameter
  combos <- expand.grid(
    ar_lag  = ar_lags,
    ma_lag  = ma_lags,
    garch_p = garch_p,
    garch_q = garch_q,
    model   = models,
    dist    = dists,
    stringsAsFactors = FALSE
  )

  cat("Total combinations to fit:", nrow(combos), "\n")
  cat("Running on", parallel::detectCores() - 1, "cores\n")
  cat("It could take couple of minutes ...\n\n")

  fit_one <- function(i) {
    arp   <- combos$ar_lag[i]
    map   <- combos$ma_lag[i]
    gp    <- combos$garch_p[i]
    gq    <- combos$garch_q[i]
    model <- combos$model[i]
    dist  <- combos$dist[i]

    spec <- ugarchspec(
      variance.model     = list(model = model, garchOrder = c(gp, gq)),
      mean.model         = list(armaOrder = c(arp, map), include.mean = TRUE),
      distribution.model = dist
    )
    # If R encounter an error it stop the cycle, tryCatch catch the error and return NULL, then continue
    fit <- tryCatch(
      withTimeout(
        ugarchfit(spec = spec, data = y, solver = "hybrid", solver.control = list(trace = 0)),
        timeout   = timeout,
        onTimeout = "silent"
      ),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    # Store only not null + right class
    if (!is.null(fit) && inherits(fit, "uGARCHfit")) { 
      ic <- infocriteria(fit)
      return(data.frame(
        ar_p    = arp,
        ma_q    = map,
        garch_p = gp,
        garch_q = gq,
        model   = model,
        dist    = dist,
        aic     = as.numeric(round(ic["Akaike", ] * length(y), 2)), # the function calculate the aic only for that obs
        bic     = as.numeric(round(ic["Bayes",  ] * length(y), 2)),
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  }

# Run more itrations simultanely with every CPU-1, instead of for cycle
  results_list <- mclapply(
    seq_len(nrow(combos)),
    fit_one,
    mc.cores = detectCores() - 1
  )
# Keep only non NULL model
  successful <- Filter(Negate(is.null), results_list)
  cat("Successful fits:", length(successful), "/", nrow(combos), "\n")
  cat("Models ordered by: ", criteria, "\n")

  results_df <- bind_rows(successful) |> arrange(.data[[criteria]])
  return(results_df)
}

results = best_garch(train)
cat("\nFitted models:\n")
print(head(results))

arma_best_order <- c(results$ar_p[1], results$ma_q[1])
garch_best_order  <- c(results$garch_p[1], results$garch_q[1])
best_dist  <- results$dist[1]
best_mod  <- results$model[1]

cat("\nBest ARMA-GARCH model orders:\nARMA:  ", arma_best_order, "\nGARCH: ", garch_best_order, "\n") 

# Fit best model 
garch_spec <- ugarchspec(variance.model = list(model = best_mod, garchOrder = garch_best_order),
                         mean.model = list(armaOrder = arma_best_order, include.mean = TRUE),
                         distribution.model = best_dist)

garch_fit <- ugarchfit(spec = garch_spec, data = train, solver = "hybrid")

cat("\nBest ARMA-GARCH model: \n")
show(garch_fit)

#### Residuals ####
cat("\n============================================================\n")
cat("\nARMA-GARCH residual check:\n")

# The ouput of garch_fit also has LB, ARCH  and APGoF tests on WEIGHTED std residuals
# output <- capture.output(show(garch_fit))
# cat(output[53:76], sep = "\n")
# The Adjusted Pearson Goodnes of Fit test tests whether the standardised residuals 
# follow the fitted distribution across different bin groupings. If YES p>0.05

# After removing both mean (ARIMA) and variance (GARCH) effects the residuals should be iid.

garch_residual <- function(model_fit) {

  # Get standardized residuals
  std_res <- as.numeric(residuals(model_fit, standardize = TRUE))
  
  lb_ret <- Box.test(std_res, lag = 10, type = "Ljung-Box")
  lb_sq <- Box.test(std_res^2, lag = 10, type = "Ljung-Box")
  arch_t <- ArchTest(std_res, lags = 5)
  
  cat(sprintf("\nLjung-Box on GARCH std resid (lag 10) p-value: %.4f (want > 0.05)\n", lb_ret$p.value))
  cat(sprintf("\nLjung-Box on GARCH std resid² (lag 10) p-value: %.4f (want > 0.05)\n", lb_sq$p.value))
  cat(sprintf("\nARCH LM on std resid p-value: %.4f (want > 0.05)\n", arch_t$p.value))
}

print(garch_residual(garch_fit))

#### Validation ####
cat("\n============================================================\n")
cat("\nModel validation:\n")

# Forecast on test data using rolling forecast (window increases on 1 week every iteration). 
# Refit the model every 1 iteration (1 week) to keep model useful
# Perform the rolling estimation and forecast
# refit.every = 1 means we re-estimate the model every single iteration
# If window.size too small the model will fail to converge

c1 <- makeCluster(detectCores() - 1)

roll <- ugarchroll(spec = garch_spec, 
                   data = ds, 
                   n.ahead = 1, 
                   forecast.length = n_test, 
                   refit.every = 1, 
                   refit.window = "moving", 
                   window.size = 300, 
                   solver = "hybrid",
                   calculate.VaR = TRUE,
                   VaR.alpha = c(0.05, 0.25, 0.75, 0.95), 
                   cluster = c1)

# In case of non-convergence problems, re-run with other solvers
roll <- resume(roll, solver = "gosolnp")
roll <- resume(roll, solver = "nloptr")
roll <- resume(roll, solver = "solnp")

# Extract the density predictions and realized ds
roll_df <- as.data.frame(roll)

# Extract the theoretical quantiles
quants <- as.data.frame(roll, which = "VaR")

# Build the forecast_df for better reading
forecast_df <- data.frame(
  pred_mu  = roll_df$Mu,
  pred_var = roll_df$Sigma^2,
  pred_vol = roll_df$Sigma,
  ci_q05  = quants[, "alpha(5%)"],
  ci_q25  = quants[, "alpha(25%)"],
  ci_q75  = quants[, "alpha(75%)"],
  ci_q95  = quants[, "alpha(95%)"],
  actual   = roll_df$Realized
)
rownames(forecast_df) <- rownames(roll_df)

cat("\nForecasting results:\n")
print(tail(round(forecast_df, 4)))

cat("\n")

### Mean evaluation ###
errors <- forecast_df$actual - forecast_df$pred_mu
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
mape <- mean(abs(errors / forecast_df$actual))
dir_acc <- mean(sign(forecast_df$actual) == sign(forecast_df$pred_mu))

cat("\nMean evaluation metrics:\n")
cat(sprintf("MAE: %.4f  |  MSE: %.4f  |  RMSE: %.4f  |  MAPE: %.4f  |  Dir. Acc: %.2f%%\n", 
            mae, mse, rmse, mape, dir_acc * 100))

### Volatility evaluation ###
# Use ds^2 as a proxy for realized volatility
proxy_var <- forecast_df$pred_mu^2
var_mae <- mean(abs(proxy_var - forecast_df$pred_var))
var_mse <- mean((proxy_var - forecast_df$pred_var)^2)
qlike <- mean(log(forecast_df$pred_var) + proxy_var / forecast_df$pred_var)

cat("\nVolatility evaluation metrics:\n")
cat(sprintf("MAE: %.4f  |  MSE: %.4f  |  QLIKE: %.4f\n", var_mae, var_mse, qlike))

### Distribution check ###
# Check if actual value fall into 90%
# 90% interval: q05 to q95
coverage_90 <- mean(forecast_df$actual >= forecast_df$ci_q05 & 
                    forecast_df$actual <= forecast_df$ci_q95)

# 50% interval: q25 to q75
coverage_50 <- mean(forecast_df$actual >= forecast_df$ci_q25 & 
                    forecast_df$actual <= forecast_df$ci_q75)

cat("\nDistribution evaluation, % of forecast that fall into that interval:\n")
cat(sprintf("90%% coverage: %.4f\n", coverage_90)) # should be ~0.90
cat(sprintf("50%% coverage: %.4f\n", coverage_50)) # should be ~0.50

forecast_df$date <- as.Date(rownames(forecast_df))

p_ds  <- ggplot(forecast_df, aes(x = date)) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95, fill = "90% Confidence interval"), alpha = 0.3) +
  geom_line(aes(y = actual, color = "Actual"), linewidth = 0.8) +
  geom_line(aes(y = pred_mu, color = "Forecasted"), linewidth = 0.8) +
  scale_color_manual(values = c("Actual" = "black", "Forecasted" = "blue")) +
  scale_fill_manual(values = c("90% Confidence interval" = "blue")) +
  labs(title = "DS: Actual vs Forecasted", x = "", y = "DS") +
  scale_x_date(date_breaks = "1 year", date_label = "%Y") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom")
print(p_ds)

#### Spread ####
cat("\n============================================================\n")
cat("\nConvert back to spread:\n")

get_spreads <- function() {
  fc_spread <- numeric(n_test)
  lower_spread <- numeric(n_test)
  upper_spread <- numeric(n_test)
  
  for (i in 1:n_test) {
    # precedent spread value
    prev_spread <- as.numeric(spread[n_train + i - 1]) 
    
    fc <- forecast_df$pred_mu[i] 
    fc_spread[i] <- prev_spread + fc
    
    ci_lower <- forecast_df$ci_q05[i] 
    ci_upper <- forecast_df$ci_q95[i] 
    
    lower_spread[i] <- prev_spread + ci_lower
    upper_spread[i] <- prev_spread + ci_upper
  }
  
  actual_spread <- as.numeric(spread[(n_train + 1):(n_train + n_test)])
  
  df <- data.frame(
    actual_spread = actual_spread,
    fc_spread = fc_spread,
    ci_q05 = lower_spread,
    ci_q95 = upper_spread
  )
  rownames(df) <- rownames(forecast_df)
  return(df)
}

df_spread <- get_spreads()
cat("\n")
cat("Forecast spread vs actual spread:\n")
print(tail(df_spread))

## Mean evaluation ###
errors <- df_spread$actual_spread - df_spread$fc_spread
mae <- mean(abs(errors))
mse <- mean(errors^2)
rmse <- sqrt(mse)
mape <- mean(abs(errors / df_spread$actual_spread))
dir_acc <- mean(sign(df_spread$actual_spread) == sign(df_spread$fc_spread))

cat("\nMean evaluation metrics:\n")
cat(sprintf("MAE: %.4f  |  MSE: %.4f  |  RMSE: %.4f  |  MAPE: %.4f  |  Dir. Acc: %.2f%%\n", 
            mae, mse, rmse, mape, dir_acc * 100))

## Distribution check ##
# Check if actual value fall into 90%
# 90% interval: q05 to q95
coverage_90 <- mean(df_spread$actual_spread >= df_spread$ci_q05 & 
                    df_spread$actual_spread <= df_spread$ci_q95)

cat("\nDistribution evaluation, % of forecast that fall into that interval:\n")
cat(sprintf("90%% coverage: %.4f\n", coverage_90)) 

df_spread$date <- as.Date(rownames(df_spread))

p_spreads  <- ggplot(df_spread, aes(x = date)) +
  geom_line(aes(y = actual_spread, color = "Actual"), linewidth = 0.8) +
  geom_line(aes(y = fc_spread, color = "Forecasted"), linewidth = 0.8) +
  geom_ribbon(aes(ymin = ci_q05, ymax = ci_q95), fill = "blue", alpha = 0.3) +
  scale_color_manual(values = c("Actual" = "black", "Forecasted" = "blue")) +
  labs(title = "Spread: Actual vs Forecasted", x = "", y = "Basis points") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = scales::breaks_width(20)) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom")

print(p_spreads)


#### Future forecast ####
cat("\n============================================================\n")
cat("\nForecasting future value:\n")

# Forecast future value 
data_fit = tail(ds, n_test) 
last_date  <- last(index(data_fit))

future_fit <- ugarchfit(spec = garch_spec, data = data_fit, solver = "hybrid")

future_fc  <- ugarchforecast(future_fit, n.ahead = 1)

mean_fc  <-  future_fc  |> fitted()  |> as.numeric()
sigma_fc  <-  future_fc |> sigma()  |> as.numeric()
dist_fc  <- future_fit@model$modeldesc$distribution
df  <- as.numeric(coef(future_fit)['shape'])
probs <- c(0.05, 0.25, 0.75, 0.95)

# take normal distribution and shift to the right one (*)
quantiles <- mean_fc + sigma_fc * qdist(dist_fc, p = probs, mu = 0, sigma = 1, shape = df)

future_fc <- data.frame(
  ds = mean_fc,
  vol = sigma_fc,
  q_05 = quantiles[1],
  q_95 = quantiles[4]
)

spread_conversion <- function(value, q_05, q_95, data){
	last_spread <- as.numeric(tail(data, 1))
	last_date <- data |> index() |> last()
	spread <- last_spread + value
	q_05 <- last_spread + q_05
	q_95 <- last_spread + q_95
	
	df <- data.frame(
	 spread = spread,
	 q_05 = q_05,
	 q_95 = q_95)
	
}

future_spread <- spread_conversion(future_fc$ds, future_fc$q_05, future_fc$q_95, spread)

cat("\nForecasted values:\n")

cat("\nds on last date:", as.character(last_date), ": ", tail(ds, 1))
cat("\nds forecasted for next week:", as.character(last_date + 7), ":\n")
print(future_fc)

cat("\nSpread on last date: ", as.character(last_date), ": ", tail(spread, 1))
cat("\nSpread forecasted for next week:", as.character(last_date + 7), ":\n")
print(future_spread)

#### Simulation forecast ####
cat("\n============================================================\n")
cat("\nSimulate different path for next value:\n")

# Point forecast doesn't think about the shock (that i modelled with GARCH)
# With forecasting I basically only used the ARMA part
# Now I simulate lot of path then take the mean, that will be my forecsat 

simulation <- ugarchsim(
  future_fit,                        
  n.sim   = 1,              # horizon
  m.sim   = 10000,           # n. simulations (paths)     
  startMethod = "sample"    # sample=last value, unconditional= mean
)

mean <- simulation |> fitted() |> mean()
vol <- simulation |> sigma() |> mean()
q_90 <- quantile(fitted(simulation), probs = c(0.05, 0.95))

cat("\nTry to forecast value simulating 10.000 possible paths:\n")

cat("\nds simulated over 10.000 paths for next week", as.character(last_date + 7), ":\n")
cat("\nMean DLA: ", mean)
cat("\nMean Vol:", vol)
cat("\nQuantiles 5% | 95%: ", q_90, "\n")

# Conversion to spread
sim_spread <- spread_conversion(value = mean, q_05 = unname(q_90[1]), q_95 = unname(q_90[2]), data = spread)
cat("\nSpread simulated over 10.000 paths for next week", as.character(last_date + 7), ":\n")
print(sim_spread)


# Print execution time
end_time <- proc.time()
execution_time <- end_time - start_time
cat("\n Total execution time: ", execution_time[3], "seconds")
