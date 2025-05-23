#########################################
# Load and prepare data
#########################################
install.packages("readr")
install.packages("forecast", type = "binary")
install.packages("tseries", type = "binary")
install.packages("gridExtra")
install.packages("dplyr")
install.packages("car", type = "binary")

library(readr)
library(forecast)   # for ARIMA, diagnostics
library(tseries)    # for ADF test
library(gridExtra)
library(dplyr)
library(ggplot2)
library(car)        # for vif
library(tidyverse)
library(lmtest)       # for bptest and Box.test

#------------- Load Ireland macro data -----------
data <- read_csv("ireland_data.csv")
View(data)

#------------- Filtering out NA, relevant years with complete inflation + regressors -----------
df <- subset(data,  !is.na(infl) & !is.na(unemp) & !is.na(strate))
df_orig <- df  # save original filtered data
View(df_orig)

#########################################
# Question 1. Descriptive Analysis
#########################################
summary(df[, c("year", "infl", "unemp", "strate")])

# Plot inflation
p1 <- ggplot(df, aes(x = year, y = infl)) + 
  geom_line(color = "darkred") + 
  ggtitle("Inflation Rate in Ireland (YoY)") +
  ylab("Inflation (%)") + xlab("Year")

print(p1)

# Plot unemployment
p2 <- ggplot(df, aes(x = year, y = unemp)) + 
  geom_line(color = "steelblue") +
  ggtitle("Unemployment Rate in Ireland") +
  ylab("%") + xlab("Year")

print(p2)

# Plot interest rate
p3 <- ggplot(df, aes(x = year, y = strate)) + 
  geom_line(color = "darkgreen") +
  ggtitle("Short-Term Interest Rate in Ireland") +
  ylab("%") + xlab("Year")

grid.arrange(p1, p2, p3, ncol = 1)

#------------- Descriptive Trend ---------------

# From the line plots:
# - Inflation in Ireland was highly volatile in the 1970s and 1980s, reaching over 20%.
# - Inflation has stabilized since the 1990s and remained below 5%.
# - Unemployment reached its high around the 1980s and again around the 2008 crisis, before falling consistently.
# - Short-term interest rates were high in the 1980s and have fallen consistently, reaching near zero after 2010.

#########################################
# Question 2. Model Specification
#########################################

#------------- Model 1: Static Linear Regression (OLS) ---------------
model1 <- lm(infl ~ unemp + strate, data = df_orig)
summary(model1)

#------------- Model 2: Dynamic with 1 lag of all variables ----------
df_lag <- df_orig %>%
  mutate(
    infl_lag1 = lag(infl, 1),
    unemp_lag1 = lag(unemp, 1),
    strate_lag1 = lag(strate, 1)
  )

model2 <- lm(infl ~ infl_lag1 + unemp_lag1 + strate_lag1, data = df_lag)
summary(model2)


#------------- Model 3: ADL (Auto Distributed Lag) model -------------

df_adl <- df_orig %>%
  mutate(
    infl_lag1 = lag(infl, 1),
    unemp_lag1 = lag(unemp, 1),
    strate_lag1 = lag(strate, 1)
  )

model_adl <- lm(infl ~ infl_lag1 + unemp + unemp_lag1 + strate + strate_lag1, data = df_adl)
summary(model_adl)

#------------- Model 4: ARIMA model  ----------------------------------
# Later for differencing
df_diff <- df_orig %>% 
  mutate(diff_infl = c(NA, diff(infl)))

#-------------------------------
# STEP 1: Check for stationarity
#-------------------------------

# Visual inspection of inflation time series
plot(df_orig$infl, type = "l", main = "Inflation Time Series", ylab = "Inflation Rate")

# ADF Test (Augmented Dickey-Fuller)
adf_result <- adf.test(df_orig$infl)
print(adf_result)
# If p-value < 0.05 → Stationary; if > 0.05 → Non-stationary

#------------- Results ---------------

# Dickey-Fuller = -2.1381  
# Lag order = 4  
# p-value = 0.5193  
# since p-value is more than 0.05, we fail to reject the null hypothesis
# Inflation based on CPI is non stationary.

#-------------------------------
# STEP 2: ACF and PACF of original inflation
#-------------------------------

# Plot ACF and PACF to guide AR lag selection
acf(df_orig$infl, main = "ACF of Inflation")
pacf(df_orig$infl, main = "PACF of Inflation")

#------------- Interpretation ---------------

# ACF decays very slowly and this implies inflation may follow a random walk or contain a unit root. 
# This relsult is consistent with our ADF test result and supports differencing the series before modeling.
# Furtehrmore, in PACF, we could see the sharp drop after lag 1 or 2, then the bars become very small.
# This is a chcharacteristic of an AR process and the PACF cutting off at lag 2 suggests that an AR(2) model 
# (or higher) might be a good fit after differencing, or we might try an ARIMA(2,1,0).

#-------------------------------
# STEP 3: First-Difference the Inflation Series
#-------------------------------
# ADF test on differenced inflation
adf.test(na.omit(df_diff$diff_infl))

#------------- Results ---------------

# Dickey-Fuller = -4.8859  
# Lag order = 4  
# p-value = 0.01 
# since p-value is less than 0.05, we reject the null hypothesis and conclude that 
# differenced Inflation series is now stationary.

#-------------------------------
# STEP 6: Plot ACF and PACF of Differenced Series
#-------------------------------
acf(na.omit(df_diff$diff_infl), main = "ACF of Differenced Inflation")
pacf(na.omit(df_diff$diff_infl), main = "PACF of Differenced Inflation")

#------------ Interpretation ---------------

# From AIC plotting, only the first lag is significantly different from 0 and it suggests a short memory process.
# ALso, the sharp drop after lag 1 is typical of a Moving Average (MA) process and implies a MA(1) component in the model.
# However, there were no clear strong spikes since most bars are small and inside the significance bands.
# it Could suggest that there’s no strong AR pattern.
# Combined with the ACF, this supports, ARIMA(0,1,1) or possibly ARIMA(1,1,1) as a comparison

#-------------------------------
# STEP 7: ARIMA model comparisons using AIC
#-------------------------------
model_arima011 <- arima(df_diff$infl, order = c(0, 1, 1))  # ARIMA(0,1,1)
model_arima111 <- arima(df_diff$infl, order = c(1, 1, 1))  # ARIMA(1,1,1)

# Compare AIC (BIC is not aligned with our aim since we're seeking a model to forecase inflation )
AIC(model_arima011)
AIC(model_arima111)

#------------ Results ---------------

# AIC (Akaike Information Criterion) value of ARIMA(0,1,1) was 477.2666
# the value of ARIMA(1,1,1) was 473.4808
# Since lower AIC value is considered as better model fit, while penalizing complexity.
# Thus, ARIMA(1,1,1) fits the inflation data better than ARIMA(0,1,1).

##################################################################################
# Question 3. Diagnostic check
# Purpose : Test if the model satisfies the assumptions for time series models
# - Residuals are white noise (no autocorrelation)
# - Homoskedasticity (constant variance)
# - Normality of residuals
##################################################################################

# 'model_arima111' (ARIMA(1,1,1)) is our preferred model from Question 2

#-------------------------------
# Step 1. Residual diagnostics (checkresiduals)
#-------------------------------
checkresiduals(model_arima111)
# Includes: ACF plot, residual time series plot, histogram + Ljung-Box test
# The Ljung-Box test checks for autocorrelation in residuals and its null hypothesis is residuals are white noise.
# then p-value was 0.2748 > 0.05, so we fail to reject null hypothesis and conclude residuals are white noise


#-------------------------------
# Step 2. ACF of residuals (manual check)
#-------------------------------
acf(resid(model_arima111), main = "ACF of Residuals (ARIMA(1,1,1))")
# The bars are inside the blue bands which means residuals have no significant autocorrelation

#-------------------------------
# Step 3. Normality of residuals
#-------------------------------
hist(resid(model_arima111), main = "Histogram of Residuals", xlab = "Residuals")
shapiro.test(resid(model_arima111))


#------------- Interpretation

# The residual histogram appears approximately normal with mild right-skewness and one outlier, 
# consistent with the Shapiro-Wilk test result (p = 0.01064). 
# However, the distribution is close enough to normal for ARIMA forecasting purposes, and does not invalidate the model.



#------------- Final conclusion ---------------

# ARIMA(1,1,1) passes all residual checks:
# - Residuals are not autocorrelated (ACF + Ljung-Box)
# - Residuals are approximately normal (histogram + Shapiro)
# - Model is valid for forecasting inflation


##############

# 5.Construct and evaluate 1-step ahead point forecasts for the last 10 years of the sample.
# 6.Compare your forecasts with those resulting from an AR (2) model.

##############

# From previous analysis we assume that model2 is our preffered model

#####
#Step 1: Create rolling 1-step ahead forecasts (last 10 years) using your preferred model
#####

# Sort data just in case

df_lag <- df_lag[order(df_lag$year), ]

# Remove NAs due to lags
df_lag <- df_lag %>% filter(!is.na(infl_lag1) & !is.na(unemp_lag1) & !is.na(strate_lag1))

# Store forecasted vs actual values
n <- 10
n_total <- nrow(df_lag)
forecasts_model2 <- numeric(n)
actuals <- df_lag$infl[(n_total - n + 1):n_total]

# Rolling forecast: fit model on expanding window
for (i in 1:n) {
  train_end <- n_total - n + i - 1
  train_data <- df_lag[1:train_end, ]
  
  model_temp <- lm(infl ~ infl_lag1 + unemp_lag1 + strate_lag1, data = train_data)
  
  newdata <- df_lag[train_end + 1, ]
  forecasts_model2[i] <- predict(model_temp, newdata = newdata)
}
 

#####
#Step 2: Do the same for the AR(2) model
#####

forecasts_ar2 <- numeric(n)

for (i in 1:n) {
  train_end <- n_total - n + i - 1
  train_data <- df_lag$infl[1:train_end]
  
  model_ar_temp <- arima(train_data, order = c(2, 0, 0))
  forecasts_ar2[i] <- forecast(model_ar_temp, h = 1)$mean
}

####
#Step 3: Compare forecast performance
####

mae_model2 <- mean(abs(forecasts_model2 - actuals))
mae_ar2 <- mean(abs(forecasts_ar2 - actuals))

rmse_model2 <- sqrt(mean((forecasts_model2 - actuals)^2))
rmse_ar2 <- sqrt(mean((forecasts_ar2 - actuals)^2))

cat("Model 2 - MAE:", mae_model2, " RMSE:", rmse_model2, "\n")
cat("AR(2) Model - MAE:", mae_ar2, " RMSE:", rmse_ar2, "\n")

####
#Step 4: Visual comparison
####

years_test <- df_lag$year[(n_total - n + 1):n_total]
df_forecasts <- data.frame(
  Year = years_test,
  Actual = actuals,
  Model2 = forecasts_model2,
  AR2 = forecasts_ar2
)

ggplot(df_forecasts, aes(x = Year)) +
  geom_line(aes(y = Actual, color = "Actual")) +
  geom_line(aes(y = Model2, color = "Model 2 (Dynamic)")) +
  geom_line(aes(y = AR2, color = "AR(2)")) +
  labs(title = "1-Step Ahead Forecasts: Last 10 Years",
       y = "Inflation Rate (%)", color = "Legend") +
  theme_minimal()

#Used for comparing forecasts

mae_model2 <- mean(abs(forecasts_model2 - actuals))
rmse_model2 <- sqrt(mean((forecasts_model2 - actuals)^2))

mae_ar2 <- mean(abs(forecasts_ar2 - actuals))
rmse_ar2 <- sqrt(mean((forecasts_ar2 - actuals)^2))

# Print
cat("Model 2 - MAE:", mae_model2, " RMSE:", rmse_model2, "\n")
cat("AR(2) Model - MAE:", mae_ar2, " RMSE:", rmse_ar2, "\n")

#Comment on model comparison 

#To evaluate forecasting performance, we compared 1-step ahead predictions from two models: a dynamic regression model (model 2) and an AR(2) model. Using data from the last 10 years of the sample, we computed standard forecast accuracy metrics. Model 2 yielded a Mean Absolute Error (MAE) of 1.69 and a Root Mean Squared Error (RMSE) of 2.47, while the AR(2) model had a slightly lower RMSE of 2.31 but a higher MAE of 1.83. This suggests that although the AR(2) model performs slightly better in penalizing large forecast errors (lower RMSE), the dynamic regression model produces more accurate forecasts on average (lower MAE). Depending on the chosen loss function, model 2 may be preferred due to its ability to better capture average forecast behavior, especially given its use of relevant macroeconomic predictors (lagged unemployment and interest rates) that provide additional explanatory power.









#### New code with comments

###Step 1: Identify forecast period

# Keep only complete observations (no NAs)
df_clean <- na.omit(df_orig)

# Identify the last 10 years of data
n <- nrow(df_clean)
forecast_years <- df_clean$year[(n-9):n]
print(forecast_years)

###Step 2: Rolling 1-step ahead forecasts with ARIMA(1,1,1)


  # Ensure df_clean is defined and has enough data
  df_clean <- df_orig  # or use your filtered df
  n <- nrow(df_clean)  # total number of observations
  
  # Sanity check
  if (n < 20) stop("Not enough data for a 10-year rolling forecast.")
  
  # Identify last 10 years
  forecast_years <- df_clean$year[(n - 9):n]
  actual <- df_clean$infl[(n - 9):n]
  
  # Empty vector to store forecasts
  forecast_arima <- numeric(10)
  
  # Rolling 1-step ahead ARIMA(1,1,1) forecast
  for (i in 1:10) {
    train_end <- n - 10 + (i - 1)
    train_data <- df_clean$infl[1:train_end]
    
    model <- arima(train_data, order = c(1, 1, 1))
    forecast_result <- predict(model, n.ahead = 1)
    forecast_arima[i] <- forecast_result$pred
  }
  
  # Results
  results_arima <- data.frame(
    Year = forecast_years,
    Actual = actual,
    Forecast = forecast_arima,
    Error = actual - forecast_arima
  )
  print(results_arima)
  
###Step 3: Perform the same rolling 1-step ahead forecast using the AR(2) model

  # Initialize vector for AR(2) forecasts
  forecast_ar2 <- numeric(10)
  
  # Rolling 1-step ahead AR(2) forecast
  for (i in 1:10) {
    train_end <- n - 10 + (i - 1)
    train_data <- df_clean$infl[1:train_end]
    
    model_ar2 <- arima(train_data, order = c(2, 0, 0))  # AR(2) model
    forecast_result_ar2 <- predict(model_ar2, n.ahead = 1)
    forecast_ar2[i] <- forecast_result_ar2$pred
  }
  
  # Results for AR(2)
  results_ar2 <- data.frame(
    Year = forecast_years,
    Actual = actual,
    Forecast = forecast_ar2,
    Error = actual - forecast_ar2
  )
  
  print(results_ar2)
  
  # Step 4: Compare Forecast Accuracy of ARIMA(1,1,1) and AR(2) Models
  
  # Calculate RMSE and MAE for ARIMA(1,1,1)
  rmse_arima <- sqrt(mean((results_arima$Actual - results_arima$Forecast)^2))
  mae_arima <- mean(abs(results_arima$Actual - results_arima$Forecast))
  
  # Calculate RMSE and MAE for AR(2)
  rmse_ar2 <- sqrt(mean((results_ar2$Actual - results_ar2$Forecast)^2))
  mae_ar2 <- mean(abs(results_ar2$Actual - results_ar2$Forecast))
  
  # Create a summary table
  accuracy_comparison <- data.frame(
    Model = c("ARIMA(1,1,1)", "AR(2)"),
    RMSE = c(rmse_arima, rmse_ar2),
    MAE = c(mae_arima, mae_ar2)
  )
  
  # Print with title
  cat("\nForecast Accuracy Comparison (RMSE and MAE):\n")
  print(accuracy_comparison)
  

  #Just to get information for the comments
  
  print(model_arima111)
  
  summary(model_ar2)

  # -------------------------------------------------------------------------
  # Forecast Construction and Comparison
  # -------------------------------------------------------------------------
  # Using the preferred ARIMA(1,1,1) model, with estimated parameters 
  # ar1 = -0.2671 (s.e. 0.3302) and ma1 = 0.5356 (s.e. 0.2843), we constructed 
  # 1-step ahead point forecasts for the last 10 years of the sample. This rolling 
  # forecast procedure involved re-estimating the model each year using all data 
  # up to that year, ensuring that the model adapted over time. 
  #
  # For comparison, forecasts were also generated from an AR(2) model estimated 
  # on the original inflation series, which had parameter estimates ar1 = 1.1068 
  # (s.e. 0.1211) and ar2 = -0.2731 (s.e. 0.1209), with an intercept of 4.6623.
  #
  # Evaluating forecast accuracy over the last decade:
  # - The AR(2) model achieved a slightly lower RMSE of 2.31 compared to 2.34 for 
  #   the ARIMA(1,1,1) model, suggesting marginally better fit in terms of squared errors.
  # - However, the ARIMA(1,1,1) model outperformed in terms of MAE (1.74 vs. 1.81), 
  #   indicating better median forecast accuracy.
  #
  # The ARIMA model’s differencing step addresses the non-stationarity found in 
  # the inflation series, supporting its theoretical suitability despite the AR(2) 
  # model’s competitive performance.
  #
  # In conclusion, while both models exhibit similar forecasting capabilities for 
  # inflation in Ireland, the ARIMA(1,1,1) model provides a more theoretically sound 
  # framework by accounting for unit roots and integrating autoregressive and moving 
  # average components, justifying its role as the preferred forecasting model.
  # -------------------------------------------------------------------------
  
















