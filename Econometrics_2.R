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

# Load Ireland macro data
setwd("/Users/edy/Downloads/Bocconi_ESS_1y_2nd/Econometrics/Part2/Assignment")
data <- read_csv("ireland_data.csv")
View(data)
# Filtering out NA, relevant years with complete inflation + regressors
df <- subset(data,  !is.na(infl) & !is.na(unemp) & !is.na(strate))
View(df)
#########################################
# 1. Descriptive Analysis
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

#########################################
# 2. Model Specification
#########################################

#-------------------------------
# STEP 1: Check for stationarity
#-------------------------------

# Visual inspection of inflation time series
plot(df$infl, type = "l", main = "Inflation Time Series", ylab = "Inflation Rate")

# ADF Test (Augmented Dickey-Fuller)
adf_result <- adf.test(df$infl)
print(adf_result)
# If p-value < 0.05 → Stationary; if > 0.05 → Non-stationary

#-------------------------------
# STEP 2: ACF and PACF for lag structure
#-------------------------------

# Plot ACF and PACF to guide AR lag selection
acf(df$infl, main = "ACF of Inflation")
pacf(df$infl, main = "PACF of Inflation")
# PACF cuts off at lag 2? Suggests AR(2)

#-------------------------------
# STEP 3: AR model comparisons using AIC
#-------------------------------

# Try AR(1), AR(2), AR(3)
aic1 <- AIC(arima(df$infl, order = c(1, 0, 0)))
aic2 <- AIC(arima(df$infl, order = c(2, 0, 0)))
aic3 <- AIC(arima(df$infl, order = c(3, 0, 0)))

cat("AIC for AR(1):", aic1, "\n")
cat("AIC for AR(2):", aic2, "\n")
cat("AIC for AR(3):", aic3, "\n")
# Lowest AIC = preferred model

# Fit selected AR(2) model
model_ar2 <- arima(df$infl, order = c(2, 0, 0))
summary(model_ar2)

#-------------------------------
# STEP 4: Include Explanatory Variables (Static and Dynamic)
#-------------------------------

# Model 1: Static Linear Regression (OLS)
model1 <- lm(infl ~ unemp + strate, data = df)
summary(model1)

# Model 2: Dynamic with 1 lag of all variables
df_lag <- df %>%
  mutate(
    infl_lag1 = lag(infl, 1),
    unemp_lag1 = lag(unemp, 1),
    strate_lag1 = lag(strate, 1)
  )

model2 <- lm(infl ~ infl_lag1 + unemp_lag1 + strate_lag1, data = df_lag)
summary(model2)

#-------------------------------
# STEP 5: Try ADL (Auto Distributed Lag) model
#-------------------------------

df_adl <- df %>%
  mutate(
    infl_lag1 = lag(infl, 1),
    unemp_lag1 = lag(unemp, 1),
    strate_lag1 = lag(strate, 1)
  )

model_adl <- lm(infl ~ infl_lag1 + unemp + unemp_lag1 + strate + strate_lag1, data = df_adl)
summary(model_adl)

#-------------------------------
# STEP 6: First Differences (optional, if inflation is non-stationary)
#-------------------------------

df <- df %>%
  mutate(diff_infl = c(NA, diff(infl)))

model_diff <- lm(diff_infl ~ lag(diff_infl, 1), data = df)
summary(model_diff)

#-------------------------------
# STEP 7: Residual Diagnostics
#-------------------------------

# ACF of residuals for AR(2) model
acf(resid(model_ar2), main = "ACF of Residuals (AR2)")

# Ljung-Box test (for autocorrelation)
Box.test(resid(model_ar2), lag = 10, type = "Ljung-Box")

# Check residuals from regression models
checkresiduals(model_ar2)
checkresiduals(model2)


# Model1 Static Linear Regression (OLS)
model1<- lm(infl ~ unemp + strate, data = df)
summary(model1)

# Model2 Include Lags (Dynamic Model)
# Start fresh: reload your filtered data

df_lag <- df %>%
  # Keep only rows where all needed variables are present
  #filter(year >= 1970 & !is.na(infl) & !is.na(unemp) & !is.na(strate)) %>%
  mutate(
    infl_lag1 = lag(infl, 1),
    unemp_lag1 = lag(unemp, 1),
    strate_lag1 = lag(strate, 1)
  )
summary(df_lag)  

# Now estimate the model
model2 <- lm(infl ~ infl_lag1 + unemp_lag1 + strate_lag1, data = df_lag)
summary(model2)


# Model3 AR(2) model (for comparison later)
model_ar2 <- arima(df$infl, order = c(2, 0, 0))
model_ar2

# ACF and PACF for inflation
acf(df$infl, main = "ACF of Inflation")
pacf(df$infl, main = "PACF of Inflation")

# Compare AR models with different orders
AIC(arima(df$infl, order = c(1,0,0)))
AIC(arima(df$infl, order = c(2,0,0)))
AIC(arima(df$infl, order = c(3,0,0)))

##################################################################################
# 3. Diagnostic check
# Purpose : Test if the model satisfies the key OLS assumptions
# - Linearity
# - Homoskedasticity (constant variance)
# - No autocorrelation in residuals
# - Normality of errors
# - No influential outliers
# - No multicollinearity
##################################################################################

# Assume 'model2' is your preferred model from Question 2
# model2 <- lm(infl ~ infl_lag1 + unemp_lag1 + strate_lag1, data = df_lag)

#-------------------------------
# 1. Residual plots: check linearity, homoskedasticity, outliers
#-------------------------------
par(mfrow = c(2, 2))
plot(model2)  # Top-left: residuals vs fitted (linearity & variance)
              # Top-right: Q-Q plot (normality)
              # Bottom-left: Scale-location (variance)
              # Bottom-right: Cook's distance (influence)

#-------------------------------
# 2. Autocorrelation of residuals
#-------------------------------
acf(residuals(model2), main = "ACF of Residuals (Model 2)")
Box.test(residuals(model2), lag = 10, type = "Ljung-Box")
# p-value > 0.05 suggests no significant autocorrelation

#-------------------------------
# 3. Normality of residuals
#-------------------------------
hist(residuals(model2), main = "Histogram of Residuals", xlab = "Residuals")
shapiro.test(residuals(model2))
# p-value > 0.05 suggests residuals are approximately normal

#-------------------------------
# 4. Homoskedasticity (Breusch-Pagan test)
#-------------------------------
bptest(model2)
# p-value > 0.05 suggests constant variance (good)

#-------------------------------
# 5. Multicollinearity (Variance Inflation Factor)
#-------------------------------
vif(model2)
# VIF < 5 suggests no serious multicollinearity problem

#-------------------------------
# 6. Influential outliers (Cook's Distance)
#-------------------------------
plot(cooks.distance(model2), type = "h", main = "Cook's Distance")
abline(h = 4/length(residuals(model2)), col = "red", lty = 2)
# Points above the red line are influential
# You can print them:
which(cooks.distance(model2) > 4/length(residuals(model2)))
