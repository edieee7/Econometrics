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
library(forecast)
library(tseries)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(car)

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

##################################################################################
# 3. Diagnostic check
# Purpose : we're testing whether our model satisfies the key assumptions:
# Linearity, Homoskedasticity (constant variance), No autocorrelation in residuals
# Normality of errors, No influential outliers
##################################################################################
# 1. Residual checks
par(mfrow = c(2, 2))
plot(model2)

# 2. Autocorrelation of residuals
acf(residuals(model2))

# 3.Normality of residuals
hist(residuals(model2), main = "Histogram of Residuals", xlab = "Residuals")
shapiro.test(residuals(model2))

# 4. check for multicolinearity
library(car)
vif(model2) 
# shows the Variance Inflation Factor (VIF) 
# less than 5 means there are no serious miscollinearity problem