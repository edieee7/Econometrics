*1. Describe the data 
use "assignment_data_group_24.dta", clear
des

tab sector, plot sort 
tabstat ys - n, s(mean median sd skew k range min max) 
tabstat ys - n,  s(mean median sd skew k range min max) by (sector) columns (statistics) 
*--------------------------------------------------------------------------------
*2. Estimate an employment equation through an OLS regression with dependent variable n and regressors w and k. 
reg n w k
*--------------------------------------------------------------------------------
*3. Test the significance of coefficients on w and k and comment. 

test w
test k

*--------------------------------------------------------------------------------
*4. By using ys and/or sectoral indicators test the null of homogeneity across sectors 

regress n w k i.sector
testparm i.sector
*--------------------------------------------------------------------------------
*5. If you find evidence of sectoral heterogeneity, find a convenient way to accommodate it in your regression model. 


*-------------------------------------------------------------------------------
*6. Test the null of conditional homoskedasticity and comment. 

* BP test
estat hettest

* --> Null hypothesis is homoskedasticity
* --> p-value was 0.4159 (lower than 0.05), meaning fail to reject H0, homoskedasticity

* White test
estat imtest, white
*similar to BP test, but it checks if the variance of residuals depends on independent variables in a nonlinear way.

* --> p-value was 0.0825, meaning fail to reject H0, homoskedasticity

* Graphically check
rvfplot, yline(0)

* Graph shows residuals in y-axis and fitted value in x-axis. 
* It can be observed that residuals randomly scattered and it does not show any clustering, patterns, and systematic increasing or decreasing. Variance looks constant. 
* It confirms there are no heteroskedasticity.

*--------------------------------------------------------------------------------
*7. If you find evidence of heteroskedasticity, find a way to address the issue and rerun regression analysis with the correction.

*Since both BP tests and White test indicates there is no heteroskedasticity at the 5% level, so it is not necessary to adjust more. 
*In case of adjusting heteroskedasticity, HC1, HC2, HC3 improve on the standard White heteroskedasticity-consistent standard errors. 

regress n w k i.sector, robust
*which is equal to regress n w k i.sector, vce(hc1)
regress n w k i.sector, vce(hc2)
regress n w k i.sector, vce(hc3)

*As a result, standar error increaed from HC1 to HC3, but p-values were unchanged.
