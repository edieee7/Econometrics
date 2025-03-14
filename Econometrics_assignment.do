*1. Describe the data 
use "assignment_data_group_24.dta", clear
des


*2. Estimate an employment equation through an OLS regression with dependent variable n and regressors w and k. 

reg n w k

*3. Test the significance of coefficients on w and k and comment. 

test w
test k



*4. By using ys and/or sectoral indicators test the null of homogeneity across sectors 

regress n w k i.sector
testparm i.sector

*5. If you find evidence of sectoral heterogeneity, find a convenient way to accommodate it in your regression model. 



*6. Test the null of conditional homoskedasticity and comment. 

estat hettest
estat imtest, white


*7. If you find evidence of heteroskedasticity, find a way to address the issue and rerun regression analysis with the correction.

regress n w k, robust
