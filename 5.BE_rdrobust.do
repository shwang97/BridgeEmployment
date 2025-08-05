/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the analysis using local polynomial regression
*/

clear
set more off
macro define DTA "LLE"

use $DTA/BE_clean, replace
keep if year>=2012 & year<=2021

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {	

 forvalues i = 1/500 {
 preserve 
 bsample
        * Compute RD estimates for emp in both age groups
        rdrobust emp year_c if dum==1 & `x'!=. , deriv(1) h(6) ///
            p(1) q(2) vce(cluster year_c) kernel(tri)
        scalar coef1 = e(tau_cl)
        rdrobust emp year_c if dum==0 & `x'!=., deriv(1) h(6)  ///
            p(1) q(2) vce(cluster year_c) kernel(tri)
        scalar coef2 = e(tau_cl)
        scalar coef_first = coef1 - coef2

        * Compute RD estimates for the outcome variable in both age groups
        rdrobust `x' year_c if dum==1 & `x'!=., deriv(1) h(6)  ///
            p(1) q(2) vce(cluster year_c) kernel(tri)
        scalar coef3 = e(tau_cl)
	
        rdrobust `x' year_c if dum==0 & `x'!=., deriv(1) h(6) ///
            p(1) q(2) vce(cluster year_c) kernel(tri)
        scalar coef4 = e(tau_cl)
        scalar coef_reduced = coef3 - coef4

        * Calculate the ratio coefficient (the effect on the outcome relative to emp)
        scalar result_`x'_`i' = coef_reduced / coef_first
restore 
  }
}

clear
set obs 500 
  
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {	

  gen result_`x'=. 
  forval i = 1/500 {
  	  replace result_`x'=result_`x'_`i' if _n ==`i'
  }
}

save $DTA/BE_boot, replace


******* plot ********
use $DTA/BE_boot, replace

local i =0
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln Triglyceride_ln chol_HDL_ln  {	
local i = `i'+1 
sum result_`x', detail
scalar coef_rd_`i' = r(mean)
scalar lower_rd_`i' = r(p5)
scalar upper_rd_`i' = r(p95)

gen coef_rd_`i' = coef_rd_`i'
gen lower_rd_`i' = lower_rd_`i'
gen upper_rd_`i' = upper_rd_`i'
}

collapse (mean) coef_rd_* lower_rd_* upper_rd_*

gen id =1 
reshape long coef_rd_ lower_rd_ upper_rd_, i(id)
save $DTA/rdestimates, replace

use $DTA/BE_main_1, replace
merge 1:1 _j using $DTA/rdestimates


forval x = 1/2{
gen upper_`x' = coef`x'_ +1.96*se`x'_ 
gen lower_`x' = coef`x'_ -1.96*se`x'_ 	
}

cap drop _j1 
cap drop _j2 

replace _j = _j+1 
replace _j = _j+1 if _j>=9 

gen _j1 = _j-0.2 
gen _j2 = _j+0.2 


tw (scatter  coef1_ _j1, color("33 145 130") msize(small)) (rcap upper_1 lower_1 _j1, lcolor("33 145 130")) (scatter  coef_rd_ _j2, color("68 1 84") msize(small)) (rcap upper_rd lower_rd _j2, lcolor("68 1 84")),ylabel(-1(0.5)1, nogrid format(%4.1f)) xsize(16) ysize(12) xlabel(1 "{bf: Clinical}" 2 "Obesity" 3 "Hypertension" 4 "ISH" 5 "Diabetes" 6 "Hypercholesterolemia" 7 "Hypertriglyceridemia" 8 "Hypo-HDL"  9 "{bf:Subclinical}" 10 "BMI" 11 "Mid-BP" 12 "Pulse pressure"  13 "Glucose" 14 "Cholestrol" 15 "Triglyceride" 16 "HDL-cholesterol" , nogrid angle(45) labsize(small)) ytitle("") yline(0) legend(order(1 "Parametric" 3 "Nonparametric") pos(1) col(3) ) ytitle("Effect estimates") plotregion(lcolor(black)) 

graph export $DTA/nonparametric.png, width(1600) height(1200) as(png) replace






