/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the robustness checks

Robustness checks include:
1. different bandwidth 
2. different age range for control group 
5. MTE framework 
local polynomial/romano wolf correction in seperate files
*/

clear
set more off
macro define DTA "LLE"
use $DTA/BE_clean, replace

******* different bandwidth choice *******
preserve
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
forval i = 4/6 { 

ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist if year_c>=-`i'& year_c<=`i'-1
scalar coef_`x'_`i' = _b[emp]
scalar se_`x'_`i'= _se[emp]
}
}

clear 
set obs 1
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
forval i =  4/6 {
	gen coef_`x'_`i' = coef_`x'_`i'
	gen se_`x'_`i' = se_`x'_`i'
}		
}

gen id = 1
reshape long coef_obesity_ se_obesity_ coef_high_pressure_ se_high_pressure_ coef_ISH_ se_ISH_ coef_diabetes_ se_diabetes_  coef_high_chol_ se_high_chol_ coef_high_TG_ se_high_TG_ coef_low_HDL_ se_low_HDL_ coef_BMI_ln_ se_BMI_ln_ coef_Midbp_ln_ se_Midbp_ln_ coef_PP_ln_ se_PP_ln_ coef_glucose_ln_ se_glucose_ln_ coef_cholestrol_ln_ se_cholestrol_ln_ coef_chol_HDL_ln_ se_chol_HDL_ln_ coef_Triglyceride_ln_ se_Triglyceride_ln_, i(id)

save $DTA/bandwidth, replace
restore

use $DTA/bandwidth, replace

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln  {

gen upper_`x'= coef_`x'_ + 1.96 * se_`x'_
gen lower_`x'= coef_`x'_ - 1.96 * se_`x'_

}

//// graph for clinical outcomes
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL {
tw (rcap upper_`x' lower_`x' _j, lcolor(navy)) (scatter coef_`x'_ _j, color(navy) lwidth(thick)) (function y=0, range(3.5 6.6) lpattern(dash) lcolor(black)), ylabel(-1(0.5)1, nogrid) xlabel(3.5 " " 4 "4" 5 "5" 6 "6" 6.5 " ", nogrid notick) xsize(12) ysize(12) legend(off) xtitle("Bandwidth", ) ytitle("Effect estimates", ) title("{bf:`x'}", size(medium)) 

graph save $DTA/`x'_bw.gph, replace

}

graph combine $DTA/obesity_bw.gph  $DTA/high_pressure_bw.gph $DTA/ISH_bw.gph $DTA/diabetes_bw.gph $DTA/high_chol_bw.gph $DTA/high_TG_bw.gph $DTA/low_HDL_bw.gph, rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(32)

//graph use $DTA/bw.gph
graph export $DTA/bw.png, as(png) width(3000) height(3200) replace

//// graph for clinical outcomes
foreach x in BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln  {
tw (rcap upper_`x' lower_`x' _j, lcolor(navy)) (scatter coef_`x'_ _j, color(navy) lwidth(thick)) (function y=0, range(3.5 6.6) lpattern(dash) lcolor(black)), ylabel(-1(0.5)1, nogrid) xlabel(3.5 " " 4 "4" 5 "5" 6 "6" 6.5 " ", nogrid notick) xsize(12) ysize(12) legend(off) xtitle("Bandwidth", ) ytitle("Effect estimates", ) title("{bf:`x'}", size(medium)) 
graph save $DTA/`x'_bw.gph, replace

}

graph combine $DTA/BMI_ln_bw.gph  $DTA/Midbp_ln_bw.gph $DTA/PP_ln_bw.gph $DTA/glucose_ln_bw.gph $DTA/cholestrol_ln_bw.gph $DTA/chol_HDL_ln_bw.gph $DTA/Triglyceride_ln_bw.gph, rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(32)

******* different age range for control group******
preserve
keep if year>=2012 & year<=2021

foreach x in  obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln{

forval i =40/55 {

ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist if age>=`i' , robust
scalar coef_`x'_`i' = _b[emp]
scalar se_`x'_`i'= _se[emp]
}
}

clear 
set obs 1
foreach x in  obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
forval i =  40/55  {
	gen coef_`x'_`i' = coef_`x'_`i'
	gen se_`x'_`i' = se_`x'_`i'
}		
}

gen id = 1
reshape long coef_obesity_ se_obesity_ coef_high_pressure_ se_high_pressure_ coef_ISH_ se_ISH_ coef_diabetes_ se_diabetes_  coef_high_chol_ se_high_chol_ coef_high_TG_ se_high_TG_ coef_low_HDL_ se_low_HDL_ coef_BMI_ln_ se_BMI_ln_ coef_Midbp_ln_ se_Midbp_ln_ coef_PP_ln_ se_PP_ln_ coef_glucose_ln_ se_glucose_ln_ coef_cholestrol_ln_ se_cholestrol_ln_ coef_chol_HDL_ln_ se_chol_HDL_ln_ coef_Triglyceride_ln_ se_Triglyceride_ln_, i(id)

save $DTA/BE_agerange, replace

restore

use $DTA/BE_agerange, replace

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln{

gen upper_`x'= coef_`x'_ + 1.96 * se_`x'_
gen lower_`x'= coef_`x'_ - 1.96 * se_`x'_

}

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL {
tw (rarea upper_`x' lower_`x' _j, color(navy*0.3)) (line coef_`x'_ _j, lcolor(navy) lwidth(thick)) (function y=0, range(40 55) lpattern(dash) lcolor(black)), ylabel(-1(0.5)1, nogrid) xlabel(, nogrid) xsize(12) ysize(12) legend(off) xtitle("Age range", ) ytitle("Effect estimates", ) title("{bf:`x'}", size(medium)) 

graph save $DTA/`x'_age.gph, replace

}

graph combine $DTA/obesity_age.gph  $DTA/high_pressure_age.gph $DTA/ISH_age.gph $DTA/diabetes_age.gph $DTA/high_chol_age.gph $DTA/high_TG_age.gph $DTA/low_HDL_age.gph, rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(32)

//graph use $DTA/bw_age.gph
graph export $DTA/bw_age.png, as(png) width(3000) height(3200) replace

******* MTE framework *******
** Note: MTE cannot be assessed using 2SLS, as it examine whether the effect differs by first stage residual

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL {
preserve 
keep if year>=2012 & year<=2021

reg emp i.cut#c.year_c#i.dum i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust
cap drop resid
predict resid, resid 

reg `x' emp resid c.resid#c.emp i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust

matrix b = e(b)
scalar b_emp = b[1, "emp"]
scalar b_int = b[1, "c.resid#c.emp"]

matrix V = e(V)
scalar v_emp      = V["emp", "emp"]
scalar v_int      = V["c.resid#c.emp", "c.resid#c.emp"]
scalar cov_emp_int= V["emp", "c.resid#c.emp"]

**** Create a grid for the latent resistance U_d from 0.05 to 0.95
**** Compute the empirical CDF for resid to form U_d.
summarize resid, detail
scalar u_low = r(p5)
scalar u_high = r(p95)

**** Create a grid for U_d based on its empirical 5th and 95th percentiles
clear
set obs 100
gen ud_grid = u_low + (_n - 1)*(u_high - u_low)/99
gen MTE = b_emp + b_int * ud_grid
gen se_MTE = sqrt(v_emp + (ud_grid^2)*v_int + 2*ud_grid*cov_emp_int)


gen lower = MTE - 1.96*se_MTE
gen upper = MTE + 1.96*se_MTE

cap drop _j
gen _j = _n/100

twoway ///
    (rarea upper lower _j, sort color(navy*0.3)) ///
    (line MTE _j, sort lcolor(navy) lwidth(thick)) ///
    (function y=0, range(0 1) lcolor(black) lpattern(dash)), ///
    xtitle("U{sub:d}") ytitle("MTE") ///
    xlabel(0(0.2)1, nogrid) legend(off) ylabel(-1(0.5)1, nogrid) xsize(12) ysize(12) title("{bf:`x'}", size(medium))

graph save `x'_MTE.gph, replace
restore
}

graph combine obesity_MTE.gph  high_pressure_MTE.gph ISH_MTE.gph diabetes_MTE.gph high_chol_MTE.gph high_TG_MTE.gph low_HDL_MTE.gph, rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(32)

graph use $DTA/MTE.gph
graph export $DTA/BE_MTE.png, as(png) width(3000) height(3200) replace



