/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the causal decomposition analysis
*/

clear
set more off
macro define DTA "LLE"
use $DTA/BE_clean, replace
keep if year>=2012 & year<=2021

**** High SES = college or more (including some college/2yearcollege)
**** low SES = HS or less
cap drop college 
gen college =1 if educ>=1&educ<=5
replace college=0 if educ>=6&educ<=8

replace sex = sex-1

cap drop emp_int
gen emp_int = emp * college 
keep if college !=. 

***** first stage relevance ****
ivreg2 obesity (emp_int emp = 1.cut#c.year_c#1.dum 1.cut#c.year_c#1.dum#1.college) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.cut#i.college c.year_c#i.college i.cut#c.year_c#i.college i.cut#i.dum#i.college c.year_c#i.dum#i.college i.dum#i.college i.college i.sex i.sex#i.college, first

***** Causal decomposition analysis******
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln Triglyceride_ln chol_HDL_ln  {	
* ------------ group -> mediator
reg emp college sex if dum==1 & `x'!=., 
est store emp_m           

* ------------ mediator -> outcomes
quietly reg emp i.cut#c.year_c#i.dum i.cut#c.year_c#i.dum#i.college i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.cut#i.college c.year_c#i.college i.cut#c.year_c#i.college i.cut#i.dum#i.college c.year_c#i.dum#i.college i.dum#i.college i.college i.sex i.sex#i.college if `x' !=. 
cap drop emp_hat
predict emp_hat if e(sample)

quietly reg emp_int i.cut#c.year_c#i.dum i.cut#c.year_c#i.dum#i.college i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.cut#i.college c.year_c#i.college i.cut#c.year_c#i.college i.cut#i.dum#i.college c.year_c#i.dum#i.college i.dum#i.college i.college i.sex i.sex#i.college  if `x' !=. 
cap drop emp_int_hat
predict emp_int_hat if e(sample)

reg `x' emp_hat emp_int_hat ///
    i.cut c.year_c i.cut#c.year_c ///
    i.cut#i.dum c.year_c#i.dum i.dum ///
    i.cut#i.college c.year_c#i.college ///
    i.cut#c.year_c#i.college ///
    i.cut#i.dum#i.college c.year_c#i.dum#i.college ///
    i.dum#i.college i.college ///
    i.sex i.sex#i.college if `x'!=., 
est store m_y    

* ------------ "group -> outcomes" 
reg `x' college sex if dum==1 &`x'!=., 
*** total disparity 
scalar coef_`x'_1 = _b[college]
scalar se_`x'_1 = _se[college]
est store emp_y
         
** For SE calculation based on the delta method            		      
suest emp_m m_y emp_y, vce(robust)

*** PIE = diff in prevalence * baseline effect
nlcom PIE: _b[emp_m_mean:college] * _b[m_y_mean:emp_hat]
scalar coef_`x'_2 = r(table)[1, 1]
scalar se_`x'_2 = r(table)[2, 1]

*** INTref = baseline prev * diff in impact
sum sex if dum==1 &`x'!=.
scalar sex_m = r(mean)
nlcom INTref:  (_b[emp_m_mean:_cons]+_b[emp_m_mean:sex]*sex_m)* _b[m_y_mean:emp_int_hat]
scalar coef_`x'_3 = r(table)[1, 1]
scalar se_`x'_3 = r(table)[2, 1]

*** INTmed = diff in prev * diff in impact
nlcom INTmed: _b[emp_m_mean:college] * _b[m_y_mean:emp_int_hat]
scalar coef_`x'_4 = r(table)[1, 1]
scalar se_`x'_4 = r(table)[2, 1]

***** CDE = TE - PE (portion eliminated) = TE - (PIE+INTmed+INTref)
nlcom CDE: _b[emp_y_mean:college] - (_b[emp_m_mean:college] * _b[m_y_mean:emp_hat] + _b[emp_m_mean:_cons]  * _b[m_y_mean:emp_int_hat] +   _b[emp_m_mean:college] * _b[m_y_mean:emp_int_hat])  
scalar coef_`x'_5 = r(table)[1, 1]
scalar se_`x'_5 = r(table)[2, 1]
}

clear 
set obs 1
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
forval i = 1/5 {
	gen coef_`x'_`i' = coef_`x'_`i'
	gen se_`x'_`i' = se_`x'_`i'
}		
}

gen id = 1
reshape long coef_obesity_ se_obesity_ coef_high_pressure_ se_high_pressure_ coef_ISH_ se_ISH_ coef_diabetes_ se_diabetes_  coef_high_chol_ se_high_chol_ coef_high_TG_ se_high_TG_ coef_low_HDL_ se_low_HDL_ coef_BMI_ln_ se_BMI_ln_ coef_Midbp_ln_ se_Midbp_ln_ coef_PP_ln_ se_PP_ln_ coef_glucose_ln_ se_glucose_ln_ coef_cholestrol_ln_ se_cholestrol_ln_ coef_chol_HDL_ln_ se_chol_HDL_ln_ coef_Triglyceride_ln_ se_Triglyceride_ln_, i(id)

save $DTA/BE_CD, replace

***** graph *****
use $DTA/BE_CD, replace
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
capture drop p_`x'
gen test_`x' = coef_`x'_ / se_`x'_
gen p_`x' = 2*(1 - normprob(abs(test_`x')))
}

foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
gen upper_`x' = coef_`x'_ + 1.96 * se_`x'_
gen lower_`x' = coef_`x'_ - 1.96 * se_`x'_
}


foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
gen ml_`x' = "†" if p_`x'<0.10
replace ml_`x' = "*" if p_`x' <0.05
replace ml_`x' = "**" if p_`x' <0.01
replace ml_`x' = "***" if p_`x' <0.001
}

replace _j = _j-1

//// graph for clinical outcomes
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL {
tw (bar coef_`x'_ _j if _j==0, color("68 1 84") mlab(ml_`x') mlabsize(4) mlabcolor(black)mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==1, color("59 82 139") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==2, color("33 145 140") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==3, color("94 201 98") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==4, color("253 231 37") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)), xsize(10) ysize(12) legend(off) xtitle("") ytitle("") xlabel(0 "Total" 1 "PIE" 2 "INTref" 3 "INTmed"  4 "CDE" ,  angle(45) nogrid) ylabel(-0.3(0.1)0.3, nogrid) title("{bf:`x'}", size(medium)) 

graph save `x'_graph.gph, replace
}

graph combine obesity_graph.gph  high_pressure_graph.gph ISH_graph.gph diabetes_graph.gph high_chol_graph.gph high_TG_graph.gph low_HDL_graph.gph, rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(36)

graph use $DTA/decomp1.gph
graph export $DTA/BE_CD1.png, as(png) width(3000) height(3600) replace

//// graph for subclinical outcomes
foreach x in BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln {
tw (bar coef_`x'_ _j if _j==0, color("68 1 84") mlab(ml_`x') mlabsize(4) mlabcolor(black)mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==1, color("59 82 139") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==2, color("33 145 140") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==3, color("94 201 98") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)) (bar coef_`x'_ _j if _j==4, color("253 231 37") mlab(ml_`x') mlabsize(4) mlabcolor(black) mlabcolor(black) mlabp(12) mlabg(-1)), xsize(10) ysize(12) legend(off) xtitle("") ytitle("") xlabel(0 "Total" 1 "PIE" 2 "INTref" 3 "INTmed"  4 "CDE" ,  angle(45) nogrid) ylabel(-0.2(0.1)0.2, nogrid) title("{bf:`x'}", size(medium)) 

graph save `x'_graph.gph, replace
}

graph combine BMI_ln_graph.gph Midbp_ln_graph.gph PP_ln_graph.gph glucose_ln_graph.gph cholestrol_ln_graph.gph Triglyceride_ln_graph.gph chol_HDL_ln_graph.gph , rows(3) cols(3)  imargin(0 0 0 0)  xsize(30) ysize(36)

graph use $DTA/decomp2.gph
graph export $DTA/BE_CD2.png, as(png) width(3000) height(3600) replace

