/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the main analysis
*/

clear
set more off
macro define DTA "LLE"
use $DTA/BE_clean, replace

keep if year>=2012 & year<=2021

****** Employment trend *****
preserve

collapse (mean) emp , by(year age_ca)
tw  (scatter emp  year if age_ca==3, color(black)) (qfit emp year if year<2017 & age_ca==3, lcolor(cranberry) range(2012 2017)) (qfit emp year if year>=2017 & age_ca==3, lcolor(cranberry)), xline(2017) xlabel(2012(1)2021, nogrid) ylabel(0.2(0.1)0.5, nogrid) legend(off) xtitle("Year") ytitle("Employment rate") xsize(12) ysize(12) plotregion(lcolor(black)) title("(a) Aged 65 or more")
graph export $DTA/FS_1.png, as(png) replace 	


tw  (scatter emp  year if age_ca==2, color(black)) (lfit emp year if year<2017 & age_ca==2, lcolor(cranberry) range(2012 2017)) (lfit emp year if year>=2017 & age_ca==2, lcolor(cranberry)), xline(2017) xlabel(2012(1)2021, nogrid) ylabel(0.6(0.1)0.9, nogrid) legend(off) xtitle("Year") ytitle("Employment rate") xsize(12) ysize(12) plotregion(lcolor(black)) title("(b) Aged 40-64")
graph export $DTA/FS_2.png, as(png) replace 	

restore

******* density graphs ******* 
foreach x in BMI Midbp PP glucose cholestrol chol_HDL Triglyceride {
kdensity `x', lcolor(navy) xsize(12) ysize(12) note("") xlabel(, nogrid) ylabel(, nogrid) title("`x'", size(medium)) xtitle("")

graph save $DTA/kd_`x'.gph, replace	
}	

graph combine ///
    $DTA/kd_BMI.gph            ///
    $DTA/kd_Midbp.gph          ///
    $DTA/kd_PP.gph             ///
    $DTA/kd_glucose.gph        ///
    $DTA/kd_cholestrol.gph     ///
    $DTA/kd_Triglyceride.gph  ///
    $DTA/kd_chol_HDL.gph,       ///
    rows(3) cols(3) imargin(0 0 0 0) xsize(36) ysize(36)

//graph use $DTA/kdensity.gph

graph export $DTA/kdensity.png, as(png) width(3600) height(3600) replace
 
******* First stage equation  ******
*** Note: using ivreg2 to calculate F-stat and LM statistic
/// RKD 
ivreg2 obesity (emp = c.year_c#1.cut) c.year_c i.cut  if dum==1, robust first

/// RKD + DID 
ivreg2 obesity (emp = c.year_c#1.cut#1.dum) c.year_c i.cut i.dum c.year_c#i.cut c.year_c#i.dum i.dum#i.cut , robust first 

******** balance check *******
preserve 
tab edu, g(edu)
tab region, g(region)

local v = 0 
foreach x in sex edu1 edu2 edu3 marital rural region1 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 height_std HP_hist HL_hist IHD_hist STR_hist DM_hist  {
local v = `v'+1

//// slope change among aged 65+ 
reg `x' i.cut#c.year_c i.cut c.year_c  if dum==1, robust
scalar coef1_`v' = r(table)[1, 2]
gen coef1_`v' = coef1_`v'
scalar se1_`v' = r(table)[2, 2]
gen se1_`v' = se1_`v'

/// the difference in slope change between aged 65+ and aged 40-64
reg `x' i.cut#c.year_c#i.dum i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum, robust
scalar coef2_`v' = r(table)[1, 4]
gen coef2_`v' = coef2_`v'
scalar se2_`v' = r(table)[2, 4]
gen se2_`v' = se2_`v'

}

//// graph for balance check
collapse (mean) coef1_* se1_* coef2_* se2_*
gen id =1 
reshape long coef1_ se1_ coef2_ se2_, i(id)

forval x = 1/2 { 
gen upper`x' = coef`x'_ + 1.96 *se`x'_ 
gen lower`x' = coef`x'_ - 1.96 *se`x'_	
}

replace _j = _j +1
replace _j = _j +1 if _j>=8
replace _j = _j +1 if _j>=20
replace _j = _j +1 if _j>=22


gen _j1 = _j-0.1
gen _j2 = _j+0.1 

graph set window fontface "Times New Roman"

tw (scatter coef1 _j1, color("122 209 81")) (rcap upper1 lower1 _j1, lcolor("122 209 81")) (scatter coef2 _j2, color("65 68 135")) (rcap upper2 lower2 _j2, lcolor("65 68 135")) ,  ylabel(-0.1(0.02)0.1, nogrid labsize(small)) xsize(16) ysize(12) xlabel(1 "{bf:Sociodemographic}" 2 "Women" 3 "Middle school or less" 4 "High school" 5 "College" 6 "Married" 7 "Rural" 8 "{bf:Region}" 9 "Seoul" 10 "Busan" 11 "Daegu" 12 "Incheon" 13 "Daejeon" 14 "Ulsan" 15 "Gyeonggi" 16 "Gangwon" 17 "Choongcheong" 18 "Jeonnam" 19 "Gyungsang" 20 "{bf:Physical}" 21 "Height" 22 "{bf:Medical history (family)}" 23 "High pressure" 24 "Hyperlipidemia" 25 "Ischemic heart disease" 26 "Stroke" 27 "Diabetes", nogrid angle(45) labsize(small)) ytitle("") yline(0) legend(order(1 "RKD" 3 "RKD+DiD") pos(1) col(2)) xsc(r(0.5, .)) ytitle("Effect estimates") plotregion(lcolor(black)) 

graph export $DTA/BE_balance.png, as(png) width(1600) height(1200) replace
restore

***** parellel trend assumption (event study)*****

preserve 
///// using 2016 as a reference period
reg emp ib(2016).year##i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust

forval i = 1/10 {
local v = (`i'+6)*2  
scalar coef_`i'= r(table)[1, `v']
scalar se_`i'= r(table)[2, `v']	
}


clear 
set obs 1 
forval i=1/10 {
	gen coef_`i'=coef_`i'
	gen se_`i'=se_`i'	
}
gen id = 1 
reshape long coef_ se_, i(id)
gen upper = coef_ + 1.96*se_ 
gen lower = coef_ - 1.96*se_ 

gen _j2 = _j+2011
labmask _j, val(_j2)

graph set window fontface "times new roman"
tw  (rcap upper lower _j, color("65 68 135")) (scatter coef_ _j, color("65 68 135")) , xlabel(1(1)10, nogrid valuelabels angle(45)) ylabel(-0.3(0.1)0.3, nogrid) legend(off) xsize(12) ysize(12) xtitle("Year") ytitle("Effect estimates") yline(0) plotregion(lcolor(black))

graph export $DTA/BE_parallel.png, as(png) width(1200) height(1200) replace

restore

***** Main analysis (biomarkers) *****
preserve
local v = 0 
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln Triglyceride_ln chol_HDL_ln  {
local v = `v'+1 

* baseline model
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum, robust
scalar coef1_`v' = _b[emp]
scalar se1_`v'= _se[emp]

* covariates adjusted
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust

scalar coef2_`v' = _b[emp]
scalar se2_`v'= _se[emp]
}

clear 
set obs 1 
forval x = 1/14 {
	forval i = 1/2 {
		gen coef`i'_`x' = coef`i'_`x'
		gen se`i'_`x' = se`i'_`x'
	}
}
gen id =1 
reshape long coef1_ se1_ coef2_ se2_, i(id)
save $DTA/BE_main_1, replace

/// main graph (biomarkers)
use $DTA/BE_main_1, replace
forval x = 1/2 {
gen upper_`x' = coef`x'_ +1.96*se`x'_ 
gen lower_`x' = coef`x'_ -1.96*se`x'_ 	
}
cap drop _j1 
cap drop _j2 

replace _j = _j+1 
replace _j = _j+1 if _j>=9

gen _j1 = _j-0.1
gen _j2 = _j+0.1 

tw (scatter  coef1_ _j1, color("33 145 130") msize(small)) (rcap upper_1 lower_1 _j1, lcolor("33 145 130")) (scatter  coef2_ _j2, color("68 1 84") msize(small)) (rcap upper_2 lower_2 _j2, lcolor("68 1 84")),ylabel(-1(0.5)1, nogrid format(%4.1f)) xsize(16) ysize(12) xlabel(1 "{bf: Clinical}" 2 "Obesity" 3 "Hypertension" 4 "ISH" 5 "Diabetes" 6 "Hypercholesterolemia" 7 "Hypertriglyceridemia" 8 "Hypo-HDL"  9 "{bf:Subclinical}" 10 "BMI" 11 "Mid-BP" 12 "Pulse pressure"  13 "Glucose" 14 "Cholestrol" 15 "Triglyceride" 16 "HDL-cholesterol" , nogrid angle(45) labsize(small)) ytitle("") yline(0) legend(order(1 "Baseline" 3 "Covariates adj.") pos(1) col(3) ) ytitle("Effect estimates") plotregion(lcolor(black)) 

graph export $DTA/main1.png, width(1600) height(1200) as(png) replace 

restore 

*************** Reduced form equations ***************
preserve 

local v = 0 
foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL BMI_ln Midbp_ln PP_ln glucose_ln cholestrol_ln chol_HDL_ln Triglyceride_ln  {
local v = `v'+1 

//// slope change among 65+
reg `x' i.cut#c.year_c i.cut c.year_c i.cut#c.year_c  if dum == 1, robust
scalar coef_1_`v' = _b[1.cut#c.year_c]
scalar se_1_`v' = _se[1.cut#c.year_c]

//// slope change among 65+
reg `x' i.cut#c.year_c i.cut c.year_c i.cut#c.year_c  if dum == 0, robust
scalar coef_2_`v' = _b[1.cut#c.year_c]
scalar se_2_`v' = _se[1.cut#c.year_c]

/// diff in slople change 
reg `x' i.cut#c.year_c#i.dum i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum, robust 
scalar coef_3_`v' = _b[1.cut#c.year_c#1.dum]
scalar se_3_`v' = _se[1.cut#c.year_c#1.dum]

//// average partial effect 
margins , expression(2 * _b[1.cut#c.year_c#1.dum])
scalar coef_4_`v'  = r(table)[1, 1]
scalar se_4_`v' = r(table)[2, 1]
}

clear 
set obs 1
forval x = 1/14 {
	forval i = 1/4 {
		gen coef_`i'_`x'=coef_`i'_`x'
		gen se_`i'_`x'=se_`i'_`x'
	}
}

gen id =1 
reshape long coef_1_ se_1_ coef_2_ se_2_  coef_3_ se_3_  coef_4_ se_4_ , i(id)

forval x = 1/4 {
gen test_`x' = coef_`x'_ / se_`x'_
gen p_`x' = 2*(1 - normprob(abs(test_`x')))	
}

forval x = 1/4 {
	gen upper_`x' = coef_`x'_+1.96*se_`x'_
	gen lower_`x' = coef_`x'_-1.96*se_`x'_
}

save $DTA/reduced, replace
restore 

***** Main analysis (health behaviors) *****

preserve
local v = 0 
foreach x in move exercise_1 exercise_2 sedentary_std alchol_freq_std alchol_number_std alchol_high_freq_std smoking smoking_number_std insurance med_serv  {
local v = `v'+1 

* baseline model
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum
scalar coef1_`v' = _b[emp]
scalar se1_`v'= _se[emp]

* covariates adjusted
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist

scalar coef2_`v' = _b[emp]
scalar se2_`v'= _se[emp]
}

clear 
set obs 1 
forval x = 1/11 {
	forval i = 1/2 {
		gen coef`i'_`x' = coef`i'_`x'
		gen se`i'_`x' = se`i'_`x'
	}
}
gen id =1 
reshape long coef1_ se1_ coef2_ se2_, i(id)
save $DTA/BE_main_2, replace

/// main graph (health behaviors) 
use $DTA/BE_main_2, replace

forval x = 1/2 {
gen upper_`x' = coef`x'_ +1.96*se`x'_ 
gen lower_`x' = coef`x'_ -1.96*se`x'_ 	
}
cap drop _j1 
cap drop _j2 
cap drop _j 

gen _j = _n

replace _j = _j+1 
replace _j= _j+1  if _j>=6
replace _j = _j+1 if _j>=12

gen _j1 = _j-0.2 
gen _j2 = _j+0.2 

tw  (scatter  coef1_ _j, color("33 145 130") msize(small)) (rcap upper_1 lower_1 _j, lcolor("33 145 130")) (scatter  coef2_ _j2, color("68 1 84") msize(small)) (rcap upper_2 lower_2 _j2, lcolor("68 1 84")),ylabel(-4(1)3, nogrid) xsize(16) ysize(12) xlabel(1 "{bf:Physical activities}" 2 "Movement" 3 "Vigorous activity" 4 "Moderate activity" 5 "Sedentary hours" 6 "{bf:Health behaviors}" 7 "Alcohol frequency" 8 "Alcohol amount" 9 "Binge drinking" 10 "Smoking status" 11 "Cigarettes per day" 12 "{bf:Medical access}" 13 "Sponsored/private ins." 14 "Health examination" , nogrid angle(45) labsize(small)) ytitle("") yline(0) legend(order(1 "Baseline" 3 "Covariates adj.") pos(1) col(3)) ytitle("Effect estimates") plotregion(lcolor(black)) 

graph export $DTA/Main2.png, width(1600) height(1200) as(png) replace 
restore

***** Main analysis (dietary patterns) *****
local v = 0 
foreach x in HEI_std HEI_BR_std HEI_CEREAL_std HEI_TFRUIT_std HEI_FFRUIT_std HEI_TVEG_std HEI_VEG_std HEI_PROTF_std HEI_DAIRY_std HEI_SFA_std HEI_NA_std HEI_CHO_std HEI_FAT_std HEI_EN_std {
local v = `v'+1 

* baseline model
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum , robust
scalar coef1_`v' = _b[emp]
scalar se1_`v'= _se[emp]

* covariates adjusted
ivreg2 `x' (emp = i.cut#c.year_c#i.dum) i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum i.sex i.edu i.marital i.rural i.region c.height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust

scalar coef2_`v' = _b[emp]
scalar se2_`v'= _se[emp]
}

preserve
clear 
set obs 1 
forval x = 1/14 {
	forval i = 1/2 {
		gen coef`i'_`x' = coef`i'_`x'
		gen se`i'_`x' = se`i'_`x'
	}
}
gen id =1 
reshape long coef1_ se1_ coef2_ se2_, i(id)
save $DTA/BE_main_3, replace

//// graph (dietary patterns)
use $DTA/main3, replace
forval x = 1/2 {
gen upper_`x' = coef`x'_ +1.96*se`x'_ 
gen lower_`x' = coef`x'_ -1.96*se`x'_ 	
}
cap drop _j1 
cap drop _j2 
cap drop _j 

gen _j = _n

gen _j1 = _j-0.1
gen _j2 = _j+0.1

* Total score/ Have breakfast/ mixed grains/ Fresh fruits / Vegetables / Meat, fish, eggs and beans / Milk and milk products / Saturated fatty acid / Sodium / Carbohydrate / 
tw  (scatter  coef1_ _j1, color("33 145 130") msize(small)) (rcap upper_1 lower_1 _j1, lcolor("33 145 130")) (scatter  coef2_ _j2, color("68 1 84") msize(small)) (rcap upper_2 lower_2 _j2, lcolor("68 1 84")),ylabel(-4(2)4, nogrid) xsize(16) ysize(12) xlabel(1 "Total score" 2 "Breakfast" 3 "Mixed grains" 4 "Total fruits" 5 "Fresh fruits" 6 "Total vegetables" 7 "Vegetables" 8 "Meat, fish, eggs, and beans" 9 "Milk and milk products" 10 "Saturated fatty acid" 11 "Sodium" 12 "Carbohydrate" 13 "FAT" 14 "Energy intake", nogrid angle(45) labsize(small)) ytitle("") yline(0) legend(order(1 "Baseline" 3 "Covariates adj.") pos(1) col(3)) ytitle("Effect estimates") plotregion(lcolor(black)) 

graph export $DTA/Main3.png, width(1600) height(1200) as(png) replace 
restore
