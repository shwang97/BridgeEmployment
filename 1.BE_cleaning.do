/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for preparing data for the analysis
*/

clear
set more off
macro define DTA "LLE"

use $DTA/HN_10, replace 
forval x=11/22 {
append using $DTA/HN_`x'	
}

***** baseline characteristics ******
gen age_ca =1 if age>=25 & age<= 39
replace age_ca =2  if age>=40 & age<=64
replace age_ca =3 if age>=65

gen emp = 1 if EC1_1 ==1 
replace emp =0 if EC1_1 == 2

gen height = HE_ht
replace region = region +1 if year<=2015 & region>=8
recode region (8=11)
recode region (12=11)(14=13)(16=15)

gen marital = 2-marri_1 if marri_1!=9 
gen rural = town_t -1 

* medical history (family)
foreach x in HP HL IHD STR DM {
gen `x'_hist=0 if HE_`x'fh1~=.& HE_`x'fh2~=. & HE_`x'fh3~=.
replace `x'_hist =1 if HE_`x'fh1==1|HE_`x'fh2==1|HE_`x'fh3==1 
}

****** outcome variables ******* 

* biomarkers
gen glucose= HE_glu
gen glycated = HE_HbA1c
gen cholestrol = HE_chol
gen chol_HDL = HE_HDL_st2
gen Triglyceride = HE_TG 

gen BMI = HE_BMI 
gen Midbp= (HE_dbp +HE_sbp)/2 
gen PP = HE_sbp - HE_dbp 

gen obesity = BMI>=25 if BMI~=. 
gen high_pressure = HE_sbp>140 | HE_dbp>90 if HE_sbp~=. & HE_dbp~=. /// hypertension
gen ISH = HE_sbp >=140  & HE_dbp<90 if HE_sbp!=. & HE_sbp!=. ///Isolated systolic hypertension

gen diabetes =0 if glucose!=. & HE_fst!=. 
replace diabetes = 1 if glucose>=126 & glucose!=. & HE_fst>=8 &HE_fst!=. 
replace diabetes =1 if glucose>=200 & glucose!=. & HE_fst<8 &HE_fst!=. 

gen high_chol = cholestrol>=240 if cholestrol~=. /// high cholestrol
 
gen high_TG =0 if Triglyceride!=. & HE_fst!=.
replace high_TG = 1 if Triglyceride>=150 & Triglyceride!=. & HE_fst>=8 &HE_fst!=. 
replace high_TG = 1 if Triglyceride>=175 & Triglyceride!=. & HE_fst<8 &HE_fst!=. 

gen low_HDL=0 if chol_HDL!=. 
replace low_HDL =1 if chol_HDL<40 & chol_HDL!=. 


* health behaviors
forval x = 7/8 {
local v = `x'-6
gen exercise_`v' = BE3_`x'1==1 | BE3_`x'5==1 if BE3_`x'1~=88 & BE3_`x'1~=99 & BE3_`x'5~=88 & BE3_`x'5~=99
}

gen sedentary = BE8_1 if BE8_1~=88 & BE8_1~=99
gen move = 2-BE3_91 if BE3_91 ~=8 & BE3_91 ~=9
gen alchol_freq = BD1_11 if BD1_11~=8 &BD1_11~=9

forval x = 1/2 {
gen alchol_high_`x' = BD2_3`x' if  BD2_3`x'~=9
replace alchol_high_`x'=0 if BD2_3`x'==8
}

gen alchol_number = BD2_1 if BD2_1~=9
replace alchol_number =0 if BD2_1 ==8

gen alchol_high_freq =.
replace alchol_high_freq = alchol_high_1 if sex==1 
replace alchol_high_freq = alchol_high_2 if sex==2 

cap drop smoking
gen smoking = BS3_1 if BS3_1~=9 
recode smoking (1/2=1) (3/8=0)

gen smoking_number= BS3_2 if BS3_2~=999
replace smoking_number=0 if BS3_2==888

gen compmed = tins 
recode compmed (20=1)(10=0)(30=0)(99=0)

gen private = npins if npins~=99 & npins~=9
replace private = 2-private

gen insurance = 1 if compmed ==1 | private==1
replace insurance =0 if compmed==0 & private==0

gen med_serv = BH1 if BH1~=9
replace med_serv = 2-med_serv
rename med_serv med_serv_std

****** RKD variables ****
gen cut= 1 if year>=2017
replace cut =0 if year<2017
gen year_c = year-2017

gen dum = 1 if age_ca ==3 
replace dum = 0 if age_ca==2

**** sample selection
foreach x in sex edu marital rural region height HP_hist HL_hist IHD_hist STR_hist DM_hist {
keep if `x'!=. 	
}
keep if dum!=. 

***** log transformation and standardization 
egen height_std = std(height) 
foreach x in BMI Midbp PP glucose glycated cholestrol chol_HDL Triglyceride {
egen `x'_std = std(`x') 
}

foreach x in sedentary alchol_freq alchol_number alchol_high_freq smoking_number {
cap drop `x'_std
egen `x'_std = std(`x') 
}

foreach x in HEI HEI_BR HEI_CEREAL HEI_TFRUIT HEI_FFRUIT HEI_TVEG HEI_VEG HEI_PROTF HEI_DAIRY HEI_SFA HEI_NA HEI_CHO HEI_FAT HEI_EN {
egen `x'_std = std(`x') 
}

save $DTA/BE_clean, replace


