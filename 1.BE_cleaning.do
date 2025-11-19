/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for preparing data for the analysis
*/

clear
set more off
macro define DTA ""

use $DTA/HN_10, replace 
forval x=11/22 {
append using $DTA/HN_`x'	
}

*---------------------------------------------------------------*
* Baseline Characteristics
*---------------------------------------------------------------*

* Age categories
gen age_ca = .
replace age_ca = 1 if inrange(age, 25, 39)
replace age_ca = 2 if inrange(age, 40, 64)
replace age_ca = 3 if age >= 65

* Employment indicator
gen emp = .
replace emp = 1 if EC1_1 == 1
replace emp = 0 if EC1_1 == 2

* Height
gen height = HE_ht

* Region adjustments (pre-2015)
replace region = region + 1 if year <= 2015 & region >= 8
recode region (8 = 11)
recode region (12 = 11) (14 = 13) (16 = 15)

* Marital status
gen marital = 2 - marri_1 if marri_1 != 9

* Rural indicator
gen rural = town_t - 1

*---------------------------------------------------------------*
* Family medical history indicators
* Variables: HP = Hypertension, HL = Hyperlipidemia, 
* IHD = Ischemic heart disease, STR = Stroke, DM = Diabetes
*---------------------------------------------------------------*
foreach x in HP HL IHD STR DM {
    gen `x'_hist = 0 if !missing(HE_`x'f h1, HE_`x'f h2, HE_`x'f h3)
    replace `x'_hist = 1 if HE_`x'fh1 == 1 | HE_`x'fh2 == 1 | HE_`x'fh3 == 1
}

*---------------------------------------------------------------*
* Outcome Variables
*---------------------------------------------------------------*

*-------------------------------*
* Subclinical Symptoms
*-------------------------------*
gen glucose      = HE_glu
gen glycated     = HE_HbA1c
gen cholestrol   = HE_chol        // total cholesterol
gen chol_HDL     = HE_HDL_st2     // HDL cholesterol
gen Triglyceride = HE_TG

gen BMI = HE_BMI
gen PP  = HE_sbp - HE_dbp         // pulse pressure

*-------------------------------*
* Clinical Symptoms
*-------------------------------*

* Obesity (BMI >= 25)
gen obesity = BMI >= 25 if !missing(BMI)

* Hypertension (SBP > 140 or DBP > 90)
gen high_pressure = (HE_sbp > 140 | HE_dbp > 90) ///
                    if !missing(HE_sbp, HE_dbp)

* Isolated systolic hypertension (SBP >= 140 and DBP < 90)
gen ISH = (HE_sbp >= 140 & HE_dbp < 90) ///
          if !missing(HE_sbp, HE_dbp)

* Diabetes
gen diabetes = 0 if !missing(glucose, HE_fst)
replace diabetes = 1 if glucose >= 126 & !missing(glucose) ///
                      & HE_fst >= 8   & !missing(HE_fst)
replace diabetes = 1 if glucose >= 200 & !missing(glucose) ///
                      & HE_fst < 8    & !missing(HE_fst)

* High total cholesterol
gen high_chol = (cholestrol >= 240) if !missing(cholestrol)

* High triglycerides
gen high_TG = 0 if !missing(Triglyceride, HE_fst)
replace high_TG = 1 if Triglyceride >= 150 & !missing(Triglyceride) ///
                      & HE_fst >= 8       & !missing(HE_fst)
replace high_TG = 1 if Triglyceride >= 175 & !missing(Triglyceride) ///
                      & HE_fst < 8        & !missing(HE_fst)

* Low HDL
gen low_HDL = 0 if !missing(chol_HDL)
replace low_HDL = 1 if chol_HDL < 40 & !missing(chol_HDL)

*---------------------------------------------------------------*
* Health-Related Behaviors
*---------------------------------------------------------------*

* Walking frequency
recode BE3_31 (99 = .)
replace BE3_31 = BE3_31 - 1

cap drop walking
gen walking = BE3_31

* Sedentary time 
gen sedentary = BE8_1 if BE8_1 != 88 & BE8_1 != 99

* Alcohol consumption frequency 
gen alchol_freq = BD1_11 if BD1_11 != 8 & BD1_11 != 9

* High-risk drinking frequency (sex-specific thresholds)
forvalues x = 1/2 {
    gen alchol_high_`x' = BD2_3`x' if BD2_3`x' != 9
    replace alchol_high_`x' = 0 if BD2_3`x' == 8
}
* Sex-specific high-risk drinking variable
gen alchol_high_freq = .
replace alchol_high_freq = alchol_high_1 if sex == 1
replace alchol_high_freq = alchol_high_2 if sex == 2

* Number of drinks per occasion
gen alchol_number = BD2_1 if BD2_1 != 9
replace alchol_number = 0 if BD2_1 == 8

* Smoking Behavior
cap drop smoking
gen smoking = BS3_1 if BS3_1 != 9
recode smoking (1/2 = 1) (3/8 = 0)   

gen smoking_number = BS3_2 if BS3_2 != 999
replace smoking_number = 0 if BS3_2 == 888

* Sponsored insurance
gen compmed = tins
recode compmed (20 = 1) (10 = 0) (30 = 0) (99 = 0)

* Private insurance 
gen private = npins if npins != 99 & npins != 9
replace private = 2 - private   

* Any insurance coverage (sponsored or private)
gen insurance = 1 if compmed == 1 | private == 1
replace insurance = 0 if compmed == 0 & private == 0

* Medical examination
gen med_serv = BH1 if BH1 != 9
replace med_serv = 2 - med_serv   
rename med_serv med_serv_std

*---------------------------------------------------------------*
* RKD Variables
*---------------------------------------------------------------*

* Treatment cutoff indicator (2017 policy change)
gen cut = (year >= 2017)
replace cut = 0 if year < 2017

* Centered running variable (years relative to 2017)
gen year_c = year - 2017

* Eligibility dummy: age 65+ vs. age 40-64
gen dum = .
replace dum = 1 if age_ca == 3      // 65+
replace dum = 0 if age_ca == 2      // 40-64 (comparison group)

*---------------------------------------------------------------*
* Sample Selection
* 1. Restrict to survey years 2012-2021
* 2. Exclude observations with missing RKD eligibility dummy (dum)
* 3. Exclude individuals with missing pre-treatment characteristics
*---------------------------------------------------------------*

* Study period
keep if inrange(year, 2012, 2021)

* RKD eligibility dummy must be observed
keep if !missing(dum)

* Exclude observations with missing pre-treatment characteristics
foreach var in sex edu marital rural region height ///
                HP_hist HL_hist IHD_hist STR_hist DM_hist {
    keep if !missing(`var')
}

*---------------------------------------------------------------*
* Log Transformations and Standardization
*---------------------------------------------------------------*

* Standardize height
egen height_std = std(height)

* Log-transform continuous clinical/subclinical measures
foreach var in BMI PP glucose glycated cholestrol chol_HDL Triglyceride {
    gen `var'_ln = log(`var')
}

* Log-transform blood pressure measures from original variables
foreach var in sbp dbp {
    gen `var'_ln = log(HE_`var')
}

* Log-transform behavioral variables (with +1 to handle zeros)
foreach var in walking sedentary alchol_freq alchol_number ///
                alchol_high_freq smoking_number HEI {
    cap drop `var'_ln
    gen `var'_ln = log(`var' + 1)
}

save $DTA/, replace


