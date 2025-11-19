/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the main analysis
*/

******* Descriptive statistics *****
clear
set more off
macro define DTA "LLE"
use $DTA/, replace

*---------------------------------------------------------------*
* Define variable groups for descriptive tables
*---------------------------------------------------------------*

* Clinical & subclinical health outcomes
local health_vars ///
    i.obesity i.high_pressure i.ISH i.diabetes i.high_chol i.high_TG ///
    BMI HE_sbp HE_dbp PP ///
    glucose cholestrol chol_HDL Triglyceride ///
    
* Health-related behaviors
local behavior_vars ///
    walking sedentary i.pa_aerobic ///
    alchol_freq alchol_number alchol_high_freq ///
    i.smoking smoking_number i.insurance i.med_serv HEI

* Baseline characteristics
local baseline_vars ///
    i.emp i.sex i.edu i.marital i.rural height ///
    i.HP_hist i.HL_hist i.IHD_hist i.STR_hist i.DM_hist

*---------------------------------------------------------------*
* Descriptive Statistics: Clinical & Subclinical Health Outcomes
*---------------------------------------------------------------*

* Means & frequencies
dtable `health_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvfrequency)) ///
    continuous(, statistic(mean)) ///
    export($DTA/, replace)

* SDs & percentages
dtable `health_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvpercent)) ///
    continuous(, statistic(sd)) ///
    export($DTA/, replace)

*---------------------------------------------------------------*
* Descriptive Statistics: Health-Related Behaviors
*---------------------------------------------------------------*

* Means & frequencies
dtable `behavior_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvfrequency)) ///
    continuous(, statistic(mean)) ///
    export($DTA/, replace)

* SDs & percentages
dtable `behavior_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvpercent)) ///
    continuous(, statistic(sd)) ///
    export($DTA/, replace)

*---------------------------------------------------------------*
* Descriptive Statistics: Baseline Characteristics
*---------------------------------------------------------------*

* Means & frequencies
dtable `baseline_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvfrequency)) ///
    continuous(, statistic(mean)) ///
    export($DTA/, replace)

* SDs & percentages
dtable `baseline_vars', ///
    by(dum) ///
    nformat(%16.2fc mean sd) ///
    factor(, statistic(fvpercent)) ///
    continuous(, statistic(sd)) ///
    export($DTA/, replace)

*---------------------------------------------------------------*
* Figure3: Employment Rate by Calendar Year 
*---------------------------------------------------------------*

use $DTA/BE_clean, replace

* Locals for margins grid and number of points
local pre_grid  -5(0.1)0
local post_grid 0(0.1)5
local npts      51

*---------------------------------------------------------------*
* 1. Age 65+ (age_ca == 3)
*---------------------------------------------------------------*

quietly reg emp c.year_c#1.cut c.year_c i.cut if age_ca==3, robust

* Pre-2017 (cut = 0)
quietly margins, at(year_c=(`pre_grid') cut=0)

forvalues i = 1/`npts' {
    scalar coef1_1_`i'  = r(table)[1, `i']
    scalar lower1_1_`i' = r(table)[5, `i']
    scalar upper1_1_`i' = r(table)[6, `i']
}

* Post-2017 (cut = 1)
quietly margins, at(year_c=(`post_grid') cut=1)

forvalues i = 1/`npts' {
    scalar coef1_2_`i'  = r(table)[1, `i']
    scalar lower1_2_`i' = r(table)[5, `i']
    scalar upper1_2_`i' = r(table)[6, `i']
}

*---------------------------------------------------------------*
* 2. Age 40-64 (age_ca == 2)
*---------------------------------------------------------------*

quietly reg emp c.year_c#1.cut c.year_c i.cut if age_ca==2, robust

* Pre-2017 (cut = 0)
quietly margins, at(year_c=(`pre_grid') cut=0)

forvalues i = 1/`npts' {
    scalar coef2_1_`i'  = r(table)[1, `i']
    scalar lower2_1_`i' = r(table)[5, `i']
    scalar upper2_1_`i' = r(table)[6, `i']
}

* Post-2017 (cut = 1)
quietly margins, at(year_c=(`post_grid') cut=1)

forvalues i = 1/`npts' {
    scalar coef2_2_`i'  = r(table)[1, `i']
    scalar lower2_2_`i' = r(table)[5, `i']
    scalar upper2_2_`i' = r(table)[6, `i']
}

*---------------------------------------------------------------*
* 3. Collapse to yearly means and create fine grid (0.1 steps)
*---------------------------------------------------------------*

collapse (mean) emp, by(year age_ca cut)

preserve
    keep if year==2017              
    replace cut = 0
    replace emp = .
    tempfile add2017
    save `add2017'
restore

* Append the additional 2017 observations
append using `add2017'

* Expand each observation into 10 sub-steps (0.1 increments)
expand 10

* Within-group index (09)
bysort year age_ca cut: gen seq = _n - 1

* Fractional year: 2012, 2012.1, …
gen year_g = year + seq*0.1

sort age_ca year_g cut

*---------------------------------------------------------------*
* 4. Map stored scalars to coef/lower/upper for each group
*---------------------------------------------------------------*

* Helper indices to match scalars:
* d1: age group code; d2: pre/post (cut); d3: grid index
gen d1 = 4 - age_ca
gen d2 = cut + 1
gen d3 = (year-2012)*10 + seq + 1 if cut==0
replace d3 = (year-2017)*10 + seq + 1 if cut==1

gen coef  = .
gen lower = .
gen upper = .

foreach a in 1 2 {          // a = 1: age_ca==3; a = 2: age_ca==2
    foreach c in 1 2 {      // c = 1: cut==0; c = 2: cut==1
        forvalues y = 1/`npts' {
            replace coef  = coef`a'_`c'_`y'  if d1==`a' & d2==`c' & d3==`y'
            replace lower = lower`a'_`c'_`y' if d1==`a' & d2==`c' & d3==`y'
            replace upper = upper`a'_`c'_`y' if d1==`a' & d2==`c' & d3==`y'
        }
    }
}

sort age_ca cut year_g

* Keep only points within observation period
replace emp = . if year!=year_g
keep if year_g<=2021

*---------------------------------------------------------------*
* 5. Plot: Aged 65+ (Panel a)
*---------------------------------------------------------------*

gen ub1 = 0.5 
gen lb1 = 0.2
gen c   = 2017

twoway ///
    (rarea upper lower year_g if age_ca==3 & cut==0, ///
        fcolor(navy*0.3) lwidth(none)) ///
    (rarea upper lower year_g if age_ca==3 & cut==1, ///
        fcolor(navy*0.3) lwidth(none)) ///
    (lfit emp year if year<2017 & age_ca==3, ///
        lcolor(navy) range(2012 2017) lwidth(0.4)) ///
    (lfit emp year if year>=2017 & age_ca==3, ///
        lcolor(navy) lwidth(0.4)) ///
    (rspike ub1 lb1 c, ///
        lpattern(dash) lcolor(black) lwidth(thin)) ///
    (scatter emp year if age_ca==3, ///
        color(black)), ///
    xlabel(2012(1)2021, nogrid) ///
    ylabel(0.2(0.1)0.5, nogrid) ///
    legend(off) ///
    xtitle("Year") ///
    ytitle("Employment rate") ///
    xsize(12) ysize(12) ///
    plotregion(lcolor(black) margin(1 1 0 0)) ///
    title("(a) Aged 65 or more")

graph export "$DTA/.png", as(png) replace

*---------------------------------------------------------------*
* 6. Plot: Aged 40-64 (Panel b)
*---------------------------------------------------------------*

gen ub2 = 0.9
gen lb2 = 0.6

twoway ///
    (rarea upper lower year_g if age_ca==2 & cut==0, ///
        fcolor(navy*0.3) lwidth(none)) ///
    (rarea upper lower year_g if age_ca==2 & cut==1, ///
        fcolor(navy*0.3) lwidth(none)) ///
    (lfit emp year if year<2017 & age_ca==2, ///
        lcolor(navy) range(2012 2017) lwidth(0.4)) ///
    (lfit emp year if year>=2017 & age_ca==2, ///
        lcolor(navy) lwidth(0.4)) ///
    (rspike ub2 lb2 c, ///
        lpattern(dash) lcolor(black) lwidth(thin)) ///
    (scatter emp year if age_ca==2, ///
        color(black)), ///
    xline(2017) ///
    xlabel(2012(1)2021, nogrid) ///
    ylabel(0.6(0.1)0.9, nogrid) ///
    legend(off) ///
    xtitle("Year") ///
    ytitle("Employment rate") ///
    xsize(12) ysize(12) ///
    plotregion(lcolor(black) margin(1 1 0 0)) ///
    title("(b) Aged 40 to 64")

graph export "$DTA/.png", as(png) replace

*---------------------------------------------------------------*
* Table 2. First-Stage Equation: Checking KP-LM and KP rk F Statistics
* Note: Pseudo outcome variable (only first-stage of emp matters)
*---------------------------------------------------------------*

* Create pseudo outcome variable 
cap drop y_pseudo
gen double y_pseudo = rnormal()

* Common covariates
local covars i.sex i.edu i.marital i.rural i.region ///
             c.height HP_hist HL_hist IHD_hist STR_hist DM_hist

*---------------------------------------------------------------*
* 1. RKD Specifications (dum == 1, e.g., aged 65+)
*---------------------------------------------------------------*

* Baseline RKD first stage (parsimonious)
ivreg2 y_pseudo (emp = c.year_c#1.cut) ///
       c.year_c i.cut ///
       if dum==1, robust first

* RKD first stage with additional covariates
ivreg2 y_pseudo (emp = c.year_c#1.cut) ///
       c.year_c i.cut `covars' ///
       if dum==1, robust first

*---------------------------------------------------------------*
* 2. RKD + DiD Specifications
*---------------------------------------------------------------*

* RKD + DiD first stage (no additional covariates)
ivreg2 y_pseudo (emp = c.year_c#1.cut#1.dum) ///
       c.year_c i.cut i.dum ///
       c.year_c#i.cut c.year_c#i.dum i.dum#i.cut, ///
       robust first

* RKD + DiD first stage with additional covariates
ivreg2 y_pseudo (emp = c.year_c#1.cut#1.dum) ///
       c.year_c i.cut i.dum ///
       c.year_c#i.cut c.year_c#i.dum i.dum#i.cut ///
       `covars', ///
       robust first


*---------------------------------------------------------------*
* Figure 4: Reduced-Form Graphs for Clinical and Subclinical Outcomes
*---------------------------------------------------------------*

* Outcomes to loop over
local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln

* Grid for margins
local pre_grid  -5(0.1)0
local post_grid 0(0.1)5
local npts      51

* Year labels for x-axis
local yearlabs  2012 "12" 2013 "13" 2014 "14" 2015 "15" 2016 "16" ///
                2017 "17" 2018 "18" 2019 "19" 2020 "20" 2021 "21"
*---------------------------------------------------------------*
* 1. Loop over each outcome
*---------------------------------------------------------------*

foreach x of local outcomes {

    clear
    set more off
    use $DTA/, replace

    *-----------------------------------------------------------*
    * Margins for age 65+ (age_ca == 3)
    *-----------------------------------------------------------*

    quietly reg `x' c.year_c#1.cut c.year_c i.cut if age_ca==3, robust

    * Pre-2017 (cut = 0)
    quietly margins, at(year_c=(`pre_grid') cut=0)

    forvalues i = 1/`npts' {
        scalar coef1_1_`i'  = r(table)[1, `i']
        scalar lower1_1_`i' = r(table)[5, `i']
        scalar upper1_1_`i' = r(table)[6, `i']
    }

    * Post-2017 (cut = 1)
    quietly margins, at(year_c=(`post_grid') cut=1)

    forvalues i = 1/`npts' {
        scalar coef1_2_`i'  = r(table)[1, `i']
        scalar lower1_2_`i' = r(table)[5, `i']
        scalar upper1_2_`i' = r(table)[6, `i']
    }

    *-----------------------------------------------------------*
    * Margins for age 40-64 (age_ca == 2)
    *-----------------------------------------------------------*

    quietly reg `x' c.year_c#1.cut c.year_c i.cut if age_ca==2, robust

    * Pre-2017 (cut = 0)
    quietly margins, at(year_c=(`pre_grid') cut=0)

    forvalues i = 1/`npts' {
        scalar coef2_1_`i'  = r(table)[1, `i']
        scalar lower2_1_`i' = r(table)[5, `i']
        scalar upper2_1_`i' = r(table)[6, `i']
    }

    * Post-2017 (cut = 1)
    quietly margins, at(year_c=(`post_grid') cut=1)

    forvalues i = 1/`npts' {
        scalar coef2_2_`i'  = r(table)[1, `i']
        scalar lower2_2_`i' = r(table)[5, `i']
        scalar upper2_2_`i' = r(table)[6, `i']
    }

    *-----------------------------------------------------------*
    * Collapse to yearly means and create 0.1-year grid
    *-----------------------------------------------------------*

    collapse (mean) `x', by(year age_ca cut)

    preserve
        keep if year==2017              // keeps both age_ca==2 and 3
        replace cut = 0
        replace `x' = .
        tempfile add2017
        save `add2017'
    restore

    * Append additional 2017 obs with cut=0
    append using `add2017'

    * Expand each observation into 10 sub-steps (0.1 increments)
    expand 10

    * Within-group index (09)
    bysort year age_ca cut: gen seq = _n - 1

    * Fractional year (2012, 2012.1, …)
    gen year_g = year + seq*0.1

    sort age_ca year_g cut

    * Helper indices to match scalars:
    gen d1 = 4 - age_ca
    gen d2 = cut + 1
    gen d3 = (year-2012)*10 + seq + 1 if cut==0
    replace d3 = (year-2017)*10 + seq + 1 if cut==1

    gen coef  = .
    gen lower = .
    gen upper = .

    foreach a in 1 2 {
        foreach c in 1 2 {
            forvalues y = 1/`npts' {
                replace coef  = coef`a'_`c'_`y'  if d1==`a' & d2==`c' & d3==`y'
                replace lower = lower`a'_`c'_`y' if d1==`a' & d2==`c' & d3==`y'
                replace upper = upper`a'_`c'_`y' if d1==`a' & d2==`c' & d3==`y'
            }
        }
    }

    sort age_ca cut year_g

    replace `x' = . if year!=year_g
    keep if year_g<=2021

    gen c = 2017

    *-----------------------------------------------------------*
    * Panel (a): Aged 65+ (age_ca == 3)
    *-----------------------------------------------------------*

    quietly summarize `x' if inrange(year,2012,2021) & age_ca==3
    local ymin = r(min)
    local ymax = r(max)

    * Add padding (currently 100% of range, as in your code)
    local pad  = (`ymax' - `ymin') * 1
    local low  = round(`ymin' - `pad', 0.01)
    local high = round(`ymax' + `pad', 0.01)
    local gap  = (`high' - `low')/5

    gen ub1 = `high'
    gen lb1 = `low'

    tw ///
        (rarea upper lower year_g if age_ca==3 & cut==0, ///
            fcolor(navy*0.3) lwidth(none)) ///
        (rarea upper lower year_g if age_ca==3 & cut==1, ///
            fcolor(navy*0.3) lwidth(none)) ///
        (lfit `x' year if year<2017 & age_ca==3, ///
            lcolor(navy) range(2012 2017) lwidth(0.4)) ///
        (lfit `x' year if year>=2017 & age_ca==3, ///
            lcolor(navy) lwidth(0.4)) ///
        (rspike ub1 lb1 c, ///
            lpattern(dash) lcolor(black) lwidth(thin)) ///
        (scatter `x' year if age_ca==3, ///
            color(black)), ///
        xlabel(`yearlabs', nogrid labsize(small)) ///
        ylabel(`low'(`gap')`high', nogrid format(%4.2f) labsize(small)) ///
        legend(off) ///
        xtitle("") ytitle("") ///
        xsize(12) ysize(12) ///
        plotregion(lcolor(black) margin(1 1 0 0)) ///
        title("(a) Aged 65 or more", size(small)) ///
        graphregion(margin(0 1 0 0))

    graph save $DTA/`x'_1.gph, replace

    *-----------------------------------------------------------*
    * Panel (b): Aged 40-64 (age_ca == 2)
    *-----------------------------------------------------------*

    quietly summarize `x' if inrange(year,2012,2021) & age_ca==2
    local ymin = r(min)
    local ymax = r(max)

    * Match center, preserve overall range from panel (a)
    local mid   = (`ymax' + `ymin')/2
    local range = (`high' - `low')/2
    local low   = round(`mid' - `range', 0.01)
    local high  = round(`mid' + `range', 0.01)

    gen ub2 = `high'
    gen lb2 = `low'

    tw ///
        (rarea upper lower year_g if age_ca==2 & cut==0, ///
            fcolor(navy*0.3) lwidth(none)) ///
        (rarea upper lower year_g if age_ca==2 & cut==1, ///
            fcolor(navy*0.3) lwidth(none)) ///
        (lfit `x' year if year<2017 & age_ca==2, ///
            lcolor(navy) range(2012 2017) lwidth(0.4)) ///
        (lfit `x' year if year>=2017 & age_ca==2, ///
            lcolor(navy) lwidth(0.4)) ///
        (rspike ub2 lb2 c, ///
            lpattern(dash) lcolor(black) lwidth(thin)) ///
        (scatter `x' year if age_ca==2, ///
            color(black)), ///
        xlabel(`yearlabs', nogrid labsize(small)) ///
        ylabel(`low'(`gap')`high', nogrid format(%4.2f) labsize(small)) ///
        legend(off) ///
        xtitle("") ytitle("") ///
        xsize(12) ysize(12) ///
        plotregion(lcolor(black) margin(1 1 0 0)) ///
        title("(b) Aged 40 to 64", size(small)) ///
        graphregion(margin(0 1 0 0))

    graph save $DTA/`x'_2.gph, replace
}

*---------------------------------------------------------------*
* 2. Combine per-outcome graphs into panels
*---------------------------------------------------------------*

foreach x of local outcomes {
    graph combine $DTA/`x'_1.gph $DTA/`x'_2.gph, ///
        col(2) xsize(20) ysize(12) ///
        title("`x'")
    graph save $DTA/`x'.gph, replace
}

* Clinical outcomes panel
graph combine ///
    $DTA/obesity.gph ///
    $DTA/high_pressure.gph ///
    $DTA/ISH.gph ///
    $DTA/diabetes.gph ///
    $DTA/high_chol.gph ///
    $DTA/high_TG.gph ///
    $DTA/low_HDL.gph, ///
    col(2) xsize(44) ysize(48)

graph save $DTA/clinical.gph, replace
graph use  $DTA/clinical.gph
graph export $DTA/clinical.png, as(png) width(1100) height(1200) replace

* Subclinical outcomes panel
graph combine ///
    $DTA/BMI_ln.gph ///
    $DTA/sbp_ln.gph ///
    $DTA/dbp_ln.gph ///
    $DTA/PP_ln.gph ///
    $DTA/glucose_ln.gph ///
    $DTA/cholestrol_ln.gph ///
    $DTA/Triglyceride_ln.gph ///
    $DTA/chol_HDL_ln.gph, ///
    col(2) xsize(44) ysize(48)

graph save $DTA/subclinical.gph, replace
graph use  $DTA/subclinical.gph
graph export $DTA/subclinical.png, as(png) width(1100) height(1200) replace


*---------------------------------------------------------------*
* Table 3. Reduced-Form Equations for Clinical/Subclinical Outcomes
*---------------------------------------------------------------*

clear
set more off

macro define DTA "LLE"
use $DTA/, replace

* Outcomes for reduced-form regressions
local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    chol_HDL_ln Triglyceride_ln

local v = 0

*---------------------------------------------------------------*
* 1. Loop over outcomes: store slope changes & APE
*---------------------------------------------------------------*

foreach x of local outcomes {

    local v = `v' + 1

    * Slope change among 65+ (dum == 1)
    quietly reg `x' i.cut#c.year_c i.cut c.year_c ///
        if dum==1, robust

    scalar coef_1_`v' = _b[1.cut#c.year_c]
    scalar se_1_`v'   = _se[1.cut#c.year_c]

    * Slope change among 40-64 (dum == 0)
    quietly reg `x' i.cut#c.year_c i.cut c.year_c ///
        if dum==0, robust

    scalar coef_2_`v' = _b[1.cut#c.year_c]
    scalar se_2_`v'   = _se[1.cut#c.year_c]

    * Difference-in-differences in slope change (RKD + DiD)
    quietly reg `x' i.cut#c.year_c#i.dum i.cut c.year_c ///
        i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum, robust

    scalar coef_3_`v' = _b[1.cut#c.year_c#1.dum]
    scalar se_3_`v'   = _se[1.cut#c.year_c#1.dum]

    * Average of integrated slope
    quietly margins, expression(2 * _b[1.cut#c.year_c#1.dum])

    scalar coef_4_`v' = r(table)[1, 1]
    scalar se_4_`v'   = r(table)[2, 1]
}

*---------------------------------------------------------------*
* 2. Collect scalars into a dataset
*---------------------------------------------------------------*

clear
set obs 1

local n_outcomes : word count `outcomes'

forvalues x = 1/`n_outcomes' {
    forvalues i = 1/4 {
        gen coef_`i'_`x' = coef_`i'_`x'
        gen se_`i'_`x'   = se_`i'_`x'
    }
}

gen id = 1

reshape long coef_1_ se_1_ ///
              coef_2_ se_2_ ///
              coef_3_ se_3_ ///
              coef_4_ se_4_, i(id)

forvalues x = 1/4 {
    gen test_`x' = coef_`x'_ / se_`x'_
    gen p_`x'    = 2*(1 - normprob(abs(test_`x')))
}

forvalues x = 1/4 {
    gen upper_`x' = coef_`x'_ + 1.96*se_`x'_
    gen lower_`x' = coef_`x'_ - 1.96*se_`x'_
}

save $DTA/, replace

*---------------------------------------------------------------*
* Figure 5: Main IV Estimates
*---------------------------------------------------------------*

use $DTA/, replace

* List of outcomes
local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln ///
    walking_ln sedentary_ln pa_aerobic ///
    alchol_freq_ln alchol_number_ln alchol_high_freq_ln ///
    smoking smoking_number_ln insurance med_serv HEI_ln

* Number of outcomes
local n_outcomes : word count `outcomes'

* Common covariates
local covars i.sex i.edu i.marital i.rural i.region ///
             c.height HP_hist HL_hist IHD_hist STR_hist DM_hist

local v = 0

*---------------------------------------------------------------*
* 1. Loop over outcomes: baseline & covariate-adjusted IV
*---------------------------------------------------------------*

foreach x of local outcomes {

    local v = `v' + 1

    * Baseline model
    quietly ivreg2 `x' (emp = i.cut#c.year_c#i.dum) ///
        i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum, ///
        robust

    scalar coef1_`v' = _b[emp]
    scalar se1_`v'   = _se[emp]

    * Covariate-adjusted model
    quietly ivreg2 `x' (emp = i.cut#c.year_c#i.dum) ///
        i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum ///
        `covars', robust

    scalar coef2_`v' = _b[emp]
    scalar se2_`v'   = _se[emp]
}

*---------------------------------------------------------------*
* 2. Collect scalars into a dataset
*---------------------------------------------------------------*

clear
set obs 1

forvalues x = 1/`n_outcomes' {
    forvalues i = 1/2 {
        gen coef`i'_`x' = coef`i'_`x'
        gen se`i'_`x'   = se`i'_`x'
    }
}

gen id = 1

reshape long coef1_ se1_ coef2_ se2_, i(id)

save $DTA/BE_main_1, replace

*---------------------------------------------------------------*
* 3. Main forest plot (all outcomes)
*---------------------------------------------------------------*

use $DTA/BE_main_1, replace

* Confidence intervals
forvalues x = 1/2 {
    gen upper_`x' = coef`x'_ + 1.96*se`x'_
    gen lower_`x' = coef`x'_ - 1.96*se`x'_
}

* Clean up any previous indexing vars
cap drop _j1
cap drop _j2

* _j is created by reshape; shift indices to create gaps between blocks
replace _j = _j + 1
replace _j = _j + 2 if _j >= 9
replace _j = _j + 2 if _j >= 19
replace _j = _j + 2 if _j >= 24
replace _j = _j + 2 if _j >= 31
replace _j = _j + 2 if _j >= 35

* Jitter for baseline vs adjusted
gen _j1 = _j - 0.1
gen _j2 = _j + 0.1

* Reference x-range for category boxes
gen x1 = -2
gen x2 =  2
gen x  =  0

* Category boundaries (in y / index space)
gen r1_l =  1.25
gen r1_u =  8.75

gen r2_l = 10.25
gen r2_u = 18.75

gen r3_l = 20.25
gen r3_u = 23.75

gen r4_l = 25.25
gen r4_u = 30.75

gen r5_l = 32.25
gen r5_u = 34.75

gen r6_l = 36.25
gen r6_u = 37.75

tw ///
    (rspike upper_1 lower_1 _j1, lcolor("34 168 132*0.6") horizontal) ///
    (scatter _j1 coef1_,  color("34 168 132") msize(small)) ///
    (rspike upper_2 lower_2 _j2, lcolor("65 68 135*0.6") horizontal) ///
    (scatter _j2 coef2_,  color("65 68 135") msize(small)) ///
    ///
    (rspike r1_l r1_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r1_l r1_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 1.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 8.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r2_l r2_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r2_l r2_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 10.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 18.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r1_l r1_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r2_l r2_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    ///
    (rspike r3_l r3_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r3_l r3_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 20.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 23.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r4_l r4_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r4_l r4_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 25.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 30.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r5_l r5_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r5_l r5_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 32.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 34.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r6_l r6_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r6_l r6_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 36.25, range(-2 2) lcolor(black) lwidth(vthin)) ///
    (function y = 37.75, range(-2 2) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r3_l r3_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r4_l r4_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r5_l r5_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r6_l r6_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    , ///
    xlabel(-2(0.5)2, labsize(small) nogrid format(%4.1f)) ///
    xsize(20) ysize(30) ///
    ylabel( ///
        0.75  " " ///
        2     "Obesity" ///
        3     "Hypertension" ///
        4     "ISH" ///
        5     "Diabetes" ///
        6     "Hypercholesterolemia" ///
        7     "Hypertriglyceridemia" ///
        8     "Hypo-HDL" ///
        11    "BMI" ///
        12    "SBP" ///
        13    "DBP" ///
        14    "Pulse pressure" ///
        15    "Glucose" ///
        16    "Cholestrol" ///
        17    "Triglyceride" ///
        18    "HDL-cholesterol" ///
        21    "Days walking" ///
        22    "Sedentary hours" ///
        23    "Moderate/vigorous" ///
        26    "Alcohol frequency" ///
        27    "Alcohol amount" ///
        28    "Binge drinking" ///
        29    "Smoking status" ///
        30    "Cigarettes per day" ///
        33    "Sponsored/private ins." ///
        34    "Health examination" ///
        37    "Total score", ///
        nogrid notick labsize(small)) ///
    ytitle("") ///
    legend(order(2 "Baseline" 4 "Covariates adj.") pos(1) col(1) ///
           size(small) colgap(1)) ///
    xtitle("Effect estimates", size(small)) ///
    yscale(reverse) xscale(noline) yscale(noline) ///
    plotregion(margin(zero)) ///
    text(0.75  0 "(a) Clinical",          size(small)) ///
    text(9.75  0 "(b) Subclinical",       size(small)) ///
    text(19.75 0 "(c) Physical activities", size(small)) ///
    text(24.75 0 "(d) Health behaviors",  size(small)) ///
    text(31.75 0 "(e) Medical access",    size(small)) ///
    text(35.75 0 "(f) Dietary patterns",  size(small))

graph export $DTA/.png, width(2000) height(3000) as(png) replace




