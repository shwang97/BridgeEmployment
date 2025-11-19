/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the supplementary figures/tables
*/

clear
set more off
macro define DTA "LLE"

*---------------------------------------------------------------*
* Supplementary Figure 2. Density Graphs for Biomarkers
*---------------------------------------------------------------*
use $DTA/, replace

* Variables for density plots
local densvars ///
    BMI HE_sbp HE_dbp PP glucose cholestrol Triglyceride chol_HDL

*---------------------------------------------------------------*
* 1. Kernel density plots for each variable
*---------------------------------------------------------------*

foreach x of local densvars {
    kdensity `x', ///
        lcolor(navy) ///
        xsize(12) ysize(12) ///
        note("") ///
        xlabel(, nogrid) ///
        ylabel(, nogrid) ///
        title("`x'", size(medium)) ///
        xtitle("")

    graph save $DTA/kd_`x'.gph, replace
}

*---------------------------------------------------------------*
* 2. Combine all density graphs into one panel
*---------------------------------------------------------------*

graph combine ///
    $DTA/kd_BMI.gph          ///
    $DTA/kd_HE_sbp.gph       ///
    $DTA/kd_HE_dbp.gph       ///
    $DTA/kd_PP.gph           ///
    $DTA/kd_glucose.gph      ///
    $DTA/kd_cholestrol.gph   ///
    $DTA/kd_chol_HDL.gph     ///
    $DTA/kd_Triglyceride.gph, ///
    rows(3) cols(3) ///
    imargin(0 0 0 0) ///
    xsize(36) ysize(36)

graph save $DTA/kdensity.gph, replace
graph use  $DTA/kdensity.gph

graph export $DTA/kdensity.png, as(png) width(3600) height(3600) replace


*---------------------------------------------------------------*
* Supplementary Table 3. Comparison of Model Fit Across Polynomial Specifications
*   - BIC 
*   - AMSE from rdmse 
*---------------------------------------------------------------*

local outcomes ///
    emp obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln Triglyceride_ln ///
    chol_HDL_ln

local n_outcomes : word count `outcomes'

*---------------------------------------------------------------*
* 1. BIC comparison (linear vs quadratic in time)
*---------------------------------------------------------------*

use $DTA/, replace

foreach z in 0 1 {    // z = 0: dum==0; z = 1: dum==1

    local v = 0

    foreach x of local outcomes {

        local v = `v' + 1

        * Linear in time: c.year_c##i.cut
        quietly reg `x' c.year_c##i.cut if dum==`z', robust
        estimates store m1

        * Quadratic in time: c.year_c##c.year_c##i.cut
        quietly reg `x' c.year_c##c.year_c##i.cut if dum==`z', robust
        estimates store m2

        quietly estimates stats m1 m2

        forvalues i = 1/2 {
            scalar BIC`z'_`i'_`v' = r(S)[`i', 6]
        }
    }
}

clear
set obs 1

foreach z in 0 1 {
    forvalues v = 1/`n_outcomes' {
        forvalues i = 1/2 {
            gen BIC`z'_`i'_`v' = BIC`z'_`i'_`v'
        }
    }
}

gen id = 1

reshape long BIC0_1_ BIC1_1_ BIC0_2_ BIC1_2_, i(id)

save $DTA/, replace

*---------------------------------------------------------------*
* 2. AMSE comparison using rdmse 
*---------------------------------------------------------------*

use $DTA/, replace

foreach z in 0 1 {

    local v = 0

    foreach x of local outcomes {

        local v = `v' + 1

        * Local polynomial of order 1
        quietly rdmse `x' year_c if dum==`z', ///
            deriv(1) h(5) b(5) p(1)
        scalar AMSE`z'_1_`v' = r(amse_cl)

        * Local polynomial of order 2
        quietly rdmse `x' year_c if dum==`z', ///
            deriv(1) h(5) b(5) p(2)
        scalar AMSE`z'_2_`v' = r(amse_cl)
    }
}

clear
set obs 1

foreach z in 0 1 {
    forvalues v = 1/`n_outcomes' {
        forvalues i = 1/2 {
            gen AMSE`z'_`i'_`v' = AMSE`z'_`i'_`v'
        }
    }
}

gen id = 1

reshape long AMSE0_1_ AMSE1_1_ AMSE0_2_ AMSE1_2_, i(id)

save $DTA/, replace

*---------------------------------------------------------------*
* Supplementary Figure 3: Balance Check RKD vs RKD+DiD (Pre-treatment Covariates)
*---------------------------------------------------------------*

use $DTA/, replace

tab edu,    gen(edu)
tab region, gen(region)

local balvars ///
    sex ///
    edu1 edu2 edu3 ///
    marital rural ///
    region1 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 ///
    height_std ///
    HP_hist HL_hist IHD_hist STR_hist DM_hist

local n_balvars : word count `balvars'

local v = 0

foreach x of local balvars {

    local v = `v' + 1

    * Slope change among aged 65+
    quietly reg `x' i.cut#c.year_c i.cut c.year_c if dum==1, robust
    scalar coef1_`v' = _b[1.cut#c.year_c]
    scalar se1_`v'   = _se[1.cut#c.year_c]

    * Difference in slope change between 65+ and 40-64
    quietly reg `x' i.cut#c.year_c#i.dum ///
                  i.cut c.year_c i.cut#c.year_c ///
                  i.cut#i.dum c.year_c#i.dum i.dum, robust
    scalar coef2_`v' = _b[1.cut#c.year_c#1.dum]
    scalar se2_`v'   = _se[1.cut#c.year_c#1.dum]
}

clear
set obs 1

forvalues v = 1/`n_balvars' {
    forvalues i = 1/2 {
        gen coef`i'_`v' = coef`i'_`v'
        gen se`i'_`v'   = se`i'_`v'
    }
}

gen id = 1

reshape long coef1_ se1_ coef2_ se2_, i(id)

forvalues x = 1/2 {
    gen upper`x' = coef`x'_ + 1.96*se`x'_
    gen lower`x' = coef`x'_ - 1.96*se`x'_
}

save $DTA/, replace

* graph
use $DTA/, replace

replace _j = _j + 1
replace _j = _j + 2 if _j >= 8
replace _j = _j + 2 if _j >= 21
replace _j = _j + 2 if _j >= 24

gen _j1 = _j - 0.1
gen _j2 = _j + 0.1

graph set window fontface "Times New Roman"

gen r1_l =  1.25
gen r1_u =  7.75

gen r2_l =  9.25
gen r2_u = 20.75

gen r3_l = 22.25
gen r3_u = 23.75

gen r4_l = 25.25
gen r4_u = 30.75

gen x1 = -0.1
gen x2 =  0.1
gen x  =  0

tw ///
    (rspike upper1 lower1 _j1, ///
        lcolor("34 168 132*0.6") horizontal) ///
    (scatter _j1 coef1_, ///
        color("34 168 132")) ///
    (rspike upper2 lower2 _j2, ///
        lcolor("65 68 135*0.6") horizontal) ///
    (scatter _j2 coef2, ///
        color("65 68 135")) ///
    ///
    (rspike r1_l r1_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r1_l r1_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 1.25, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    (function y = 7.75, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r2_l r2_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r2_l r2_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 9.25,  range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    (function y = 20.75, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r3_l r3_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r3_l r3_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 22.25, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    (function y = 23.75, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r4_l r4_u x1, lcolor(black) lwidth(vthin)) ///
    (rspike r4_l r4_u x2, lcolor(black) lwidth(vthin)) ///
    (function y = 25.25, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    (function y = 30.75, range(-0.1 0.1) lcolor(black) lwidth(vthin)) ///
    ///
    (rspike r1_l r1_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r2_l r2_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r3_l r3_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    (rspike r4_l r4_u x, lcolor(black) lwidth(vthin) lpattern(dash)) ///
    , ///
    xlabel(-0.1(0.05)0.1, nogrid labsize(small)) ///
    xsize(15) ysize(20) ///
    ylabel( ///
        0.75 " " ///
        2    "Women" ///
        3    "Middle school or less" ///
        4    "High school" ///
        5    "College" ///
        6    "Married" ///
        7    "Rural" ///
        10   "Seoul" ///
        11   "Busan" ///
        12   "Daegu" ///
        13   "Incheon" ///
        14   "Daejeon" ///
        15   "Ulsan" ///
        16   "Gyeonggi" ///
        17   "Gangwon" ///
        18   "Choongcheong" ///
        19   "Jeonnam" ///
        20   "Gyungsang" ///
        23   "Height" ///
        26   "High pressure" ///
        27   "Hyperlipidemia" ///
        28   "Ischemic heart disease" ///
        29   "Stroke" ///
        30   "Diabetes", ///
        nogrid labsize(small) notick) ///
    ytitle("") ///
    legend(order(2 "RKD" 4 "RKD+DiD") pos(1) col(1) size(small)) ///
    xtitle("Effect estimates", size(small)) ///
    yscale(reverse) xscale(noline) yscale(noline) ///
    plotregion(margin(zero)) ///
    text(0.75  0 "(a) Sociodemographic",         size(small)) ///
    text(8.75  0 "(b) Region",                   size(small)) ///
    text(21.75 0 "(c) Physical",                 size(small)) ///
    text(24.75 0 "(c) Medical history (family)", size(small))

graph export $DTA/.png, as(png) width(1500) height(2000) replace

*---------------------------------------------------------------*
* Supplementary Table 4: Romano-Wolf correction
*---------------------------------------------------------------*

use $DTA/, replace

gen int_1 = cut*year_c*dum   
gen int_2 = cut*year
gen int_3 = cut*dum
gen int_4 = year_c*dum

local controls ///
    int_2 int_3 int_4 ///
    cut dum year_c ///
    sex2 edu2 edu3 ///
    region2 region3 region4 region5 region6 region7 region8 region9 ///
    region10 region11 region12 region13 ///
    height HP_hist HL_hist IHD_hist STR_hist DM_hist

local clinical ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL

rwolf2 ///
    (ivregress 2sls obesity      (emp = int_1) `controls', robust) ///
    (ivregress 2sls high_pressure (emp = int_1) `controls', robust) ///
    (ivregress 2sls ISH          (emp = int_1) `controls', robust) ///
    (ivregress 2sls diabetes     (emp = int_1) `controls', robust) ///
    (ivregress 2sls high_chol    (emp = int_1) `controls', robust) ///
    (ivregress 2sls high_TG      (emp = int_1) `controls', robust) ///
    (ivregress 2sls low_HDL      (emp = int_1) `controls', robust), ///
    indepvars(emp, emp, emp, emp, emp, emp, emp) ///
    verbose reps(1000)

rwolf2 ///
    (ivregress 2sls BMI_ln          (emp = int_1) `controls', robust) ///
    (ivregress 2sls HE_sbp          (emp = int_1) `controls', robust) ///
    (ivregress 2sls HE_dbp          (emp = int_1) `controls', robust) ///
    (ivregress 2sls PP_ln           (emp = int_1) `controls', robust) ///
    (ivregress 2sls glucose_ln      (emp = int_1) `controls', robust) ///
    (ivregress 2sls cholestrol_ln   (emp = int_1) `controls', robust) ///
    (ivregress 2sls Triglyceride_ln (emp = int_1) `controls', robust) ///
    (ivregress 2sls chol_HDL_ln     (emp = int_1) `controls', robust), ///
    indepvars(emp, emp, emp, emp, emp, emp, emp, emp) ///
    verbose reps(1000)


*---------------------------------------------------------------*
* Supplementary Table 5. Hausman tests
* 1) Comparing models with vs without excluding missing
*    pre-treatment characteristics
* 2) Comparing unweighted vs weighted results
*---------------------------------------------------------------*

*===============================================================*
* Part 1. Excluding vs NOT excluding pre-treatment missingness
*===============================================================*

use $DTA/, replace

* Interactions
gen int1 = cut*dum
gen int2 = year_c*dum 
gen int3 = cut*year_c 
gen int4 = cut*year_c*dum 

* Sample indicator: 1 if all pre-treatment covariates observed
gen sample = 1
foreach x in sex edu marital rural region height HP_hist HL_hist IHD_hist STR_hist DM_hist {
    replace sample = 0 if `x'==.
}

* Outcomes
local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln

tempname mem
tempfile tbl

postfile `mem' str25 outcome ///
              str18 b1 double se1 ///
              str18 b2 double se2 ///
              double pH ///
    using `tbl', replace

*--------------------------------------------------------*
* Loop over outcomes
*--------------------------------------------------------*
foreach x of local outcomes {

    * Model 1: restricted sample (complete pre-treatment covariates)
    quietly ivregress 2sls `x' (emp = int4) ///
        cut year_c dum int1 int2 int3 if sample==1, robust
    scalar bU_emp  = _b[emp]
    scalar seU_emp = _se[emp]
    scalar pU_emp  = r(table)[4, 1]

    * Model 2: full sample (no exclusion on pre-treatment covariates)
    quietly ivregress 2sls `x' (emp = int4) ///
        cut year_c dum int1 int2 int3, robust
    scalar bW_emp  = _b[emp]
    scalar seW_emp = _se[emp]
    scalar pW_emp  = r(table)[4, 1]

    * GMM-based Hausman test: restricted vs unrestricted sample
    quietly gmm  ///
        (eq1: `x' - {xb: emp cut year_c dum int1 int2 int3 _cons}) ///
        (eq2: sample*(`x' - {xc: emp cut year_c dum int1 int2 int3 _cons})), ///
        instruments(eq1: cut year_c dum int1 int2 int3 int4) ///
        instruments(eq2: cut year_c dum int1 int2 int3 int4) ///
        onestep winitial(unadjusted, indep) vce(robust, indep)

    test _b[xb:emp] = _b[xc:emp]
    scalar pH = r(p)

    local dagger = uchar(8224)
    local sU = cond(pU_emp<0.001,"***", ///
               cond(pU_emp<0.01,"**", ///
               cond(pU_emp<0.05,"*", ///
               cond(pU_emp<0.10,"`dagger'",""))))
    local sW = cond(pW_emp<0.001,"***", ///
               cond(pW_emp<0.01,"**", ///
               cond(pW_emp<0.05,"*", ///
               cond(pW_emp<0.10,"`dagger'",""))))
    local sH = cond(pH<0.001,"***", ///
               cond(pH<0.01,"**", ///
               cond(pH<0.05,"*", ///
               cond(pH<0.10,"`dagger'",""))))

    local b1s : display %9.3f bU_emp
    local b2s : display %9.3f bW_emp
    local b1s "`b1s'`sU'"
    local b2s "`b2s'`sW'"

    post `mem' ("`x'") ("`b1s'") (seU_emp) ("`b2s'") (seW_emp) (pH)
}

postclose `mem'
use `tbl', clear

order outcome b1 se1 b2 se2 pH
label var outcome "Outcome"
label var b1      "Model 1: b"
label var se1     "Model 1: se"
label var b2      "Model 2: b"
label var se2     "Model 2: se"
label var pH      "Hausman p-value"

format se1 %9.3f
format se2 %9.3f
format pH  %6.3f

list, noobs abbreviate(24)

export excel using "$DTA/.xlsx", firstrow(variables) replace


*===============================================================*
* Part 2. Unweighted vs weighted results
*===============================================================*

use $DTA/, replace

* Interactions
gen int1 = cut*dum
gen int2 = year_c*dum 
gen int3 = cut*year_c 
gen int4 = cut*year_c*dum 

* Weight and weighted regressors
gen sqrtw   = sqrt(wt_hs)
gen w_emp   = emp*sqrtw
gen w_int4  = int4*sqrtw
gen w_cut   = cut*sqrtw
gen w_year_c = year_c*sqrtw
gen w_dum   = dum*sqrtw
gen w_int1  = int1*sqrtw
gen w_int2  = int2*sqrtw
gen w_int3  = int3*sqrtw

local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln

tempname mem
tempfile tbl

postfile `mem' str25 outcome ///
              str18 b1 double se1 ///
              str18 b2 double se2 ///
              double pH ///
    using `tbl', replace

*--------------------------------------------------------*
* Loop over outcomes
*--------------------------------------------------------*
foreach x of local outcomes {

    * Model 1: unweighted
    quietly ivregress 2sls `x' (emp = int4) ///
        cut year_c dum int1 int2 int3, robust
    scalar bU_emp  = _b[emp]
    scalar seU_emp = _se[emp]
    scalar pU_emp  = r(table)[4, 1]

    * Model 2: weighted
    quietly ivregress 2sls `x' (emp = int4) ///
        cut year_c dum int1 int2 int3 [pw=wt_hs], robust
    scalar bW_emp  = _b[emp]
    scalar seW_emp = _se[emp]
    scalar pW_emp  = r(table)[4, 1]

    * Weighted outcome for GMM
    gen w_`x' = `x'*sqrtw

    * GMM-based Hausman test: unweighted vs weighted
    quietly gmm  ///
        (eq1: `x'   - {xb:  emp   cut year_c dum int1 int2 int3 _cons}) ///
        (eq2: w_`x' - {xc: w_emp w_cut w_year_c w_dum w_int1 w_int2 w_int3 sqrtw}), ///
        instruments(eq1: cut year_c dum int1 int2 int3 int4) ///
        instruments(eq2: w_cut w_year_c w_dum w_int1 w_int2 w_int3 sqrtw w_int4) ///
        onestep winitial(unadjusted, indep) vce(robust, indep)

    test _b[xb:emp] = _b[xc:w_emp]
    scalar pH = r(p)

    local dagger = uchar(8224)
    local sU = cond(pU_emp<0.001,"***", ///
               cond(pU_emp<0.01,"**", ///
               cond(pU_emp<0.05,"*", ///
               cond(pU_emp<0.10,"`dagger'",""))))
    local sW = cond(pW_emp<0.001,"***", ///
               cond(pW_emp<0.01,"**", ///
               cond(pW_emp<0.05,"*", ///
               cond(pW_emp<0.10,"`dagger'",""))))
    local sH = cond(pH<0.001,"***", ///
               cond(pH<0.01,"**", ///
               cond(pH<0.05,"*", ///
               cond(pH<0.10,"`dagger'",""))))

    local b1s : display %9.3f bU_emp
    local b2s : display %9.3f bW_emp
    local b1s "`b1s'`sU'"
    local b2s "`b2s'`sW'"

    post `mem' ("`x'") ("`b1s'") (seU_emp) ("`b2s'") (seW_emp) (pH)
}

postclose `mem'
use `tbl', clear

order outcome b1 se1 b2 se2 pH
label var outcome "Outcome"
label var b1      "Model 1: b"
label var se1     "Model 1: se"
label var b2      "Model 2: b"
label var se2     "Model 2: se"
label var pH      "Hausman p-value"

format se1 %9.3f
format se2 %9.3f
format pH  %6.3f

list, noobs abbreviate(24)

export excel using "$DTA/.xlsx", firstrow(variables) replace

*---------------------------------------------------------------*
* Supplementary Table 6: Results under Different Bandwidth Choices
*---------------------------------------------------------------*

clear
set more off
macro define DTA "LLE"
use $DTA/, replace

local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    chol_HDL_ln Triglyceride_ln

local n_outcomes : word count `outcomes'

local v = 0
foreach x of local outcomes {
    local v = `v' + 1

    forvalues i = 4/6 {
        quietly ivreg2 `x' ///
            (emp = i.cut#c.year_c#i.dum) ///
            i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum ///
            if year_c>=-`i' & year_c<=`i'-1

        scalar coef_`i'_`v' = _b[emp]
        scalar se_`i'_`v'   = _se[emp]
    }
}

clear
set obs 1

forvalues v = 1/`n_outcomes' {
    forvalues i = 4/6 {
        gen coef_`i'_`v' = coef_`i'_`v'
        gen se_`i'_`v'   = se_`i'_`v'
    }
}

gen id = 1

reshape long coef_4_ se_4_ coef_5_ se_5_ coef_6_ se_6_, i(id)

forvalues x = 4/6 {
    gen p_`x' = 2*(1 - normal(abs(coef_`x'/se_`x')))
}

save $DTA/, replace

*---------------------------------------------------------------*
* Supplementary Table 7. Results under Different Kernel Functions
* 1) Triangular kernel
* 2) Epanechnikov kernel
*---------------------------------------------------------------*

clear
set more off
macro define DTA "LLE"

local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln

local n_outcomes : word count `outcomes'


local reps = 500

* Kernels to compare
local kernels tri epa

*---------------------------------------------------------------*
* Loop over kernels: tri & epa
*---------------------------------------------------------------*

foreach kern of local kernels {

    *-----------------------------------------------------------*
    * 1. Bootstrap RKDIV estimates for each outcome
    *-----------------------------------------------------------*

    use $DTA/, replace

    * For each outcome, run bootstrap RKD-IV under kernel
    foreach x of local outcomes {

        forvalues i = 1/`reps' {
            preserve
            bsample

            * First-stage: kink in emp, difference 65+ vs 40-64
            rdrobust emp year_c if dum==1 & `x'!=., ///
                deriv(1) h(5) p(1) kernel(`kern')
            scalar coef1 = e(tau_cl)

            rdrobust emp year_c if dum==0 & `x'!=., ///
                deriv(1) h(5) p(1) kernel(`kern')
            scalar coef2 = e(tau_cl)

            scalar coef_first = coef1 - coef2

            * Reduced-form: kink in outcome, difference 65+ vs 40-64
            rdrobust `x' year_c if dum==1 & `x'!=., ///
                deriv(1) h(5) p(1) kernel(`kern')
            scalar coef3 = e(tau_cl)

            rdrobust `x' year_c if dum==0 & `x'!=., ///
                deriv(1) h(5) p(1) kernel(`kern')
            scalar coef4 = e(tau_cl)

            scalar coef_reduced = coef3 - coef4

            * Ratio: RF / FS = RKD-IV effect
            scalar result_`x'_`i' = coef_reduced / coef_first

            restore
        }
    }

    *-----------------------------------------------------------*
    * 2. Convert scalar results to a dataset (500 obs)
    *-----------------------------------------------------------*

    clear
    set obs `reps'

    foreach x of local outcomes {
        gen result_`x' = .
        forvalues i = 1/`reps' {
            replace result_`x' = result_`x'_`i' if _n==`i'
        }
    }

    save $DTA/kernel_`kern', replace
}

*---------------------------------------------------------------*
* Supplementary Figures 6 & 7.
* Effects Across Different Age Ranges in the Comparison Group
*---------------------------------------------------------------*

clear
set more off
macro define DTA "LLE"
use $DTA/, replace

* Outcomes and age thresholds
local outcomes ///
    obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
    BMI_ln sbp_ln dbp_ln PP_ln glucose_ln cholestrol_ln ///
    Triglyceride_ln chol_HDL_ln

local ages 40/55

*---------------------------------------------------------------*
* 1. Estimate effects by varying minimum age in comparison group
*---------------------------------------------------------------*

foreach x of local outcomes {
    forvalues a = `ages' {
        quietly ivreg2 `x' ///
            (emp = i.cut#c.year_c#i.dum) ///
            i.cut c.year_c i.cut#c.year_c i.cut#i.dum c.year_c#i.dum i.dum ///
            if age>=`a', robust

        scalar coef_`x'_`a' = _b[emp]
        scalar se_`x'_`a'   = _se[emp]
    }
}

*---------------------------------------------------------------*
* 2. Bring scalar results into a dataset
*---------------------------------------------------------------*

clear
set obs 1

foreach x of local outcomes {
    forvalues a = `ages' {
        gen coef_`x'_`a' = coef_`x'_`a'
        gen se_`x'_`a'   = se_`x'_`a'
    }
}

gen id = 1

reshape long ///
    coef_obesity_      se_obesity_      ///
    coef_high_pressure_ se_high_pressure_ ///
    coef_ISH_          se_ISH_          ///
    coef_diabetes_     se_diabetes_     ///
    coef_high_chol_    se_high_chol_    ///
    coef_high_TG_      se_high_TG_      ///
    coef_low_HDL_      se_low_HDL_      ///
    coef_BMI_ln_       se_BMI_ln_       ///
    coef_sbp_ln_       se_sbp_ln_       ///
    coef_dbp_ln_       se_dbp_ln_       ///
    coef_PP_ln_        se_PP_ln_        ///
    coef_glucose_ln_   se_glucose_ln_   ///
    coef_cholestrol_ln_ se_cholestrol_ln_ ///
    coef_chol_HDL_ln_  se_chol_HDL_ln_  ///
    coef_Triglyceride_ln_ se_Triglyceride_ln_, ///
    i(id)

save $DTA/BE_agerange, replace

*---------------------------------------------------------------*
* 3. Build upper / lower CIs and rename biomarkers for plotting
*---------------------------------------------------------------*

use $DTA/BE_agerange, replace

* 95% CI for each original outcome
foreach x of local outcomes {
    gen upper_`x' = coef_`x'_ + 1.96*se_`x'_
    gen lower_`x' = coef_`x'_ - 1.96*se_`x'_
}

* Create renamed biomarker variables (coef/upper/lower)
foreach y in coef upper lower {
    gen BMI_`y'          = `y'_BMI_ln
    gen SBP_`y'          = `y'_sbp_ln
    gen DBP_`y'          = `y'_dbp_ln
    gen PP_`y'           = `y'_PP_ln
    gen Glucose_`y'      = `y'_glucose_ln
    gen Cholesterol_`y'  = `y'_cholestrol_ln
    gen Triglyceride_`y' = `y'_Triglyceride_ln
    gen Hypo_HDL_`y'     = `y'_chol_HDL_ln
}


foreach x in obesity high_pressure ISH diabetes high_chol high_TG low_HDL ///
             BMI SBP DBP PP Glucose Cholesterol Triglyceride Hypo_HDL {

    tw ///
        (rarea `x'_upper `x'_lower _j, color(navy*0.3)) ///
        (line `x'_coef _j, lcolor(navy) lwidth(thick)) ///
        (function y=0, range(40 55) lpattern(dash) lcolor(black)), ///
        ylabel(-1(0.5)1, nogrid) ///
        xlabel(, nogrid) ///
        xsize(12) ysize(12) ///
        legend(off) ///
        xtitle("Age range") ///
        ytitle("Effect estimates") ///
        title("{bf:`x'}", size(medium))

    graph save $DTA/`x'_age.gph, replace
}

* Combine clinical outcomes
graph combine ///
    $DTA/obesity_age.gph      ///
    $DTA/high_pressure_age.gph ///
    $DTA/ISH_age.gph          ///
    $DTA/diabetes_age.gph     ///
    $DTA/high_chol_age.gph    ///
    $DTA/high_TG_age.gph      ///
    $DTA/low_HDL_age.gph,     ///
    rows(3) cols(3) ///
    imargin(0 0 0 0) ///
    xsize(30) ysize(32)

graph export $DTA/bw_age.png, as(png) width(3000) height(3200) replace

* Combine subclinical / biomarker outcomes
graph combine ///
    $DTA/BMI_age.gph          ///
    $DTA/SBP_age.gph          ///
    $DTA/DBP_age.gph          ///
    $DTA/PP_age.gph           ///
    $DTA/Glucose_age.gph      ///
    $DTA/Cholesterol_age.gph  ///
    $DTA/Triglyceride_age.gph ///
    $DTA/Hypo_HDL_age.gph,    ///
    rows(3) cols(3) ///
    imargin(0 0 0 0) ///
    xsize(30) ysize(32)

graph export $DTA/bw_age_sub.png, as(png) width(3000) height(3200) replace



