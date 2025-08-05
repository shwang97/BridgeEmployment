/* 
Replication Package:  Health Effects of Bridge Employment
Author: Sungsik Hwang, Shiro Furuya, and Jenna Nobles
Contact: shwang97@wisc.edu
Date: 2025-08-04
Note: The file contains the code for the Romano-Wolf correction
*/


clear
set more off
macro define DTA "LLE"
use $DTA/BE_clean, replace
keep if year>=2012 & year<=2021

gen int_1 = cut*year_c*dum 
gen int_2 = cut*year
gen int_3 = cut*dum 
gen int_4 = year_c*dum

foreach x in sex edu region {
	tab `x', g(`x')	
}

rwolf2 ///
(ivregress 2sls obesity (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls high_pressure (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls ISH (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls diabetes (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls high_chol (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls high_TG (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls low_HDL (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3  region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height HP_hist HL_hist IHD_hist STR_hist DM_hist, robust), ///
indepvars(emp, emp, emp, emp, emp, emp, emp) verbose reps(1000)

rwolf2 ///
(ivregress 2sls BMI_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c  sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls Midbp_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls PP_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c  sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls glucose_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c  sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls cholestrol_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c  sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls Triglyceride_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust) ///
(ivregress 2sls chol_HDL_ln (emp = int_1) int_2 int_3 int_4 cut dum year_c sex2 edu2 edu3 region2 region3 region4 region5 region6 region7 region8 region9 region10 region11 region12 region13 height  HP_hist HL_hist IHD_hist STR_hist DM_hist, robust), ///
indepvars(emp, emp, emp, emp, emp, emp, emp) verbose reps(1000)
