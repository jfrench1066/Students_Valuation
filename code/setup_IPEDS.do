/*

This do file creates a data set with the IPEDS variables used to construct Table 2.

by Esteban Aucejo, Jacob French, Basit Zafar

*/


*SAT/ACT  IPEDS 2015-2018
{
use "${data}ipeds/SAT-ACT/cdsfile_all_STATA_RV_5162020-864.dta", clear

*Verbal 25 percentile
gen old=satvr25+satwr25

merge m:1 old using "${data}ipeds/SAT-ACT/SAT_verbal_crosswalk.dta"

replace satvr25=new if _merge==3

drop if _merge==2
drop satwr25 _merge old new

*Verbal 75 percentile
gen old=satvr75+satwr75

merge m:1 old using "${data}ipeds/SAT-ACT/SAT_verbal_crosswalk.dta"

replace satvr75=new if _merge==3

drop if _merge==2
drop satwr75 _merge old new

*Math 25 percentile

gen old=satmt25

merge m:1 old using "${data}ipeds/SAT-ACT/SAT_math_crosswalk.dta"
drop if _merge==2

replace satmt25=new if _merge==3

drop old new _merge

*Math 75 percentile

gen old=satmt75

merge m:1 old using "${data}ipeds/SAT-ACT/SAT_math_crosswalk.dta"
drop if _merge==2

replace satmt75=new if _merge==3

drop old new _merge

tempfile sat_act_2015

save `sat_act_2015', replace

use  "${data}ipeds/SAT-ACT/cdsfile_all_STATA_RV_5162020-855.dta", clear

append using "${data}ipeds/SAT-ACT/cdsfile_all_STATA_RV_5162020-242.dta"

append using "${data}ipeds/SAT-ACT/cdsfile_all_STATA_RV_5162020-778.dta"

append using `sat_act_2015'

tempfile act_sat

save `act_sat', replace
}

********************************************************************************

use "${data}ipeds/cdsfile_all_STATA_RV_5152020-609.dta", clear //IPEDS 2018

merge 1:1 unitid using "${data}ipeds/cdsfile_all_STATA_RV_5152020-635.dta"

keep if  _merge==3

drop _merge

merge 1:m unitid using `act_sat'
keep if _merge==3

drop _merge

*Generating the variables of interest

gen female = eftotlw/eftotlt
gen black= efbkaat/ eftotlt
gen white= efwhitt/ eftotlt
gen hispanic= efhispt/ eftotlt
gen asian=efasiat/ eftotlt
gen native_american =efaiant/ eftotlt
gen pacific_islander= efnhpit/ eftotlt
gen morethanone= ef2mort/ eftotlt
gen unknown= efunknt / eftotlt
gen intl_student= efnralt/ eftotlt
gen weight_sat=satnum
gen weight_act=actnum
gen weight_all=eftotlt

#delimit ;
gen ASU=1 if 
unitid==448886|
unitid==420574|
unitid==104151|
unitid==407009
;
#delimit cr

keep if sector<=3 //keeping only 4-year colleges

gen US_1=1  //all other 4 year colleges

sort unitid

bys unitid: gen first=1 if _n==1
egen total_US_1=total(eftotlt) if first==1

egen total_ASU=total(eftotlt) if ASU==1 & first==1


bys fips: egen max=max( eftotlt) if first==1 & sector==1
gen US_2=1 if max== eftotlt & max!=.  //All flagship univiersity in the state

gen US_3=1 if first==1 & c18basic==15

egen total_US_2=total(eftotlt) if first==1 & US_2==1
egen total_US_3=total(eftotlt) if first==1 & US_3==1


gen top=.

replace top=1 if first==1 & unitid==186131 //Princeton University
replace top=1 if first==1 & unitid==166027 //Harvard University
replace top=1 if first==1 & unitid==190150 //Columbia University
replace top=1 if first==1 & unitid==166683 //MIT
replace top=1 if first==1 & unitid==130794 //Yale University
replace top=1 if first==1 & unitid==243744 //Standford University
replace top=1 if first==1 & unitid==144050 //University of Chicago
replace top=1 if first==1 & unitid==215062 //University of Pennsylvania
replace top=1 if first==1 & unitid==147767 //Northwestern University
replace top=1 if first==1 & unitid==198419 //Duke University
replace top=1 if first==1 & unitid==162928 //Johns Hopkins University

egen total_top=total(eftotlt) if first==1 & top==1

foreach var in satvr25 satvr75 satmt25 satmt75 actcm25 actcm75 weight_sat weight_act{

gen `var'_h =`var'

}


keep unitid female black white hispanic intl_student satvr25* satvr75* satmt25* satmt75* actcm25* actcm75* weight_all weight_sat* weight_act* US_* ASU total_* first top native_american pacific_islander morethanone unknown asian

save "${data}IPEDS.dta", replace

