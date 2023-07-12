/*
	This file replicates all tables/figures in the CV JPUBE paper
	
	To perfectly match the paper, this file should be run on M1 macOS 
	operating system
	
	Esteban Aucejo, Jacob French, and Basit Zafar
	
	June 2023
*/

clear all 
set more off 

set type double 

graph set window fontface "Times New Roman"
set scheme s2color

*********************
* Paths
*********************

global homefolder "/Users/jfrench/Dropbox/ASU Survey/Analysis/CV_code/cv_replication/"

cd "${homefolder}"
global data "data/"
global code "code/"

global figures "output/figures/"	
global tables "output/tables/"

global restricted_data "/Users/jfrench/CV_restricted_data/"

*********************
* Options
*********************

local setup_validation			1
local setup_indiv_est			1
local shrunk_indiv_est			1

local table1 					1	
local table2					0	//Requires Restricted Data
local table3					1	
local table4					1	
local table5					1
local table6					1	

local figure1					1			
local figure2					1			
local figure3					1			
local figure4					0	//Requires Restricted Data		
local figure5					1			
local figure6					1			
local figure7a					1			
local figure7b					1			

local tableA1 					1	
local tableA2					1	
local tableA3					1
local tableA4					1	

local figureA1a					1			
local figureA1b					1			
local figureA2					0	//Requires Restricted Data		
local figureA3					1			

************************
*	individual ests
************************

if `setup_validation'{	
	set seed 0

	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)
	
	sort id _indiv_tag
	gen scen_n_drop1 = floor(((6 - 1)+1)*runiform() + 1) if _indiv_tag == 1
	by id: replace scen_n_drop1 = scen_n_drop1[1]
	
	gen _scen_n_drop2 = scen_n_drop1 - 1 + floor(((5 - 1)+1)*runiform() + 1) if _indiv_tag == 1
	gen scen_n_drop2 = mod(_scen_n_drop2,6)  + 1
	by id: replace scen_n_drop2 = scen_n_drop2[1]

	gen cost_drop = floor(((7 - 1)+1)*runiform() + 1) if _indiv_tag == 1
	by id: replace cost_drop = scen_n_drop1[1]

	preserve	
		keep if (hypo_num == scen_n_drop1 | hypo_num == scen_n_drop2) & cost == cost_drop
		keep id hypo_num cv scen_n_drop* cost_drop
		reshape wide cv, i(id) j(hypo_num)
		save ${data}validation_real.dta, replace
	restore
	drop if (hypo_num == scen_n_drop1 | hypo_num == scen_n_drop2) & cost == cost_drop

	*these variables hold wtp values for each factor
	local meth "frac"

	gen _wtp_no_covid_`meth' = .
	gen _wtp_in_person_`meth' = .
	gen _wtp_social_`meth' = .
	gen _wtp_vaccine_`meth' = .

	gen _wtpse_no_covid_`meth' = .
	gen _wtpse_in_person_`meth' = .
	gen _wtpse_social_`meth' = .
	gen _wtpse_vaccine_`meth' = .

	
	gen _coef_no_covid_`meth' = .
	gen _coef_no_covid_social_`meth' = .
	gen _coef_in_person_`meth' = .
	gen _coef_social_`meth' = .
	gen _coef_vaccine_`meth' = .
	gen _coef_cost_`meth' = .
	gen _coef__cons_`meth' = .
		
	gen _novarflag = 0 
	gen _mliter = 0
	set seed 1
	
	gen cv_prob = cv/100
		
	levelsof id, local(ids)
	foreach i in `ids'{
		*first check if any variation in likelihood
		quietly sum cv if id == `i'
		if r(max) != r(min){ 
			clear matrix
			quietly fracreg logit cv_prob no_covid in_person social vaccine no_covid_social cost_all if id == `i', iterate(1000) vce(robust)
			quietly replace _mliter = e(ic)

			foreach factor in no_covid in_person social vaccine cost no_covid_social _cons{			
				cap quietly replace _coef_`factor'_`meth' = _b[`factor'] if id == `i'
			}
			if _b[cost_all] != 0{
				quietly cap nlcom (wtp_in_person: - _b[in_person]/_b[cost_all]) ///
					(wtp_social: - (_b[social] + _b[no_covid_social])/_b[cost_all] ) ///
					(wtp_vaccine: - _b[vaccine]/_b[cost_all]) ///
					(wtp_no_covid: - (_b[no_covid]+ _b[no_covid_social])/_b[cost_all]) ///
					, post iterate(10000)
				if _rc==0{ //if no error
					foreach factor in no_covid in_person social vaccine{			
						quietly replace _wtp_`factor'_`meth' = _b[wtp_`factor'] if id == `i'
						quietly replace _wtpse_`factor'_`meth' = _se[wtp_`factor'] if id == `i'
					} 
				}
				else{
					foreach factor in no_covid in_person social vaccine{
						quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
						quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
						quietly replace _novarflag = 1 if id == `i'
					}
				}
			}
			else{
				foreach factor in no_covid in_person social vaccine{
					quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
					quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
					quietly replace _novarflag = 1 if id == `i'
				}
			}

		}
		*if no variation, assign 0 wtp
		else{
			foreach factor in no_covid in_person social vaccine{
				quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
				quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
				quietly replace _novarflag = 1 if id == `i'
			}
		}
		display("finished: `i'")
	}	
	
	bys ResponseId: gen _tag=1 if _n==1

	
	
	*keep only wtp variables
	keep if _tag==1
	keep Q81 _wtp_* _wtpse_* _coef_* id _novarflag cost2019 check1 check2 female honors low_inc first_gen nonwhite female ResponseId
	
	ren _* *
	
	*save data
	save ${data}indiv_WTP_val.dta, replace
}

if `setup_indiv_est'{
		
	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)
	
	sort ResponseId scenario
	
	//keep Q81 ResponseId id cv no_covid in_person social vaccine no_covid_social cost_all cost2019 check1 check2 female honors low_inc first_gen nonwhite female
	
	*these variables hold wtp values for each factor
	local meth "frac"

	gen _wtp_no_covid_`meth' = .
	gen _wtp_in_person_`meth' = .
	gen _wtp_social_`meth' = .
	gen _wtp_vaccine_`meth' = .

	gen _wtpse_no_covid_`meth' = .
	gen _wtpse_in_person_`meth' = .
	gen _wtpse_social_`meth' = .
	gen _wtpse_vaccine_`meth' = .

	
	gen _coef_no_covid_`meth' = .
	gen _coef_in_person_`meth' = .
	gen _coef_social_`meth' = .
	gen _coef_vaccine_`meth' = .
	gen _coef_cost_`meth' = .
	gen _coef__cons_`meth' = .
		
	gen _novarflag = 0 
	gen _mliter = 0
	set seed 1
	
	gen cv_prob = cv/100
	
	levelsof id, local(ids)
	foreach i in `ids'{
		*first check if any variation in likelihood
		quietly sum cv if id == `i'
		if r(max) != r(min){ 
			clear matrix
			quietly fracreg logit cv_prob no_covid in_person social vaccine no_covid_social cost_all if id == `i', iterate(1000) vce(robust)
			quietly replace _mliter = e(ic)

			foreach factor in no_covid in_person social vaccine cost no_covid_social _cons{			
				cap quietly replace _coef_`factor'_`meth' = _b[`factor'] if id == `i'
			}
			if _b[cost_all] != 0{
				quietly cap nlcom (wtp_in_person: - _b[in_person]/_b[cost_all]) ///
					(wtp_social: - (_b[social] + _b[no_covid_social])/_b[cost_all] ) ///
					(wtp_vaccine: - _b[vaccine]/_b[cost_all]) ///
					(wtp_no_covid: - (_b[no_covid]+ _b[no_covid_social])/_b[cost_all]) ///
					, post iterate(10000)
				if _rc==0{ //if no error
					foreach factor in no_covid in_person social vaccine{			
						quietly replace _wtp_`factor'_`meth' = _b[wtp_`factor'] if id == `i'
						quietly replace _wtpse_`factor'_`meth' = _se[wtp_`factor'] if id == `i'
					} 
				}
				else{
					foreach factor in no_covid in_person social vaccine{
						quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
						quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
						quietly replace _novarflag = 1 if id == `i'
					}
				}
			}
			else{
				foreach factor in no_covid in_person social vaccine{
					quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
					quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
					quietly replace _novarflag = 1 if id == `i'
				}
			}

		}
		*if no variation, assign 0 wtp
		else{
			foreach factor in no_covid in_person social vaccine{
				quietly replace _wtp_`factor'_`meth' = 0 if id == `i'
				quietly replace _wtpse_`factor'_`meth' = 0 if id == `i'
				quietly replace _novarflag = 1 if id == `i'
			}
		}
		display("finished: `i'")
	}	
	
	bys ResponseId: gen _tag=1 if _n==1
	
	*keep only wtp variables
	keep if _tag==1
	keep Q81 _wtp_* _wtpse_* _coef_* id _novarflag cost2019 check1 check2 female honors low_inc first_gen nonwhite female ResponseId
	
	ren _* *
	
	*save data
	save ${data}indiv_WTP.dta, replace
	
}

if `shrunk_indiv_est'{	
	use ${data}indiv_WTP.dta, clear

	gen excludeflag_social = 0
	gen excludeflag_in_person = 0
	gen excludeflag_vaccine = 0
	gen excludeflag_no_covid = 0
	
	foreach f in no_covid in_person vaccine social{	
		
		*exclude from shrink variance estimate if non-variance in individual estimation
		replace excludeflag_`f' = 1 if novarflag == 1

		*if really big SE, set estiamte to 0, exclude from shrink variance estimate
		_pctile wtpse_`f'_frac, p(95)
		local r3 = r(r1)
		replace excludeflag_`f' = 1 if wtpse_`f'_frac>`r3' & !missing(wtp_`f'_frac)
		replace wtpse_`f'_frac = 0 if wtpse_`f'_frac>`r3' & !missing(wtp_`f'_frac)
		replace wtp_`f'_frac = 0 if wtpse_`f'_frac>`r3' & !missing(wtp_`f'_frac)
		
		*if WTP within 3% of extremes, exclude from shrink variance estimate 
		_pctile wtp_`f'_frac if excludeflag_`f'!=1, p(3 97)
		local r1 = r(r1)
		local r2 = r(r2)
		
		*if WTP within 1% of extremes, exclude from shrink variance estimate and winsorize
		_pctile wtp_`f'_frac if excludeflag_`f'!=1, p(1 99)
		local r3 = r(r1)
		local r4 = r(r2)
		
		*if really be WTP, exclude from var calc
		replace excludeflag_`f' = 1 if wtp_`f'_frac<`r1' & !missing(wtp_`f'_frac)
		
		replace wtp_`f'_frac = `r3' if wtp_`f'_frac<`r3' & !missing(wtp_`f'_frac)
		
		replace excludeflag_`f' = 1 if wtp_`f'_frac>`r2' & !missing(wtp_`f'_frac)
	
		replace wtp_`f'_frac = `r4' if wtp_`f'_frac>`r4' & !missing(wtp_`f'_frac)
				
		*error variance
		cap drop _evar
		gen _evar = wtpse_`f'_frac * wtpse_`f'_frac 
		
		sum _evar if excludeflag_`f' != 1
		local Verror = r(sum)/r(N)
		
		sum wtp_`f'_frac if excludeflag_`f' != 1
		local Vtotal = r(Var)
		local Vtrue = `Vtotal' - `Verror'
		display "VTrue:  `Vtrue'"
		local wtpbar = r(mean)
			
		gen double _adjustvalue_`f' = (`Vtrue' / (`Vtrue' + _evar))

		gen wtp_`f'_s = `wtpbar' + _adjustvalue_`f' * (wtp_`f'_frac - `wtpbar')
	}
	
	save ${data}indiv_WTP_with_shrink.dta, replace
}

************************
*	tables
************************

if `table1'{
	use ${data}surveydata.dta, clear

	keep if !missing(cv) & !missing(cost_all)
	
	bys id hypo_num: gen _tag = 1 if _n==1
	
	gen cv_g100 = cv<100 if !missing(cv)
	
	forval s = 1/6{
		sum prob_scenario if _tag == 1 & hypo_num == `s' & !missing(cost_all)
		local prob`s': display %3.1f r(mean)
		
		sum cv if hypo_num == `s' & !missing(cost_all) & cost_1 == 1, d
		local prob_return`s': display %3.1f r(mean)
		local prob_return_med`s': display %3.1f r(p50)
		
		sum cv_g100 if hypo_num == `s' & !missing(cost_all) & cost_1 == 1 
		local prob_return100`s': display %4.3f r(mean)
	}
	
	cap file close table
	local fwt "file write table"
	file open table using "${tables}table1.tex", write replace
	`fwt' "\begin{table}[htp]  \caption{Scenarios}\label{t:scenario_likeli} \centering \begin{threeparttable} \begin{tabularx}{\textwidth}{l Y Y Y Y | Y Y Y}" _n
	`fwt' " & (1) & (2) & (3) & (4) & (5) & (6) & (7) \\"_n
	`fwt' " & & & & & \\"_n
	`fwt' " & COVID Controlled & In-Person Instruction & Normal Campus Social & Vaccine & Avg. Likelihood Return & Median Likelihood Return & Share Not Certain Return ($< 100$) \\[0.8em] \hline"_n
	`fwt' " & & & & & \\"_n
	`fwt' "	Scenario 1 & 1 & 1 & 1 & 0 & `prob_return1' & `prob_return_med1' & `prob_return1001' \\[0.8em]"
	`fwt' "	Scenario 2 & 0 & 0 & 0 & 0 & `prob_return2' & `prob_return_med2' & `prob_return1002' \\[0.8em]"
	`fwt' "	Scenario 3 & 1 & 0 & 0 & 0 & `prob_return3' & `prob_return_med3' & `prob_return1003' \\[0.8em]"
	`fwt' "	Scenario 4 & 1 & 0 & 1 & 0 & `prob_return4' & `prob_return_med4' & `prob_return1004' \\[0.8em]"
	`fwt' "	Scenario 5 & 0 & 1 & 1 & 0 & `prob_return5' & `prob_return_med5' & `prob_return1005' \\[0.8em]"
	`fwt' "	Scenario 6 & 0 & 1 & 1 & 1 & `prob_return6' & `prob_return_med6' & `prob_return1006' \\[0.8em]"
	
	`fwt' "\hline" _n
	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} Table characterizes each of the 6 scenarios participants were asked to consider, along with several statistics about the likelihood respondents assigned to returning for the fall semester with costs held constant at previous levels." ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table
	
}

if `table2'{
	*Setup IPEDS data
	do "${code}setup_IPEDS.do"
	
	*Survey data
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	drop _tag
	keep if !missing(cv) & !missing(cost_all)

	gen mid_inc= hhinc * 1000

	gen freshman=(yr_asu==1)
	gen sophomore=(yr_asu==2)
	gen junior=(yr_asu==3)
	gen senior=(yr_asu==4)

	*Creating weigths
	count 
	local wh=1 
	gen  weight_all=`wh' 

	count if act!=.
	local wh=1 
	gen  weight_act=`wh' 

	count if act!=. & honors==1
	local wh=1 
	gen  weight_act_h=`wh' 


	count if sat_verb!=. & sat_math!=.
	local wh=1 
	gen  weight_sat=`wh' 


	count if sat_verb!=. & sat_math!=. & honors==1
	local wh=1 
	gen  weight_sat_h=`wh' 
	 
	replace sat_math=. if sat_verb==. //only people that reported both for similarity with the IPEDS data

	*Creating %tile for ACT and SAT
	egen satvr25= pctile(sat_verb), p(25)
	egen satvr75= pctile(sat_verb), p(75)
	egen satmt25= pctile(sat_math), p(25) 
	egen satmt75= pctile(sat_math), p(75)
	egen actcm25= pctile(act), p(25) 
	egen actcm75= pctile(act), p(75) 

	egen satvr25_h= pctile(sat_verb) if honors==1, p(25)
	egen satvr75_h= pctile(sat_verb) if honors==1, p(75)
	egen satmt25_h= pctile(sat_math) if honors==1, p(25) 
	egen satmt75_h= pctile(sat_math) if honors==1, p(75)
	egen actcm25_h= pctile(act) if honors==1, p(25) 
	egen actcm75_h= pctile(act) if honors==1, p(75) 

	keep Res female weight_all weight_act* weight_sat* satvr25* satvr75* satmt25* satmt75* actcm25* actcm75* white black hispanic asian intl_student honors first_gen mid_inc freshman sophomore junior senior

	gen US_1=0
	gen US_2=0
	gen US_3=0

	gen top=0

	gen ASU=0

	gen ASU_2=0

	*Appending IPEDS data
	append using "${data}IPEDS.dta"
	
	*Appending confidential ASU data
	append using "${restricted_data}ASU_characteristics.dta"

	local US_2 "\begin{tabular}{c} Flagship \\ Univ.\tnote{d} \end{tabular}"

	replace mid_inc=mid_inc/1000

	*Writing the table
	local x = 2
	
	cap file close table
	file open table using "${tables}table2.tex", write replace
	local fwt "file write table"


	`fwt' "\begin{table}[htp]	\centering" _n
	`fwt' "\begin{threeparttable} \caption{Summary Statistics}\label{t:samptab}" _n
	`fwt' "\footnotesize"_n
			
	`fwt' "\renewcommand{\arraystretch}{1.5} \setlength\tabcolsep{4pt}" _n
	`fwt' "\begin{tabular}{lccccc|cccc}"_n
	`fwt' "  &   (1) & (2) & (3)& (4) & (5) & (6) & (7) & (8) & (9) \\"_n
	`fwt' " & & & & & & & & &  \\"_n

	
	`fwt' "  &   \begin{tabular}{c} Survey \\ All \end{tabular} & ASU & \begin{tabular}{c} P-value\\ (1)-(2) \end{tabular} & `US_`x'' &  \begin{tabular}{c} P-value\\ (1)-(4) \end{tabular} &  \begin{tabular}{c} Survey \\ Honors \end{tabular} &  \begin{tabular}{c} P-value\\ (6)-(2) \end{tabular}& \begin{tabular}{c} Top-10 \\ Univ.\tnote{e} \end{tabular} &  \begin{tabular}{c} P-value\\ (6)-(8) \end{tabular} \\"_n


	local female_n Female
	local black_n Black
	local white_n White
	local hispanic_n Hispanic 
	local asian_n Asian
	local intl_student_n "Int. Students"
	local satvr25_n "SAT Verbal 25th \%tile" 
	local satvr75_n "SAT Verbal 75th \%tile"
	local satmt25_n "SAT Math 25th \%tile"
	local satmt75_n "SAT Math 75th \%tile"
	local actcm25_n "ACT 25th \%tile"
	local actcm75_n "ACT 75th \%tile"

	local mid_inc_n "Family Income\tnote{a,c}"
	local first_gen_n "First-generation\tnote{a,b}"

	local freshman_n "Freshman\tnote{a}"
	local sophomore_n "Sophomore\tnote{a}"
	local junior_n "Junior\tnote{a}"
	local senior_n "Senior\tnote{a}"

	`fwt' "\hline " _n


	`fwt' "\multicolumn{6}{c|}{}&\multicolumn{4}{c}{}\\"_n


	foreach var in female black white hispanic intl_student {

	preserve

	sort US_1 unitid first

	duplicates drop unitid if US_1==1, force

	foreach g in ASU US_`x'{

	mean `var' [w= weight_all ], over(`g')
	local mean_`g'_0: display %4.2f (_b[`var'@0.`g'])
	local mean_`g'_1: display %4.2f (_b[`var'@1.`g'])

	test _b[`var'@0.`g']=_b[`var'@1.`g']
	local pvalue_`g': display %4.2f r(p)
	}

	drop if ASU==0 & honors!=1

	foreach g in ASU top{


	 mean `var' [w= weight_all ], over(`g')
	local mean_`g'_0_h: display %4.2f (_b[`var'@0.`g'])
	local mean_`g'_1_h: display %4.2f (_b[`var'@1.`g'])

	test _b[`var'@0.`g']=_b[`var'@1.`g']
	local pvalue_`g'_h: display %4.2f r(p)
	}

	restore


	`fwt' "``var'_n' & `mean_ASU_0' &  `mean_ASU_1' & `pvalue_ASU' &  `mean_US_`x'_1' & `pvalue_US_`x'' & `mean_ASU_0_h' & `pvalue_ASU_h' & `mean_top_1_h' & `pvalue_top_h' \\" _n


	}


	foreach var in  first_gen {

	mean `var' , over(ASU_2)
	local mean_ASU_2_0: display %4.2f (_b[`var'@0.ASU_2])
	local mean_ASU_2_1: display %4.2f (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2: display %4.2f r(p)

	preserve 
	drop if ASU_2==0 & honors!=1

	mean `var' , over(ASU_2)
	local mean_ASU_2_0_h: display %4.2f (_b[`var'@0.ASU_2])
	local mean_ASU_2_1_h: display %4.2f (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2_h: display %4.2f r(p)

	restore

	`fwt' "``var'_n' & `mean_ASU_2_0' &  `mean_ASU_2_1' & `pvalue_ASU_2' &  - & - & `mean_ASU_2_0_h' & `pvalue_ASU_2_h' & - & - \\" _n

	}

	foreach var in mid_inc  {

	mean `var' , over(ASU_2)
	local mean_ASU_2_0: display%11.2gc (_b[`var'@0.ASU_2])
	local mean_ASU_2_1: display %11.2gc (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2: display %4.2f r(p)

	preserve 
	drop if ASU_2==0 & honors!=1

	mean `var' , over(ASU_2)
	local mean_ASU_2_0_h: display %11.2gc (_b[`var'@0.ASU_2])
	local mean_ASU_2_1_h: display %11.2gc (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2_h: display %4.2f r(p)

	restore

	`fwt' "``var'_n' & `mean_ASU_2_0' &  `mean_ASU_2_1' & `pvalue_ASU_2' &  - & - & `mean_ASU_2_0_h' & `pvalue_ASU_2_h' & - & - \\" _n

	}

	`fwt' "\multicolumn{6}{c|}{}&\multicolumn{4}{c}{}\\"_n


	foreach var in freshman sophomore junior senior  {

	mean `var' , over(ASU_2)
	local mean_ASU_2_0: display %4.2f (_b[`var'@0.ASU_2])
	local mean_ASU_2_1: display %4.2f (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2: display %4.2f r(p)

	preserve 
	drop if ASU_2==0 & honors!=1

	mean `var' , over(ASU_2)
	local mean_ASU_2_0_h: display %4.2f (_b[`var'@0.ASU_2])
	local mean_ASU_2_1_h: display %4.2f (_b[`var'@1.ASU_2])

	test _b[`var'@0.ASU_2]=_b[`var'@1.ASU_2]
	local pvalue_ASU_2_h: display %4.2f r(p)

	restore

	`fwt' "``var'_n' & `mean_ASU_2_0' &  `mean_ASU_2_1' & `pvalue_ASU_2' &  - & - & `mean_ASU_2_0_h' & `pvalue_ASU_2_h' & - & - \\" _n

	}



	`fwt' "\multicolumn{6}{c|}{}&\multicolumn{4}{c}{}\\"_n

	preserve

	foreach var in satvr25 satvr75 satmt25 satmt75{

	foreach g in ASU US_`x'{


	 mean `var' [w= weight_sat ], over(`g')
	local mean_`g'_0: display %4.0f (_b[`var'@0.`g'])
	local mean_`g'_1: display %4.0f (_b[`var'@1.`g'])

	test _b[`var'@0.`g']=_b[`var'@1.`g']
	local pvalue_`g': display %4.2f r(p)
	}


	drop if ASU==0 & honors!=1

	foreach g in ASU top{

	 mean `var'_h  [w= weight_sat_h ], over(`g')
	local mean_`g'_0_h: display %4.0f (_b[`var'_h@0.`g'])
	local mean_`g'_1_h: display %4.0f (_b[`var'_h@1.`g'])

	test _b[`var'_h@0.`g']=_b[`var'_h@1.`g']
	local pvalue_`g'_h: display %4.2f r(p)

	}


	`fwt' "``var'_n' & `mean_ASU_0' &  `mean_ASU_1' & `pvalue_ASU' &  `mean_US_`x'_1' & `pvalue_US_`x'' & `mean_ASU_0_h' & `pvalue_ASU_h' & `mean_top_1_h' & `pvalue_top_h' \\" _n
	}

	restore

	`fwt' "\multicolumn{6}{c|}{}&\multicolumn{4}{c}{}\\" _n
	preserve
	foreach var in actcm25 actcm75{

	foreach g in ASU US_`x'{


	 mean `var' [w= weight_act ], over(`g')
	local mean_`g'_0: display %4.0f (_b[`var'@0.`g'])
	local mean_`g'_1: display %4.0f (_b[`var'@1.`g'])

	test _b[`var'@0.`g']=_b[`var'@1.`g']
	local pvalue_`g': display %4.2f r(p)
	}


	drop if ASU==0 & honors!=1

	foreach g in ASU top{


	 mean `var'_h [w= weight_act_h ], over(`g')
	local mean_`g'_0_h: display %4.0f (_b[`var'_h@0.`g'])
	local mean_`g'_1_h: display %4.0f (_b[`var'_h@1.`g'])

	test _b[`var'_h@0.`g']=_b[`var'_h@1.`g']
	local pvalue_`g'_h: display %4.2f r(p)
	}

	`fwt' "``var'_n' & `mean_ASU_0' &  `mean_ASU_1' & `pvalue_ASU' &  `mean_US_`x'_1' & `pvalue_US_`x'' & `mean_ASU_0_h' & `pvalue_ASU_h' & `mean_top_1_h' & `pvalue_top_h' \\" _n
	}

	restore

	sort US_1 unitid first

	duplicates drop unitid if US_1==1, force

	count if US_1==0

	local survey: display %11.0gc r(N)

	sum total_ASU 
	local ASU: display %11.0gc r(mean)

	sum total_US_`x'
	local US: display %20.0gc r(mean)

	sum total_top
	local top:  display %20.0gc r(mean)

	count if honors==1
	local honors: display %11.0gc r(N)

	`fwt' " & & & & & & & & &  \\"_n

	`fwt' "\hline" _n
	`fwt' "\multicolumn{6}{c|}{}&\multicolumn{4}{c}{}\\"_n

	`fwt' "\multicolumn{1}{c}{\textit{Sample Size}} &  `survey' &  `ASU' & & `US' & & `honors' & &`top' & \\" _n
	 
	`fwt' " & & & & & & & & &  \\"_n

	`fwt' "\hline" _n
	`fwt' "\end{tabular}" _n	

	`fwt' "\begin{tablenotes}[flushleft]"_n
	`fwt' "\item \textit{Notes:} Data in columns (2), (3) and (8) is from IPEDS 2018. The flagship universities are the 4-year public universities with the highest number of undergraduate students in each state. Means for these columns are weighted by total number of undergraduates in each institution. ACT and SAT data are weighted averages of 2018-2015 years from IPEDS. P-value columns show the p-value of a difference in means test between the two columns indicated by the numbers in the heading.  \item[a] Data in the ASU column from a different source. This data includes everyone taking at least one class for credit during the Spring semester of 2018 and attended ASU as their first full-time university. Income and First-generation variables for the ASU data are constructed with the data of the first available year, which is not the first year of college for most of the sample. \item[b] Students with no parent with a college degree. \item[c] Family income in thousands of dollars. \item[d] The largest public universities in each state. \item[e] Top 10 universities according to the US News Ranking 2020." _n
	`fwt' "\end{tablenotes}"_n 

	`fwt' "\end{threeparttable}" _n	
	`fwt' "\end{table}" _n

	 file close table
}

if `table3'{
	use ${data}surveydata.dta, clear
	set matsize 2000

	keep if !missing(cv) & !missing(cost_all)
	
	gen cv_prob = cv/100
	
	sum act_eq if _indiv_tag == 1, d
	gen low_act = act_eq < r(p50) if !missing(act_eq)
		
	gen group1 = 1
	
	gen group2 = 1 if low_inc==1
	gen group3 = 1 if low_inc==0

	gen group4 = 1 if first_gen==1
	gen group5 = 1 if first_gen==0

	gen group6 = 1 if nonwhite==1
	gen group7 = 1 if nonwhite==0
		
	gen group8 = 1 if female==1
	gen group9 = 1 if female==0
	
	gen group10 = 1 if honors==0
	gen group11 = 1 if honors==1
	
	gen group12 = 1 if cohort==1
	gen group13 = 1 if cohort==2
	gen group14 = 1 if cohort==3
	
	local lab1 "All"

	local lab2 "Lower-Income"
	local lab3 "Higher-Income"
	local lab4 "First-gen."
	local lab5 "Second-gen."
	local lab6 "Nonwhite"
	local lab7 "White"
	local lab8 "Females"
	local lab9 "Males"
	local lab11 "Honors"
	local lab10 "Non-honors"
	local lab12 "2023+ grads"
	local lab13 "2022 grads"
	local lab14 "2021 grads"
	
	local lab_wtp_ip "In-Person"
	local lab_wtp_s "Social"
	local lab_wtp_nc "No COVID"
	local lab_wtp_v "Vaccine"
	   
	cap drop _indiv_tag
	bys ResponseId: gen _indiv_tag=1 if _n==1
	
	forval g=1/1{
		estimates clear

		quietly reg cost2019 if _indiv_tag==1 & group`g'==1
		local n`g': display %6.0fc e(N)
		local m`g': display %6.0fc _b[_cons]*1000
		local num_m`g' = _b[_cons]*1000
		
		estimates clear
		
		fracreg logit cv_prob no_covid in_person social vaccine cost_all no_covid_social ibn.id if group`g'==1, vce(robust) iterate(1000) nocons
		
		nlcom (wtp_ip: - _b[in_person]/_b[cost_all] * 1000) ///
			(wtp_s: - (_b[social] + _b[no_covid_social])/_b[cost_all] * 1000 ) ///
			(wtp_v: - _b[vaccine]/_b[cost_all] * 1000) ///
			(wtp_nc: - (_b[no_covid]+ _b[no_covid_social])/_b[cost_all] * 1000) ///
			, post
		
		foreach est in wtp_ip wtp_s wtp_nc wtp_v{
			local `est'`g': display %6.0fc _b[`est']
			quietly test `est'
			if r(p)<0.1{
				local `est'`g' = "``est'`g''*"
			}
			if r(p)<0.05{
				local `est'`g' = "``est'`g''*"
			}
			if r(p)<0.01{
				local `est'`g' = "``est'`g''*"
			}
			local `est'`g'_se: display %6.0fc _se[`est']
			local `est'`g'_se = subinstr("(``est'`g'_se')"," ","",.)
			local `est'`g'_pct: display %4.2f _b[`est']/`num_m1' * 100 
			local `est'`g'_pct = "[``est'`g'_pct'\%]"
		}
				
	}
	
	forval g=2(2)10{
		local g2 = `g' + 1
		estimates clear
		foreach i in `g' `g2'{
			
			quietly fracglm cv_prob no_covid in_person social vaccine cost_all no_covid_social ibn.id if group`i'==1, link(logit) vce(oim) iterate(1000) nocons
	
			estimates store est`i'a
			
			quietly reg cost2019 if _indiv_tag==1 & group`i'==1
			local n`i': display %6.0fc e(N)
			local m`i': display %6.0fc _b[_cons]
			local num_m`i' = _b[_cons]
		}
		
		quietly suest est`g'a est`g2'a, robust
		
		nlcom 	(wtp_ip`g': - _b[est`g'a_cv_prob:in_person]/_b[est`g'a_cv_prob:cost_all] * 1000 ) ///
			(wtp_s`g': - (_b[est`g'a_cv_prob:social] + _b[est`g'a_cv_prob:no_covid_social])/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_v`g': - _b[est`g'a_cv_prob:vaccine]/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_nc`g': - (_b[est`g'a_cv_prob:no_covid]+_b[est`g'a_cv_prob:no_covid_social])/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_ip`g2': - _b[est`g2'a_cv_prob:in_person]/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_s`g2': - (_b[est`g2'a_cv_prob:social]+_b[est`g2'a_cv_prob:no_covid_social])/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_v`g2': - _b[est`g2'a_cv_prob:vaccine]/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_nc`g2': - (_b[est`g2'a_cv_prob:no_covid]+_b[est`g2'a_cv_prob:no_covid_social])/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			, post
		
		foreach est in wtp_ip wtp_s wtp_nc wtp_v{
			forval i = `g'/`g2'{
				local `est'`i': display %6.0fc _b[`est'`i']
				quietly test `est'`i'
				if r(p)<0.1{
					local `est'`i' = "``est'`i''*"
				}
				if r(p)<0.05{
					local `est'`i' = "``est'`i''*"
				}
				if r(p)<0.01{
					local `est'`i' = "``est'`i''*"
				}
				local `est'`i'_se: display %6.0fc _se[`est'`i']
				local `est'`i'_se = subinstr("(``est'`i'_se')"," ","",.)
				local `est'`i'_pct: display %4.2f _b[`est'`i']/`num_m1' * 100 
				local `est'`i'_pct = "[``est'`i'_pct'\%]"
			}
						
			test _b[`est'`g']==_b[`est'`g2']
			local `est'`g'_p: display %4.3f r(p)
		}
	}
	
	forval g=12/12{
		local g2 = `g' + 1
		local g3 = `g' + 2
		estimates clear
		foreach i in `g' `g2' `g3'{
	
			quietly fracglm cv_prob no_covid in_person social vaccine cost_all no_covid_social ibn.id if group`i'==1, link(logit) vce(oim) iterate(1000) nocons
			estimates store est`i'a
			
			quietly reg cost2019 if _indiv_tag==1 & group`i'==1
			local n`i': display %6.0fc e(N)
			local m`i': display %6.0fc _b[_cons]
			local num_m`i' = _b[_cons]
		}
		
		quietly suest est`g'a est`g2'a est`g3'a, robust
		
		nlcom 	(wtp_ip`g': - _b[est`g'a_cv_prob:in_person]/_b[est`g'a_cv_prob:cost_all] * 1000 ) ///
			(wtp_s`g': - (_b[est`g'a_cv_prob:social] + _b[est`g'a_cv_prob:no_covid_social])/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_v`g': - _b[est`g'a_cv_prob:vaccine]/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_nc`g': - (_b[est`g'a_cv_prob:no_covid]+_b[est`g'a_cv_prob:no_covid_social])/_b[est`g'a_cv_prob:cost_all] * 1000) ///
			(wtp_ip`g2': - _b[est`g2'a_cv_prob:in_person]/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_s`g2': - (_b[est`g2'a_cv_prob:social]+_b[est`g2'a_cv_prob:no_covid_social])/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_v`g2': - _b[est`g2'a_cv_prob:vaccine]/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_nc`g2': - (_b[est`g2'a_cv_prob:no_covid]+_b[est`g2'a_cv_prob:no_covid_social])/_b[est`g2'a_cv_prob:cost_all] * 1000) ///
			(wtp_ip`g3': - _b[est`g3'a_cv_prob:in_person]/_b[est`g3'a_cv_prob:cost_all] * 1000) ///
			(wtp_s`g3': - (_b[est`g3'a_cv_prob:social]+_b[est`g3'a_cv_prob:no_covid_social])/_b[est`g3'a_cv_prob:cost_all] * 1000) ///
			(wtp_v`g3': - _b[est`g3'a_cv_prob:vaccine]/_b[est`g3'a_cv_prob:cost_all] * 1000) ///
			(wtp_nc`g3': - (_b[est`g3'a_cv_prob:no_covid]+_b[est`g3'a_cv_prob:no_covid_social])/_b[est`g3'a_cv_prob:cost_all] * 1000) ///
			, post
		
		foreach est in wtp_ip wtp_s wtp_nc wtp_v{
			foreach i in `g' `g2' `g3'{
				local `est'`i': display %6.0fc _b[`est'`i']
				quietly test `est'`i'
				if r(p)<0.1{
					local `est'`i' = "``est'`i''*"
				}
				if r(p)<0.05{
					local `est'`i' = "``est'`i''*"
				}
				if r(p)<0.01{
					local `est'`i' = "``est'`i''*"
				}
				local `est'`i'_se: display %6.0fc _se[`est'`i']
				local `est'`i'_se = subinstr("(``est'`i'_se')"," ","",.)
				local `est'`i'_pct: display %4.2f _b[`est'`i']/`num_m1' * 100 
				local `est'`i'_pct = "[``est'`i'_pct'\%]"
			}
						
			test _b[`est'`g']==_b[`est'`g2']==_b[`est'`g3']
			local `est'`g'_p: display %4.3f r(p)
		}
	}
	


	local title " WTP Estimates"


	********************
	*Table
	********************
	cap file close table
	local fwt "file write table"
	file open table using "${tables}table3.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:`name'} \setlength\tabcolsep{3pt} \centering \footnotesize \begin{threeparttable} \begin{tabularx}{0.7\textwidth}{lY@{\hskip 3.5em}YY@{\hskip 3.5em}YY} " _n
	`fwt' " &  & (1) & & (2) & \\ "_n
	`fwt' " & \\"_n

	`fwt' " & N & WTP Social & p-value & WTP In-Person & p-value \\ \hline"_n
	`fwt' " & \\"_n
	`fwt' " \multirow{3}{*}{All} & \multirow{3}{*}{`n1'} "
	foreach var in wtp_s wtp_ip{
		`fwt' "& ``var'1' & "
	}
	`fwt' " \\"_n
	`fwt' " & "
	foreach var in wtp_s wtp_ip{
		`fwt' "& ``var'1_se' & "
	}
	`fwt' " \\"_n
	`fwt' " & "
	foreach var in wtp_s wtp_ip{
		`fwt' "& ``var'1_pct' & "
	}
	`fwt' " \\"_n
	`fwt' " \\"_n
	`fwt' " \\"_n
	
	forval g=2(2)10{
		local g2 = `g' + 1
		`fwt' " \multirow{3}{*}{`lab`g''} & \multirow{3}{*}{`n`g''} "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'' & \multirow{6}{*}{``var'`g'_p'} "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'_se' & "
		}
		`fwt' " \\"_n	
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'_pct' & "
		}
		`fwt' " \\[0.5em]"_n	
		
		`fwt' " \multirow{3}{*}{`lab`g2''} & \multirow{3}{*}{`n`g2''}"
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'_se' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'_pct' & "
		}
		`fwt' " \\"_n
		`fwt' " \\"_n	
		`fwt' " \\"_n	
	}
	
	forval g=12/12{
		local g2 = `g' + 1
		local g3 = `g' + 2
		`fwt' " \multirow{3}{*}{`lab`g''} & \multirow{3}{*}{`n`g''} "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'' & \multirow{9}{*}{``var'`g'_p'} "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'_se' & "
		}
		`fwt' " \\"_n	
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g'_pct' & "
		}
		`fwt' " \\[0.5em]"_n	
		
		`fwt' " \multirow{3}{*}{`lab`g2''} & \multirow{3}{*}{`n`g2''}"
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'_se' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g2'_pct' & "
		}
		`fwt' " \\[0.5em]"_n
		
		`fwt' " \multirow{3}{*}{`lab`g3''} & \multirow{3}{*}{`n`g3''}"
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g3'' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g3'_se' & "
		}
		`fwt' " \\"_n
		
		`fwt' " & "
		foreach var in wtp_s wtp_ip{
			`fwt' "& ``var'`g3'_pct' & "
		}
		`fwt' " \\"_n
		
		`fwt' " \\"_n	
		`fwt' " \\"_n	
	}
	
	
	`fwt' " \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} Willingness-to-pay reported in dollars per year.
			Standard errors in parentheses derived via delta method.
			Willingness-to-Pay as a percent of average net costs (\textdollar`m1') displayed in brackets.
			P-value from test of equal WTP between demographic groups.
			Grads refers to student expected cohort by graduation year.
			*, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. " ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table
}

if `table4'{
	use ${data}surveydata.dta, clear
	
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace

	use ${data}indiv_WTP_with_shrink.dta, clear
	
	gen W_wtp_social = wtp_social_s * 1000
	gen W_wtp_inperson = wtp_in_person_s * 1000
	
	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
		
	gen group1 = 1
	
	gen group2 = 1 if low_inc==1
	gen group3 = 1 if low_inc==0

	gen group4 = 1 if first_gen==1
	gen group5 = 1 if first_gen==0

	gen group6 = 1 if nonwhite==1
	gen group7 = 1 if nonwhite==0
		
	gen group8 = 1 if female==1
	gen group9 = 1 if female==0
	
	gen group10 = 1 if honors==0
	gen group11 = 1 if honors==1
	
	gen group12 = 1 if cohort==1
	gen group13 = 1 if cohort==2
	gen group14 = 1 if cohort==3
	
	local lab1 "All"

	local lab2 "Lower-inc."
	local lab3 "Higher-inc."
	local lab4 "First-gen."
	local lab5 "Second-gen."
	local lab6 "Nonwhite"
	local lab7 "White"
	local lab8 "Females"
	local lab9 "Males"
	local lab11 "Honors"
	local lab10 "Non-honors"
	local lab12 "2023+ Grads"
	local lab13 "2022 Grads"
	local lab14 "2021 Grads"
	
	
	local abrev_W_wtp_social "s"
	local abrev_W_wtp_inperson "ip"
	
	sum cost2019 if group1 == 1
	local m_cost = r(mean)*1000
	local m_cost_format: display %6.0fc r(mean)*1000
	
	forval g=1/14{
		foreach wtp in W_wtp_social W_wtp_inperson{
			/*
			sum cost2019 if group`g' == 1
			local m_cost = r(mean)*1000
			*/
			
			sum `wtp' if group`g' == 1, d
			local wtp_`abrev_`wtp''_g`g'_p50: display %6.0fc r(p50)
			local wtp_`abrev_`wtp''_g`g'_p10: display %6.0fc r(p10)
			local wtp_`abrev_`wtp''_g`g'_p90: display %6.0fc r(p90)
			local wtp_`abrev_`wtp''_g`g'_sd: display %6.0fc r(sd)
			/*
			local wtp_`abrev_`wtp''_g`g'_m: display %6.0fc r(mean)
			local wtp_`abrev_`wtp''_g`g'_pct: display %4.2f r(mean)/`m_cost' * 100
			local wtp_`abrev_`wtp''_g`g'_pct = subinstr("[`wtp_`abrev_`wtp''_g`g'_pct'\%]"," ","",.)
			*/
			
			reg `wtp' if group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_m: display %6.0fc _b[_cons]
			local wtp_`abrev_`wtp''_g`g'_se: display %6.0fc _se[_cons]
			local wtp_`abrev_`wtp''_g`g'_se = subinstr("(`wtp_`abrev_`wtp''_g`g'_se')"," ","",.)
			local wtp_`abrev_`wtp''_g`g'_m = "`wtp_`abrev_`wtp''_g`g'_m'\,`wtp_`abrev_`wtp''_g`g'_se'"
					
			local wtp_`abrev_`wtp''_g`g'_pct: display %4.2f _b[_cons]/`m_cost' * 100
			local wtp_`abrev_`wtp''_g`g'_pct = subinstr("[`wtp_`abrev_`wtp''_g`g'_pct'\%]"," ","",.)

			test _cons
			local stars = ""
			if r(p) < 0.1{
				local stars = "`stars'*"
			}
			if r(p) < 0.05{
				local stars = "`stars'*"
			}
			if r(p) < 0.01{
				local stars = "`stars'*"
			}
			local wtp_`abrev_`wtp''_g`g'_m = "$\stackrel{`stars'}{`wtp_`abrev_`wtp''_g`g'_m'}$"			
			

			
			
			count if !missing(`wtp') & group`g' == 1
			local n = r(N)
			count if `wtp'>10 & !missing(`wtp') & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_g1: display %3.1f r(N)/`n' *100
			
			count if `wtp'<-10 & !missing(`wtp')  & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_l1: display %3.1f r(N)/`n' *100
			
			
		}
	}
	
	forval g = 2(2)10{
		local g2 = `g' + 1
		gen _group = 0 if group`g' == 1
		replace _group = 1 if group`g2' == 1
		
		foreach wtp in W_wtp_social W_wtp_inperson{
			ttest `wtp', by(_group)
			local p_`abrev_`wtp''_`g2': display %4.3f r(p) 
			local p_`abrev_`wtp''_`g2' = "p-val: `p_`abrev_`wtp''_`g2''"
		}
		cap drop _group
	}
	
	*anova for cohort (3 means, not 2)
	foreach wtp in W_wtp_social W_wtp_inperson{
		oneway `wtp' cohort
		local p_`abrev_`wtp''_14: display %4.3f Ftail(`r(df_m)',`r(df_r)',`r(F)')
		local p_`abrev_`wtp''_14 = "p-val: `p_`abrev_`wtp''_14'"
	} 
	
		
	*********
	*Table
	*********
	
	local title "Moments of Individual WTP Distribution"
	cap file close table
	local fwt "file write table"
	file open table using "${tables}table4.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:indiv_tab2} \setlength\tabcolsep{1pt} \centering \begin{threeparttable} \begin{tabularx}{1.4\textwidth}{lcYYYYYYYcYYYYYY} " _n
	`fwt' " & \multicolumn{7}{c}{WTP Social} && \multicolumn{7}{c}{WTP In-Person}\\ \cline{2-8} \cline{10-16}"_n
	`fwt' " & \\"_n
	`fwt' " & Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10 && Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10  \\ \hline"_n
	`fwt' " & \\"_n


	forval g=1/14{
		`fwt' " \qquad `lab`g'' "
		foreach wtp in s ip{
			foreach mom in m p50 sd p10 p90 g1 l1{
				`fwt' "& `wtp_`wtp'_g`g'_`mom''"
			}
			if "`wtp'" == "s"{
				`fwt' "& "
			}
		}
		`fwt' "\\"_n
		//`fwt' " & `wtp_s_g`g'_se' &&&&&&&& `wtp_ip_g`g'_se' &&&&&& \\"_n
		
		if "`g'"=="1"{
			`fwt' " & `wtp_s_g1_pct' &&&&&&&& `wtp_ip_g1_pct' &&&&&& \\"_n

		}
		if inlist(`g',3,5,7,9,11,14){
			`fwt' " & `p_s_`g'' &&&&&&&& `p_ip_`g'' &&&&&& \\"_n
		}
		
		if inlist(`g',1,3,5,7,9,11){
			`fwt' " & \\"_n
		}
	} 
	
	`fwt' " &  \\ \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} Bayesian shrinkage algorithm applied to estimates. 
			Willingness-to-Pay measured in dollars per year. 
			Standard error of mean estimate in parentheses. 
			P-value for equality of means between demographic groups displayed below each set of means.
			Willingness-to-Pay as a percent of average net costs (\textdollar`m_cost_format') displayed in brackets.
			Grads refers to student expected cohort by graduation year.
			Within means column  *, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. " ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table
	
}

if `table5'{
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace
	
	use ${data}indiv_WTP.dta, clear
	
	gen W_wtp_social = wtp_social_frac * 1000
	gen W_wtp_inperson = wtp_in_person_frac * 1000			
	
	*winsorize correctly consistent with shrinkage estimates 
	foreach f in inperson social{
		_pctile W_wtp_`f', p(3 97)
		local r1 = r(r1)
		local r2 = r(r2)
		sum W_wtp_`f' if W_wtp_`f' >=`r1' & W_wtp_`f'<=`r2'
		replace W_wtp_`f' = r(mean) if W_wtp_`f'<`r1' & !missing(W_wtp_`f')
		replace W_wtp_`f' = r(mean) if W_wtp_`f'>`r2' & !missing(W_wtp_`f')
	}
	
	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
	gen _know_act = !missing(act_eq)
	reg _know_act if !missing(W_wtp_social)
	local n_act_eq: display %3.1fc _b[_cons]*100
	display "`n_act_eq'"
	
 	*****************

	cap program drop regtab 
	program regtab
		// Varlist should be individual social WTP followed by In-Person WTP 
		syntax varlist [, q(integer -1) demos(string) econ_contr(string) noecon_contr(string) other_contr(string) name(string) title(string)]
		
		******************
		
		local lab_female "Female"
		local lab_low_inc "Lower-Income"
		local lab_nonwhite "Nonwhite"
		local lab_first_gen "First-gen."
		local lab_honors "Honors"
		local lab_fresh "Freshman"
		local lab_cohort1 "2023+ grads"
		local lab_cohort2 "2022 grads"
		local lab__cons "Constant"
		local lab_n "N"
		local lab_r2 "r$^{2}$"
		local lab_m "Mean"
				
		******************
		
		tokenize `varlist'
		local short_`1' "so"
		local short_`2' "ip"

		******************
		
		foreach var of varlist `econ_contr'{
			count if missing(`var')
			if r(N)>0{
				gen _miss_`var' = missing(`var')
				replace `var' = 0 if missing(`var')
				local miss_econ_contr = "`miss_econ_contr' _miss_`var'"
				local lab: var lab `var'
				lab var `var' "`lab'$^{a}$"
			}
		}
		
		foreach var of varlist `noecon_contr'{
			count if missing(`var')
			if r(N)>0{
				gen _miss_`var' = missing(`var')
				replace `var' = 0 if missing(`var')
				local miss_noecon_contr = "`miss_noecon_contr' _miss_`var'"
				local lab: var lab `var'
				lab var `var' "`lab'$^{a}$"
				sum _miss_`var' if !missing(W_wtp_social) & !missing(W_wtp_inperson)
				local n_`var': display %4.1f (1-r(mean)) * 100
			}
		}
		
		local indep1 "`demos'"
		local indep2 "`demos' `econ_contr'"
		local indep3 "`demos' `econ_contr' `noecon_contr' `miss_econ_contr' `miss_noecon_contr'"
		
		foreach var of varlist `econ_contr' `noecon_contr'{
			local lab_`var': var lab `var' 
		}
		
		
		*****************
		
		forval c = 1/3{
			foreach v in `1' `2'{
				if `q' == -1{
					reg `v' `indep`c'' `other_contr', robust
					local `short_`v''_r2_`c': display %4.3f e(r2)
				}
				else{
					qreg `v' `indep`c'', vce(robust) q(`q')
				}
				local `short_`v''_n_`c': display %6.0fc e(N)
				foreach i in `indep`c'' _cons{
					local b_`short_`v''_`i'_`c': display %6.0fc _b[`i']
					local b_`short_`v''_`i'_`c' = subinstr("`b_`short_`v''_`i'_`c''", " ","",.)
					
					local se_`short_`v''_`i'_`c': display %6.0fc _se[`i']
					local se_`short_`v''_`i'_`c' = subinstr("(`se_`short_`v''_`i'_`c'')", " ","",.)
					
					quietly test `i'
					if r(p)<0.1{
						local b_`short_`v''_`i'_`c' = "`b_`short_`v''_`i'_`c''*"
					}
					if r(p)<0.05{
						local b_`short_`v''_`i'_`c' = "`b_`short_`v''_`i'_`c''*"
					}
					if r(p)<0.01{
						local b_`short_`v''_`i'_`c' = "`b_`short_`v''_`i'_`c''*"
					}
				}
				sum `v' if e(sample) == 1
				local `short_`v''_m_`c': display %6.0fc r(mean)
			}
		}
		
		cap file close table
		local fwt "file write table"
		file open table using "${tables}`name'.tex", write replace
		`fwt' "\begin{table}[htp]\caption{`title'}\label{t:tab3} \setlength\tabcolsep{2pt} \centering \footnotesize \begin{threeparttable} \begin{tabularx}{\textwidth}{l@{\hskip 3.5em}YYY@{\hskip 3.5em}YYY} " _n
		`fwt' " & (1)  & (2)  & (3)  & (4)  & (5)  & (6) \\"_n
		`fwt' " & \multicolumn{3}{c}{WTP Social}  & \multicolumn{3}{c}{WTP In-Person} \\ \hline"_n
		
		`fwt' " \textbf{Demographic}& \\"_n
		foreach i in `demos'{
			`fwt' "\qquad `lab_`i'' "
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `b_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\ " _n
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `se_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\[0.5em] " _n
		}
		
		`fwt' " \textbf{Economic}& \\"_n
		foreach i in `econ_contr'{
			`fwt' "\qquad `lab_`i'' "
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `b_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\ " _n
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `se_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\[0.5em] " _n
		}
		
		`fwt' " \textbf{Non-Economic}& \\"_n
		foreach i in `noecon_contr'{
			`fwt' " \qquad `lab_`i'' "
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `b_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\ " _n
			foreach v in so ip{
				forval c = 1/3{
					`fwt' "& `se_`v'_`i'_`c'' "
				}
			}
			`fwt' "\\[0.5em] " _n
		}
		
		`fwt' " \textbf{Other Controls}& \\"_n
		`fwt' " \qquad College of Major "
		foreach v in so ip{
			forval c = 1/3{
				`fwt' "& Y "
			}
		}
		`fwt' "\\[0.5em] " _n
		
		`fwt' " &  \\ \hline" _n
		`fwt' " & \\"_n
		
		if `q' == -1{
			
			foreach i in n r2 m{
				`fwt' "\multicolumn{1}{r}{`lab_`i''} "
				foreach v in so ip{
					forval c = 1/3{
						`fwt' "& ``v'_`i'_`c'' "
					}
				}
				`fwt' "\\ " _n
			}
		}
		else{
			foreach i in n{
				`fwt' "\multicolumn{1}{r}{`lab_`i''} "
				foreach v in so ip{
					forval c = 1/3{
						`fwt' "& ``v'_`i'_`c'' "
					}
				}
				`fwt' "\\ " _n
			}
		}
		`fwt' " &  \\ \hline" _n

		`fwt' "\end{tabularx}"_n
		`fwt'"\begin{tablenotes}"_n
		#delimit ; 
		`fwt'"		
			\item[] \footnotesize \emph{Notes:} Estimates in dollars per year. 
				Heteroskedastic robust standard errors in parentheses. 
				Gross educ. expend. measures the total expenditure on a student's education from all sources including grants, scholarships, loans, family, etc.
				Grads refers to student expected cohort by graduation year. Dependent variable winsorized to mean above 97th and below 3rd percentiles.
				*, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. 
				
			\item[a] ACT scores known for `n_act_eq'\% of observations. Dummy for missing score also included.
		";
		#delimit cr
		`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
		file close table
	
	end

	******************
	
	local econ_contr "work work_hrs20 total_pay_an grants loans cost_above_instate "
	local noecon_contr "live_campus_c0 social_c0 study_hrs_c0 act_eq class_oh_an2"
	local other_contr "i.college_code"
	
	gen cohort1 = cohort==1
	gen cohort2 = cohort==2
	gen cohort3 = cohort==3
	
	local demos "low_inc first_gen nonwhite female honors cohort1 cohort2"
	regtab W_wtp_social W_wtp_inperson, demos(`demos') econ_contr(`econ_contr') noecon_contr(`noecon_contr') other_contr(`other_contr') title("Correlates of Individual-level WTPs: macOS numbers with ACT FE") name("table5_v3")	
	
}

if `table6'{
	******** ****************
	*setup IPEDS cross tabs
	*************************
	
	import excel using "${data}ipeds/ipeds_cross_tabs.xlsx", clear first
	ren Raceethnicity race
	replace race = substr(race, 1, strlen(race)-1) if substr(race, strlen(race), 1) ==" "
	ren E counts
	gen gender = 1 if Gender == "Women"
	replace gender = 0 if Gender == "Men"
	keep race gender counts
	keep if !missing(gender)
	drop if race=="Total"
	
	*combine asian and pacific islander
	replace race = "Asian or Pacific Islander" if race=="Asian"
	forval g=0/1{
		sum counts if race == "Native Hawaiian or Other Pacific Islander" & gender==`g'
		replace counts = counts + r(min) if race == "Asian or Pacific Islander" & gender==`g'
		
	}
	drop if race == "Native Hawaiian or Other Pacific Islander"
	
	total counts
	gen ipeds_shares = counts/_b[counts]
	drop counts
	
	save "${data}ipeds/ipeds_cross_tabs18.dta", replace
	
	use ${data}surveydata.dta, clear
	set matsize 2000

	keep if !missing(cv) & !missing(cost_all)
	
	gen race = "White" if Q27 == "White/Caucasian"
	replace race = "White" if Q27 == "White/Caucasian,Other"
	replace race = "Asian or Pacific Islander" if Q27 == "Asian/Pacific Islander"
	replace race = "Asian or Pacific Islander" if Q27 == "Asian/Pacific Islander,Other"
	replace race = "Hispanic or Latino" if Q27 == "Hispanic/Latino"
	replace race = "Hispanic or Latino" if Q27 == "Hispanic/Latino,Other"
	replace race = "American Indian or Alaska Native" if Q27 == "American Indian"
	replace race = "American Indian or Alaska Native" if Q27 == "American Indian,Other"
	replace race = "Black or African American" if Q27 == "Black/African American"
	replace race = "Black or African American" if Q27 == "Black/African American,Other"
	replace race = "Race/ethnicity unknown" if Q27 == "Other"
	replace race = "Nonresident alien" if Q24!="United States of America"
	replace race = "Two or more races" if missing(race)
	
	merge m:1 race gender using "${data}ipeds/ipeds_cross_tabs18.dta"
	drop _merge
	
	gen cv_prob = cv/100
	
	sum act_eq if _indiv_tag == 1, d
	gen low_act = act_eq < r(p50) if !missing(act_eq)
		
	gen group1 = 1
	gen group2 = prob_scenario>0
	gen group3 = 1
	gen group4 = check1 == 1 & check2 == 1
	_pctile Durationinseconds, p(4.0816, 95.9184)
	gen group5 = Durationinseconds>r(r1) & Durationinseconds<r(r2)
	gen group6 = 1
	
	*sample weights
	hashsort race gender
	by race gender: gen _tot_obs = _N
	
	gen w1 = 1
	gen w2 = 1
	gen w3 = prob_scenario
	gen w4 = 1
	gen w5 = 1
	gen w6 = ipeds_shares / _tot_obs * 100000 //100k just to rescale weights to be closer to 1
	
	local lab1 "Baseline"
	local lab2 "Subjective Prob. Scenario $>0$"
	local lab3 "Weighted by Subjective Prob."
	local lab4 "Pass Both Attention Checks"
	local lab5 "Survey Duration: 5th-95th"
	local lab6 "Weighted by Nat. Gender X Race Shares"
	
	local labwtp_ip "WTP In-Person"
	local labwtp_s "WTP Social"
	local labwtp_nc "WTP No COVID"
	local labwtp_v "WTP Vaccine"
	   
	local labm "\qquad Mean Cost"
	local labn "\qquad N"
	
	
	cap drop _indiv_tag
	bys ResponseId : gen _indiv_tag=1 if _n==1
	
	forval g=1/6{
		estimates clear

		estimates clear
		if `g' == 6{ //only reweight means for national sample weights
			quietly reg cost2019 if _indiv_tag==1 & group`g'==1 [iweight=w`g'], robust
		}
		else{
			quietly reg cost2019 if _indiv_tag==1 & group`g'==1, robust

		}
		
		local n`g': display %6.0fc e(N)
		local m`g': display %6.0fc _b[_cons]*1000
		local num_m`g' = _b[_cons]*1000
		
		estimates clear

		quietly fracreg logit cv_prob no_covid in_person social vaccine cost_all no_covid_social ibn.id if group`g'==1 [iweight=w`g'], vce(robust) iterate(1000) nocons

		nlcom (wtp_ip: - _b[in_person]/_b[cost_all] * 1000) ///
			(wtp_s: - (_b[social] + _b[no_covid_social])/_b[cost_all] * 1000 ) ///
			(wtp_v: - _b[vaccine]/_b[cost_all] * 1000) ///
			(wtp_nc: - (_b[no_covid]+ _b[no_covid_social])/_b[cost_all] * 1000) ///
			, post
		
		foreach est in wtp_ip wtp_s wtp_nc wtp_v{
			local `est'`g': display %6.0fc _b[`est']
			quietly test `est'
			if r(p)<0.1{
				local `est'`g' = "``est'`g''*"
			}
			if r(p)<0.05{
				local `est'`g' = "``est'`g''*"
			}
			if r(p)<0.01{
				local `est'`g' = "``est'`g''*"
			}
			local `est'`g'_se: display %6.0fc _se[`est']
			local `est'`g'_se = subinstr("(``est'`g'_se')"," ","",.)
			local `est'`g'_pct: display %4.2f _b[`est']/`num_m`g'' * 100 
			local `est'`g'_pct = "[``est'`g'_pct'\%]"
		}
				
	}


	local title "WTP Estimates Robustness"
	********************
	* Table
	********************
	cap file close table
	local fwt "file write table"
	file open table using "${tables}table6.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:`name'} \setlength\tabcolsep{5pt} \centering \begin{threeparttable} \begin{tabularx}{\textwidth}{lYYYYYY} " _n
	
	`fwt' " & (1) & (2) & (3) & (4) & (5) & (6) \\"_n
	`fwt' " & `lab1' & `lab2' & `lab3' & `lab4' & `lab5'  & `lab6' \\ \hline"_n
	`fwt' " & \\"_n
	foreach var in wtp_s wtp_ip  {
		`fwt' "`lab`var'' "
		forval g = 1/6{
		`fwt' " & ``var'`g''"
		}
		`fwt' " \\"_n
		
		`fwt' " "
		forval g = 1/6{
		`fwt' " & ``var'`g'_se'"
		}
		`fwt' " \\"_n	
		`fwt' " "
		forval g = 1/6{
		`fwt' " & ``var'`g'_pct'"
		}
		`fwt' " \\"_n	
		
		`fwt' " \\"_n	
	}
		
	`fwt' " \hline" _n
	`fwt' " \\"_n	
	foreach stat in n m {
		`fwt' "`lab`stat'' "
		forval g = 1/6{
		`fwt' " & ``stat'`g''"
		}
		`fwt' " \\"_n
	}
	`fwt' " \\"_n
	`fwt' " \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} WTP reported in dollars. 
			Sample restricted/reweighted according to column header.
			Standard errors in parentheses derived via delta method.
			WTP as a percent of average cost in brackets.
			*, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively." ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table

}

************************
*	figures
************************

if `figure1'{
	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)
	keep if cost_1 == 1
	
	#delimit ;
	graph box prob_scenario
		, 
		over(hypo_num, relabel(	1 `""Controlled" "In-Person" "Social"' 
			2 `""Outbreak" "Online" "Restricted"' 
			3 `""Controlled" "Online" "Restricted"' 
			4 `""Controlled" "Online" "Social"' 
			5 `""Outbreak" "In-Person" "Social"' 
			6 `""Outbreak" "In-Person" "Social" "Vaccine"' )) 
		m(1, msize(*0.25) mc(gs10) mfc(none))
		box(1, color(gs2))
		graphregion(color(white))
		ytitle("Likelihood Scenario Realized")
		name(figure1, replace)
		;
	#delimit cr
	
	foreach g in "figure1"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
}

if `figure2'{
	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)

	gen group1 = 1 if female==1
	gen group2 = 1 if female==0
	
	gen group3 = 1 if low_inc==1
	gen group4 = 1 if low_inc==0

	gen group5 = 1 if honors == 0
	gen group6 = 1 if honors ==1
	
	local lab1_1 "Outbreak Controlled | In-Person"
	local lab1_2 "Normal Social"
	local lab2_1 "Outbreak Continues | Online"
	local lab2_2 "Restricted Social"
	local lab3_1 "Outbreak Controlled | Online"
	local lab3_2 "Restricted Social"
	local lab4_1 "Outbreak Controlled | Online"
	local lab4_2 "Normal Social"
	local lab5_1 "Outbreak Continues | In-Person"
	local lab5_2 "Normal Social"
	local lab6_1 "Outbreak Continues | In-Person"
	local lab6_2 "Normal Social | Vaccine"
	
	local paneln_1 "(a)"
	local paneln_2 "(b)"
	local paneln_3 "(c)"
	local paneln_4 "(d)"
	local paneln_5 "(e)"
	local paneln_6 "(f)"
	
	foreach h in 1 2 3 4 5 6{
		
		cap drop _*
		
		gen _g = _n if _n<=10
		gen _y = .
		gen _sd = .
		forval g = 1/6{
			sum cv if cost_1 == 1 & hypo_num == `h' & group`g' == 1
			replace _y = r(mean) if _g == `g'
			replace _sd = r(sd) if _g == `g'
		}
		gen _ub = _y + _sd
		gen _lb = _y - _sd
		
		gen _x = _n if _n<=6
		foreach i in 2 4{
			replace _x = _x + 0.5 if _g >`i'
		}
		
		*pattern data for bars
		cap drop _pattx 
		cap drop _patty
		cap drop _pattid
		gen _pattx = .
		gen _patty = .
		gen _pattid = .
		local j = 1
		forval g = 1/6{
			local jp1 = `j' + 1
			sum _x if _g == `g'
			replace _pattx = r(mean) - 0.45 if _n==`j'
			replace _pattx = r(mean) + 0.45 if _n==`jp1'
			
			sum _y if _g == `g'

			replace _patty = r(mean) if _n==`j'
			replace _patty = r(mean) if _n==`jp1'
			
			replace _pattid = `g' if _n==`j'
			replace _pattid = `g' if _n==`jp1'
			
			local j = `j' + 2
		}
		
		
		
		format _y %2.0f
		colorpalette  Paired, nograph

		#delimit ;
		twoway 
			(parea _patty _pattx if _pattid==1, color(black) pattern(none))
			(parea _patty _pattx if _pattid==2, color(black) pattern(pattern2))
			(parea _patty _pattx if _pattid==3, color(black) pattern(pattern5))
			(parea _patty _pattx if _pattid==4, color(black) pattern(pattern6))
			(parea _patty _pattx if _pattid==5, color(black) pattern(pattern8))
			(parea _patty _pattx if _pattid==6, color(black) pattern(pattern9))
			
			(scatter _y _x, msym(none) mlab(_y) mlabp(12) mlabc(black))
			
			//(rcap _ub _lb _x, lcolor(gs4))
			,
			ylabel(70(10)100, gmax angle(0))
			yline(70, lcolor(black))
			graphregion(color(white))
			xscale(noline)
			legend(order(1 "Female" 3 "Low Income" 5 "Non-Honors" 2 "Male" 4 "High Income"  6 "Honors") r(2) size(*0.6) symx(*0.4) )
			//xlabel(1 "Fe" 2 "Ma" 3.5 "LoI" 4.5 "HiI" 6 "NHo" 7 "Ho", notick labsize(*0.8) labgap(*0))
			xlabel(, notick nolabels)
			xtitle("`paneln_`h''") 
			title("`lab`h'_1'" "`lab`h'_2'", size(*0.8) color(black))
			name(plot`h', replace)
			nodraw
			ytitle(" ")
			legend(off)
			;
		#delimit cr
		
	}
	
	
	#delimit ;
	grc1leg plot1 plot2 plot3 plot4 plot5 plot6
		, 
		rows(3) 
		graphregion(color(white))  
		name(figure2, replace)
		l1title("Likelihood Return at Current Costs", size(*0.8))
		;
	#delimit cr
	
	graph display figure2, xsize(8) ysize(10)
	
	
	foreach g in "figure2"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
}

if `figure3'{
	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)
	
	keep if cost_1 == 1
	keep if hypo_num == 3 | hypo_num == 4
	
	gen a = hypo_num == 4
	
	keep ResponseId a cv
	sort ResponseId
		
	forval i = 0/1{
		preserve 
			keep if a == `i'
			ren cv cv`i'
			tempfile a`i'
			save `a`i'', replace
		restore
	}
	use `a1', clear
	merge 1:1 ResponseId using `a0'
	
	gen _x = 0 if _n == 1
	replace _x = 100 if _n == 2
	gen _y  = _x
	
	#delimit ;
	twoway (scatter cv1 cv0, mcolor(gs2%60) mlwidth(none) msize(*0.6) jitter(5) jitterseed(1))
		(line _y _x, lpattern(dash) lcolor(black))
		,
		graphregion(color(white))
		legend(off)
		xtitle("Controlled | Online | Restricted")
		ytitle("Controlled | Online | Social")
		name(figure3, replace)
		;
	#delimit cr
	
	foreach g in "figure3"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}

}

if `figure4'{

	*confidential ASU data
	import delimit "${restricted_data}enrollment.csv", clear rowrange(4:) varn(4)
	gen name = first_name + " " + last_name
	ren age_today age
	gen enroll_f20 = 1 
	gen enroll_f20_ip = campus_term_start != "ONLINE"
	keep emplid asurite enroll_f20_ip enroll_f20
	ren emplid idn
	
	tempfile enroll 
	save `enroll'
	
	use ${restricted_data}surveydata_restrict.dta, clear
	destring Q135, gen(idn)
	gen asurite = substr(email,1,strpos(email,"@")-1)
	
	merge m:1 asurite using `enroll'
	drop if _merge == 2
	drop _merge
	ren enroll_f20 _enroll_f20_1
	ren enroll_f20_ip _enroll_f20_ip_1
	
	merge m:1 idn using `enroll'
	drop if _merge == 2
	drop _merge
	ren enroll_f20 _enroll_f20_2
	ren enroll_f20_ip _enroll_f20_ip_2
	
	gen enroll_f20 = (_enroll_f20_1 == 1) | (_enroll_f20_2 == 1)
	gen enroll_f20_ip = (_enroll_f20_ip_1 == 1) | (_enroll_f20_ip_2 == 1)
	
	gen cv_100 = cv/100
	
	quietly reg enroll_f20_ip cv_100 if hypo_num == 2 & cost_1 == 1
	local b: display %5.4f _b[cv]
	local se: display %5.4f _se[cv]
	local se = subinstr("(`se')"," ","",.)
	quietly test cv
	if r(p) < 0.1{
		local b = "`b'*"
	}
	if r(p) < 0.05{
		local b = "`b'*"
	}
	if r(p) < 0.01{
		local b = "`b'*"
	}
	
	binscatter enroll_f20_ip cv_100 if hypo_num == 2 & cost_1 == 1, n(30) xtitle("Stated Likelihood of Enrolling in Fall 2020") ytitle("Share Actually Enrolled In-person in Fall 2020") name(figure4, replace)
	
	*save graphs
	foreach g in figure4{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
}

if `figure5'{
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace

	use ${data}indiv_WTP_with_shrink.dta, clear
	
	gen W_wtp_social = wtp_social_s
	gen W_wtp_inperson = wtp_in_person_s 
			
	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
	foreach f in "social" "inperson"{
		replace W_wtp_`f' = -4 if W_wtp_`f'<-4 & !missing(W_wtp_`f')
		replace W_wtp_`f' = 8 if W_wtp_`f'>8 & !missing(W_wtp_`f')
		
	}
	
	keep W_wtp_social W_wtp_inperson
	ren W_wtp_social wtp1
	ren W_wtp_inperson wtp2
	gen id = _n
	reshape long wtp, i(id) j(factor)
	#delimit ;
	distplot wtp
		, 
		over(factor)
		lpattern(solid dash)
		lcolor(gs9 gs2)
		lwidth(*2 *2)
		graphregion(color(white))
		xtitle("WTP ($1,000s/year)")
		ytitle("Density")
		title(" ")
		name(figure5, replace)
		xlabel(-4(2)8)
		legend(order(1 "Social" 2 "In-Person"))
		;
	#delimit cr
	
	*save graphs
	foreach g in "figure5"{
		graph export ${figures}figure5.png, name(`g') replace width(4000)
	}
}

if `figure6'{	
	use ${data}indiv_WTP_with_shrink.dta, clear
	
	local iplb -4000
	local ipub 8000
	local slb -4000
	local sub 8000
	
	gen W_wtp_social = wtp_social_s * 1000
	gen W_wtp_inperson = wtp_in_person_s  * 1000
	
	replace W_wtp_social = `sub' if W_wtp_social>`sub' & !missing(W_wtp_social)
	replace W_wtp_social = `slb' if W_wtp_social<`slb' & !missing(W_wtp_social)
	replace W_wtp_inperson = `ipub' if W_wtp_inperson>`ipub' & !missing(W_wtp_inperson)
	replace W_wtp_inperson = `iplb' if W_wtp_inperson<`iplb' & !missing(W_wtp_inperson)
	
	local lab1_low_inc "Lower-Income"
	local lab0_low_inc "Higher-Income"
	
	local title_social "WTP for Campus Social Life"
	local title_inperson "WTP for In-Person Instruction"
	
	local short_low_inc "li"
	local short_first_gen "fg"
	local short_nonwhite "nw"
	local short_female "fe"
	local short_honors "ho"
		
	local xlab_social "`slb'(2000)`sub'"
	local xlab_inperson "`iplb'(2000)`ipub'"

	local panel_social "(a)"
	local panel_inperson "(b)"

	foreach factor in social inperson{
		foreach demo in low_inc{
			cap drop wtp
			gen wtp = W_wtp_`factor'
			forval d=0/1{
				quietly sum wtp if `demo'==`d'
				local mean`d' = r(mean)
			}
			
			ksmirnov wtp, by(`demo')
			local p: display %4.3f r(p)
			#delimit ;
			distplot wtp
				, 
				over(`demo')
				lpattern(solid dash)
				lcolor(gs2 gs9)
				graphregion(color(white))
				xtitle("WTP", size(*0.7))
				ytitle("Cumulative Density", size(*0.7))
				title("`panel_`factor''" "`title_`factor''", size(*0.4))
				legend(order(1 "`lab0_`demo''" 2 "`lab1_`demo''") r(1) symx(*0.6) forces)
				name(wtp_`factor'_`short_`demo'', replace)
				xline(`mean0', lcolor(gs2) lpattern(solid))
				xline(`mean1', lcolor(gs9) lpattern(shortdash))
				xlabel(`xlab_`factor'', labsize(*0.7))
				ylabel(, labsize(*0.7))
				note("k-s pvalue: `p'", ring(0) pos(4) size(*0.7))
				;
			#delimit cr
		}
	}

	grc1leg wtp_social_li wtp_inperson_li, iscale(*2) graphregion(color(white)) name(figure6, replace) 
	graph display figure6, xsize(10) ysize(5)
	
	
	*save graphs
	foreach g in "figure6"{
		graph export ${figures}`g'.png, name(`g') replace width(5000) height(2500)
	}
}

if `figure7a'{
	use ${data}indiv_WTP_val.dta, clear
	
	merge 1:1 id using ${data}validation_real.dta
	
	gen coef1 = .
	gen coef2 = .
	
	gen _cv1 = .
	gen _cv2 = . 
	foreach i in 1 2{
		replace _cv`i' = cv1 if scen_n_drop`i' == 1
		replace _cv`i' = cv2 if scen_n_drop`i' == 2
		replace _cv`i' = cv3 if scen_n_drop`i' == 3
		replace _cv`i' = cv4 if scen_n_drop`i' == 4
		replace _cv`i' = cv5 if scen_n_drop`i' == 5
		replace _cv`i' = cv6 if scen_n_drop`i' == 6
				
		gen lod`i' = ln(_cv`i'/(100-_cv`i'))
		
		replace coef`i' = coef_no_covid_frac + coef_in_person_frac + coef_social_frac + coef_no_covid_social_frac if scen_n_drop`i' == 1
		replace coef`i' = 0 if scen_n_drop`i' == 2
		replace coef`i' = coef_no_covid_frac if scen_n_drop`i' == 3
		replace coef`i' = coef_no_covid_frac + coef_social_frac + coef_no_covid_social_frac if scen_n_drop`i' == 4
		replace coef`i' = coef_in_person_frac + coef_social_frac if scen_n_drop`i' == 5
		replace coef`i' = coef_in_person_frac + coef_social_frac + coef_vaccine_frac if scen_n_drop`i' == 6
	}
	
	gen Dlod = lod2 - lod1  
	
	gen Dcoef = (coef2 - coef1) if !missing(Dlod)
	
	
	gen _perdicted0 = 1 if Dcoef==0 //to ensure line isn't dropped later
	sort _perdicted0
	gen _y = .
	sum Dlod

	replace _y = -5 if _n == 2
	replace _y = 5 if _n == 1
	gen _x = _y
	
	count if !missing(Dlod)
	local N: display %5.0fc r(N)

	replace Dcoef = 5 if Dcoef > 5 & !missing(Dcoef)
	replace Dcoef = -5 if Dcoef < -5 & !missing(Dcoef)
	replace Dlod = 5 if Dlod > 5 & !missing(Dlod)
	replace Dlod = -5 if Dlod < -5 & !missing(Dlod)
	count if !missing(Dlod) 
	local N: display %5.0fc r(N)
	#delimit ;
	twoway 	(line _y _x, lcolor(black) lpattern(dash))
		(scatter Dlod Dcoef, mlwidth(none) mcolor(gs8%35) msize(*0.5))	
		,
		graphregion(color(white))
		plotregion(margin(zero))
		legend(off)
		xtitle("Predicted {&Delta} Log Odds")
		ytitle("Observed {&Delta} Log Odds")
		name(figure7a, replace)
		;
	#delimit cr

	
	foreach g in "figure7a"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
}

if `figure7b'{
	
	use ${data}surveydata.dta, clear
	set matsize 2000

	keep if !missing(cv) & !missing(cost_all)
	
	merge m:1 id using ${data}validation_real.dta, nogen
	
	drop if (hypo_num == scen_n_drop1 | hypo_num == scen_n_drop2) & cost == cost_drop
	
	gen cv_prob = cv/100

	fracreg logit cv_prob no_covid in_person social vaccine cost_all no_covid_social ibn.id, vce(robust) iterate(1000) nocons
	
	foreach coef in no_covid in_person social vaccine cost_all no_covid_social{
		gen coef_`coef'_pool = _b[`coef'] if _n == 1
	}
	keep if _n == 1
	keep coef_*_pool
	gen one = 1
	save ${data}pool_coef.dta, replace
	
	********
	
	use ${data}indiv_WTP_val.dta, clear
	
	merge 1:1 id using ${data}validation_real.dta, nogen
	
	gen one = 1
	merge m:1 one using ${data}pool_coef.dta, nogen
	
	drop if novarflag==1
	
	foreach i in 1 2{
		gen _cv`i' = .
		replace _cv`i' = cv1 if scen_n_drop`i' == 1
		replace _cv`i' = cv2 if scen_n_drop`i' == 2
		replace _cv`i' = cv3 if scen_n_drop`i' == 3
		replace _cv`i' = cv4 if scen_n_drop`i' == 4
		replace _cv`i' = cv5 if scen_n_drop`i' == 5
		replace _cv`i' = cv6 if scen_n_drop`i' == 6
				
		gen lod`i' = ln(_cv`i'/(100-_cv`i'))
		
		foreach m in frac pool{
			gen coef`i'_`m' = .
			replace coef`i'_`m' = coef_no_covid_`m'+ coef_in_person_`m'+ coef_social_`m'+ coef_no_covid_social_`m' if scen_n_drop`i' == 1
			replace coef`i'_`m' = 0 if scen_n_drop`i' == 2
			replace coef`i'_`m' = coef_no_covid_`m' if scen_n_drop`i' == 3
			replace coef`i'_`m' = coef_no_covid_`m'+ coef_social_`m'+ coef_no_covid_social_`m' if scen_n_drop`i' == 4
			replace coef`i'_`m' = coef_in_person_`m'+ coef_social_`m' if scen_n_drop`i' == 5
			replace coef`i'_`m' = coef_in_person_`m'+ coef_social_`m'+ coef_vaccine_`m' if scen_n_drop`i' == 6
		}
	}
	
	gen Dlod = lod2 - lod1  
	
	gen Dcoef_frac = (coef2_frac - coef1_frac) if !missing(Dlod)
	gen Dcoef_pool = (coef2_pool - coef1_pool) if !missing(Dlod)
	
	gen D_frac = Dcoef_frac - Dlod
	gen D_pool = Dcoef_pool - Dlod
	
	drop if missing(D_frac) | missing(D_pool)
		

	foreach var in D_frac D_pool{
		replace `var' = 4 if `var'>4 & !missing(`var')
		replace `var' = -4 if `var'<-4 & !missing(`var')
	}
	
	#delimit ;
	twoway 	
		(hist D_pool, width(0.2) fcolor(gs8) lwidth(none) lcolor(gs8)  lwidth(*1.4))	
		(hist D_frac, width(0.2) fcolor(none) lcolor(black))	
		,
		graphregion(color(white))
		plotregion(margin(zero))
		xlabel(-4(1)4)
		xtitle("Predicted {&Delta} Log Odds - Observed {&Delta} Log Odds")
		title(" ")
		legend(order(2 "Individual" 1 "Pooled"))
		name(figure7b, replace)
		;
	#delimit cr
	
	foreach g in "figure7b"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}

}

************************
*	appendix tables
************************

if `tableA1'{
	use ${data}surveydata.dta, clear
	set matsize 2000
	
	keep if !missing(cv) & !missing(cost_all)
	
	gen cv_prob = cv/100
	
	gen covid = no_covid == 0 if !missing(no_covid)
	gen covid_social = covid * social
	   
	cap drop _indiv_tag
	bys ResponseId: gen _indiv_tag=1 if _n==1
	
	
	foreach d in low_inc{
		gen `d'in_person = `d' * in_person
		gen `d'social = `d' * social
		gen `d'cost_all = `d' * cost_all
		gen `d'covid_social = `d' * covid_social
		gen `d'covid = `d' * covid
		gen `d'vaccine = `d' * vaccine
	}
	
	local m = 1
	fracreg logit cv_prob covid in_person social vaccine cost_all covid_social ibn.id, vce(robust) iterate(1000) nocons
	foreach e in covid in_person social vaccine cost_all covid_social{
		local `e'`m': display %6.2fc _b[`e']
		quietly test `e'
		if r(p)<0.1{
			local `e'`m' = "``e'`m''*"
		}
		if r(p)<0.05{
			local `e'`m' = "``e'`m''*"
		}
		if r(p)<0.01{
			local `e'`m' = "``e'`m''*"
		}
		local `e'`m'_se: display %6.2fc _se[`e']
		local `e'`m'_se = subinstr("(``e'`m'_se')"," ","",.)
	}

	local m = 2
	foreach d in low_inc{
		fracreg logit cv_prob covid in_person social vaccine cost_all covid_social `d'in_person `d'social `d'cost_all `d'covid_social `d'covid `d'vaccine ibn.id, vce(robust) iterate(1000) nocons
		foreach e in covid in_person social vaccine cost_all covid_social{
			local `e'`m': display %6.2fc _b[`e']
			quietly test `e'
			if r(p)<0.1{
				local `e'`m' = "``e'`m''*"
			}
			if r(p)<0.05{
				local `e'`m' = "``e'`m''*"
			}
			if r(p)<0.01{
				local `e'`m' = "``e'`m''*"
			}
			local `e'`m'_se: display %6.2fc _se[`e']
			local `e'`m'_se = subinstr("(``e'`m'_se')"," ","",.)
		}
		
		foreach e in covid in_person social vaccine cost_all covid_social{
			local `e'`m'd: display %6.2fc _b[`d'`e']
			quietly test `d'`e'
			if r(p)<0.1{
				local `e'`m'd = "``e'`m'd'*"
			}
			if r(p)<0.05{
				local `e'`m'd = "``e'`m'd'*"
			}
			if r(p)<0.01{
				local `e'`m'd = "``e'`m'd'*"
			}
			local `e'`m'd_se: display %6.2fc _se[`d'`e']
			local `e'`m'd_se = subinstr("(``e'`m'd_se')"," ","",.)
		}
		
		local m = `m' + 1
	
	}

	local title "Pooled Model Coefficients: Income Interaction"
	
	local lab_covid "COVID"
	local lab_in_person "In-person"
	local lab_social "Social"
	local lab_vaccine "Vaccine"
	local lab_cost_all "Cost"
	local lab_covid_social "COVID x Social"
	
	local lab_low_inc "Lower-inc."
	local lab_first_gen "First-gen."
	local lab_nonwhite "Nonwhite"
	local lab_female "Female"
	local lab_honors "Honors"
	local lab_fresh "Freshman"

	foreach d in low_inc first_gen nonwhite female honors fresh{
		foreach e in covid in_person social vaccine cost_all covid_social{
			local lab_`d'`e' "`lab_`e'' x `lab_`d''"
		}
	}
	
	********************
	* Main Table
	********************
	cap file close table
	local fwt "file write table"
	file open table using "${tables}tableA1.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:`name'} \setlength\tabcolsep{3pt} \centering \footnotesize \begin{threeparttable} \begin{tabularx}{0.6\textwidth}{lYY} " _n
	`fwt' " & (1) & (2)  \\"_n
	`fwt' " & \\ \hline"_n
	`fwt' " & \\"_n

	foreach e in social in_person cost_all covid covid_social vaccine{
		`fwt' " `lab_`e''"
		forval m = 1/2{
			`fwt' " & ``e'`m''"
		}
		`fwt' "\\"_n
		forval m = 1/2{
			`fwt' " & ``e'`m'_se'"
		}
		`fwt' "\\"_n
		`fwt' " `lab_`e'' x Lower-income"
		forval m = 1/2{
			`fwt' " & ``e'`m'd'"
		}
		`fwt' "\\"_n
		forval m = 1/2{
			`fwt' " & ``e'`m'd_se'"
		}
		`fwt' "\\"_n
		`fwt' "\\"_n
	}

	`fwt' " \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} Table displays fractional response model coefficient estimates. Column 1 corresponds to the baseline pooled model.
			Column 2 present coefficients from an alternative specification where each dependent variable is interacted with a dummy for lower-income status.
			Both models include individual fixed effects.
			*, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. " ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table

}

if `tableA2'{
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace
	
	use ${data}indiv_WTP.dta, clear
	
	gen rawW_wtp_social = wtp_social_frac * 1000
	gen rawW_wtp_inperson = wtp_in_person_frac * 1000			
	
	*winsorize correctly consistent with shrinkage estimates 
	foreach f in inperson social{
		_pctile rawW_wtp_`f', p(3 97)
		local r1 = r(r1)
		local r2 = r(r2)
		sum rawW_wtp_`f' if rawW_wtp_`f' >=`r1' & rawW_wtp_`f'<=`r2'
		replace rawW_wtp_`f' = r(mean) if rawW_wtp_`f'<`r1' & !missing(rawW_wtp_`f')
		replace rawW_wtp_`f' = r(mean) if  rawW_wtp_`f'>`r2' & !missing(rawW_wtp_`f')
	}
	
	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
	gen group1 = 1
	
	gen group2 = 1 if low_inc==1
	gen group3 = 1 if low_inc==0

	gen group4 = 1 if first_gen==1
	gen group5 = 1 if first_gen==0

	gen group6 = 1 if nonwhite==1
	gen group7 = 1 if nonwhite==0
		
	gen group8 = 1 if female==1
	gen group9 = 1 if female==0
	
	gen group10 = 1 if honors==0
	gen group11 = 1 if honors==1
	
	gen group12 = 1 if cohort==1
	gen group13 = 1 if cohort==2
	gen group14 = 1 if cohort==3
	
	local lab1 "All"

	local lab2 "Lower-inc."
	local lab3 "Higher-inc."
	local lab4 "First-gen."
	local lab5 "Second-gen."
	local lab6 "Nonwhite"
	local lab7 "White"
	local lab8 "Females"
	local lab9 "Males"
	local lab11 "Honors"
	local lab10 "Non-honors"
	local lab12 "2023+ Grads"
	local lab13 "2022 Grads"
	local lab14 "2021 Grads"
	
	local abrev_W_wtp_social "s"
	local abrev_W_wtp_inperson "ip"
	local abrev_rawW_wtp_social "raws"
	local abrev_rawW_wtp_inperson "rawip"
	
	sum cost2019 if group1 == 1
	local m_cost = r(mean)*1000
	local m_cost_format: display %6.0fc r(mean)*1000
	
	forval g=1/14{
		foreach wtp in rawW_wtp_social rawW_wtp_inperson{
			/*
			sum cost2019 if group`g' == 1
			local m_cost = r(mean)*1000
			*/
			
			sum `wtp' if group`g' == 1, d
			local wtp_`abrev_`wtp''_g`g'_p50: display %6.0fc r(p50)
			local wtp_`abrev_`wtp''_g`g'_p10: display %6.0fc r(p10)
			local wtp_`abrev_`wtp''_g`g'_p90: display %6.0fc r(p90)
			local wtp_`abrev_`wtp''_g`g'_sd: display %6.0fc r(sd)
			/*
			local wtp_`abrev_`wtp''_g`g'_m: display %6.0fc r(mean)
			local wtp_`abrev_`wtp''_g`g'_pct: display %4.2f r(mean)/`m_cost' * 100
			local wtp_`abrev_`wtp''_g`g'_pct = subinstr("[`wtp_`abrev_`wtp''_g`g'_pct'\%]"," ","",.)
			*/
			
			reg `wtp' if group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_m: display %6.0fc _b[_cons]
			local wtp_`abrev_`wtp''_g`g'_se: display %6.0fc _se[_cons]
			local wtp_`abrev_`wtp''_g`g'_se = subinstr("(`wtp_`abrev_`wtp''_g`g'_se')"," ","",.)
			local wtp_`abrev_`wtp''_g`g'_m = "`wtp_`abrev_`wtp''_g`g'_m'\,`wtp_`abrev_`wtp''_g`g'_se'"
					
			local wtp_`abrev_`wtp''_g`g'_pct: display %4.2f _b[_cons]/`m_cost' * 100
			local wtp_`abrev_`wtp''_g`g'_pct = subinstr("[`wtp_`abrev_`wtp''_g`g'_pct'\%]"," ","",.)

			test _cons
			local stars = ""
			if r(p) < 0.1{
				local stars = "`stars'*"
			}
			if r(p) < 0.05{
				local stars = "`stars'*"
			}
			if r(p) < 0.01{
				local stars = "`stars'*"
			}
			local wtp_`abrev_`wtp''_g`g'_m = "$\stackrel{`stars'}{`wtp_`abrev_`wtp''_g`g'_m'}$"			
			

			
			
			count if !missing(`wtp') & group`g' == 1
			local n = r(N)
			count if `wtp'>10 & !missing(`wtp') & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_g1: display %3.1f r(N)/`n' *100
			
			count if `wtp'<-10 & !missing(`wtp')  & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_l1: display %3.1f r(N)/`n' *100
			
			
		}
	}
	
	forval g = 2(2)10{
		local g2 = `g' + 1
		gen _group = 0 if group`g' == 1
		replace _group = 1 if group`g2' == 1
		
		foreach wtp in rawW_wtp_social rawW_wtp_inperson{
			ttest `wtp', by(_group)
			local p_`abrev_`wtp''_`g2': display %4.3f r(p) 
			local p_`abrev_`wtp''_`g2' = "p-val: `p_`abrev_`wtp''_`g2''"
		}
		cap drop _group
	}
	
	*anova for cohort (3 means, not 2)
	foreach wtp in rawW_wtp_social rawW_wtp_inperson{
		oneway `wtp' cohort
		local p_`abrev_`wtp''_14: display %4.3f Ftail(`r(df_m)',`r(df_r)',`r(F)')
		local p_`abrev_`wtp''_14 = "p-val: `p_`abrev_`wtp''_14'"
	} 
	
	
	local title "Moments of Individual WTP Distribution: Unadjusted"
	cap file close table
	local fwt "file write table"
	file open table using "${tables}tableA2_v3.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:tableA2} \setlength\tabcolsep{1pt} \centering \begin{threeparttable} \begin{tabularx}{1.4\textwidth}{lcYYYYYYYcYYYYYY} " _n
	`fwt' " & \multicolumn{7}{c}{WTP Social} && \multicolumn{7}{c}{WTP In-Person}\\ \cline{2-8} \cline{10-16}"_n
	`fwt' " & \\"_n
	`fwt' " & Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10 && Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10  \\ \hline"_n
	`fwt' " & \\"_n


	forval g=1/14{
		`fwt' " \qquad `lab`g'' "
		foreach wtp in raws rawip{
			foreach mom in m p50 sd p10 p90 g1 l1{
				`fwt' "& `wtp_`wtp'_g`g'_`mom''"
			}
			if "`wtp'" == "raws"{
				`fwt' "& "
			}
		}
		`fwt' "\\"_n
		//`fwt' " & `wtp_s_g`g'_se' &&&&&&&& `wtp_ip_g`g'_se' &&&&&& \\"_n
		
		if "`g'"=="1"{
			`fwt' " & `wtp_raws_g1_pct' &&&&&&&& `wtp_rawip_g1_pct' &&&&&& \\"_n

		}
		if inlist(`g',3,5,7,9,11,14){
			`fwt' " & `p_raws_`g'' &&&&&&&& `p_rawip_`g'' &&&&&& \\"_n
		}
		
		if inlist(`g',1,3,5,7,9,11){
			`fwt' " & \\"_n
		}
	} 
	
	`fwt' " &  \\ \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt'"\item[] \footnotesize
			\emph{Notes:} Willingness-to-Pay measured in dollars per year. 
			Standard error of mean estimate in parentheses. 
			P-value for equality of means between demographic groups displayed below each set of means.
			Willingness-to-Pay as a percent of average net costs (\textdollar`m_cost_format') displayed in brackets.
			Grads refers to student expected cohort by graduation year.
			Within means column *, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. " ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table

	
}

if `tableA3'{	
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace

	use ${data}indiv_WTP_with_shrink.dta, clear
	
	count if novarflag==1
	local nnovar: display %4.0fc r(N)
	
	foreach f in social in_person{
				
		gen winsorSEflag_`f' = 0 
		_pctile wtpse_`f'_frac, p(95)
		local r3 = r(r1)
		replace  winsorSEflag_`f' = 1 if wtpse_`f'_frac>`r3' & !missing(wtp_`f'_frac)
		
		
		gen winsorflag_`f' = 0 
		_pctile wtp_`f'_frac if novarflag!=1 & winsorSEflag_`f'!=1, p(3 97)
		local r1 = r(r1)
		local r2 = r(r2)
		replace winsorflag_`f'=1 if (wtp_`f'_frac>`r2' | wtp_`f'_frac<`r1') & !missing(wtp_`f'_frac) & novarflag!=1 & winsorSEflag_`f' != 1
		
	}
	
	count if (winsorflag_social==1 | winsorflag_in_person==1)
	local nwinsor: display %4.0fc r(N)
	
	count if (winsorSEflag_social==1 | winsorSEflag_in_person==1)
	local nSEwinsor: display %4.0fc r(N)
	
	gen W_wtp_social = wtp_social_s * 1000
	gen W_wtp_inperson = wtp_in_person_s * 1000
	
		
	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId *flag*
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
	drop if novarflag==1 | winsorflag_social==1 | winsorflag_in_person==1 | winsorSEflag_social==1 | winsorSEflag_in_person==1
	
	gen group1 = 1
	
	gen group2 = 1 if low_inc==1
	gen group3 = 1 if low_inc==0

	gen group4 = 1 if first_gen==1
	gen group5 = 1 if first_gen==0

	gen group6 = 1 if nonwhite==1
	gen group7 = 1 if nonwhite==0
		
	gen group8 = 1 if female==1
	gen group9 = 1 if female==0
	
	gen group10 = 1 if honors==0
	gen group11 = 1 if honors==1
	
	gen group12 = 1 if cohort==1
	gen group13 = 1 if cohort==2
	gen group14 = 1 if cohort==3
	
	local lab1 "All"

	local lab2 "Lower-inc."
	local lab3 "Higher-inc."
	local lab4 "First-gen."
	local lab5 "Second-gen."
	local lab6 "Nonwhite"
	local lab7 "White"
	local lab8 "Females"
	local lab9 "Males"
	local lab11 "Honors"
	local lab10 "Non-honors"
	local lab12 "2023+ Grads"
	local lab13 "2022 Grads"
	local lab14 "2021 Grads"
	
	
	local abrev_W_wtp_social "s"
	local abrev_W_wtp_inperson "ip"
	local abrev_rawW_wtp_social "raws"
	local abrev_rawW_wtp_inperson "rawip"
	
	sum cost2019 if group1 == 1
	local m_cost = r(mean)*1000
	local m_cost_format: display %6.0fc r(mean)*1000
	
	forval g=1/14{
		foreach wtp in W_wtp_social W_wtp_inperson{

			
			sum `wtp' if group`g' == 1, d
			local wtp_`abrev_`wtp''_g`g'_p50: display %6.0fc r(p50)
			local wtp_`abrev_`wtp''_g`g'_p10: display %6.0fc r(p10)
			local wtp_`abrev_`wtp''_g`g'_p90: display %6.0fc r(p90)
			local wtp_`abrev_`wtp''_g`g'_sd: display %6.0fc r(sd)

			
			reg `wtp' if group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_m: display %6.0fc _b[_cons]
			local wtp_`abrev_`wtp''_g`g'_se: display %6.0fc _se[_cons]
			local wtp_`abrev_`wtp''_g`g'_se = subinstr("(`wtp_`abrev_`wtp''_g`g'_se')"," ","",.)
			local wtp_`abrev_`wtp''_g`g'_m = "`wtp_`abrev_`wtp''_g`g'_m'\,`wtp_`abrev_`wtp''_g`g'_se'"
					
			local wtp_`abrev_`wtp''_g`g'_pct: display %4.2f _b[_cons]/`m_cost' * 100
			local wtp_`abrev_`wtp''_g`g'_pct = subinstr("[`wtp_`abrev_`wtp''_g`g'_pct'\%]"," ","",.)

			test _cons
			local stars = ""
			if r(p) < 0.1{
				local stars = "`stars'*"
			}
			if r(p) < 0.05{
				local stars = "`stars'*"
			}
			if r(p) < 0.01{
				local stars = "`stars'*"
			}
			local wtp_`abrev_`wtp''_g`g'_m = "$\stackrel{`stars'}{`wtp_`abrev_`wtp''_g`g'_m'}$"			
			

			
			
			count if !missing(`wtp') & group`g' == 1
			local n = r(N)
			count if `wtp'>10 & !missing(`wtp') & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_g1: display %3.1f r(N)/`n' *100
			
			count if `wtp'<-10 & !missing(`wtp')  & group`g' == 1
			local wtp_`abrev_`wtp''_g`g'_l1: display %3.1f r(N)/`n' *100
			
			
		}
	}
	
	forval g = 2(2)10{
		local g2 = `g' + 1
		gen _group = 0 if group`g' == 1
		replace _group = 1 if group`g2' == 1
		
		foreach wtp in W_wtp_social W_wtp_inperson rawW_wtp_social rawW_wtp_inperson{
			ttest `wtp', by(_group)
			local p_`abrev_`wtp''_`g2': display %4.3f r(p) 
			local p_`abrev_`wtp''_`g2' = "p-val: `p_`abrev_`wtp''_`g2''"
		}
		cap drop _group
	}
	
	*anova for cohort (3 means, not 2)
	foreach wtp in W_wtp_social W_wtp_inperson rawW_wtp_social rawW_wtp_inperson{
		oneway `wtp' cohort
		local p_`abrev_`wtp''_14: display %4.3f Ftail(`r(df_m)',`r(df_r)',`r(F)')
		local p_`abrev_`wtp''_14 = "p-val: `p_`abrev_`wtp''_14'"
	} 
	
		
	*********
	*shrinkage and dropping outliers and zeros (short table)
	*********
	
	local title "Moments of Individual WTP Distribution: Excluding Outlier Responses"
	cap file close table
	local fwt "file write table"
	file open table using "${tables}tableA3.tex", write replace
	`fwt' "\begin{table}[htp]\caption{`title'}\label{t:indiv_tab2_R} \setlength\tabcolsep{1pt} \centering \begin{threeparttable} \begin{tabularx}{1.4\textwidth}{lcYYYYYYYcYYYYYY} " _n
	`fwt' " & \multicolumn{7}{c}{WTP Social} && \multicolumn{7}{c}{WTP In-Person}\\ \cline{2-8} \cline{10-16}"_n
	`fwt' " & \\"_n
	`fwt' " & Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10 && Mean & Med & StDev & $10^{th}$ & $90^{th}$ & \%$>$\textdollar 10 & \%$<$-\textdollar 10  \\ \hline"_n
	`fwt' " & \\"_n


	forval g=1/14{
		`fwt' " \qquad `lab`g'' "
		foreach wtp in s ip{
			foreach mom in m p50 sd p10 p90 g1 l1{
				`fwt' "& `wtp_`wtp'_g`g'_`mom''"
			}
			if "`wtp'" == "s"{
				`fwt' "& "
			}
		}
		`fwt' "\\"_n
		//`fwt' " & `wtp_s_g`g'_se' &&&&&&&& `wtp_ip_g`g'_se' &&&&&& \\"_n
		
		if "`g'"=="1"{
			`fwt' " & `wtp_s_g1_pct' &&&&&&&& `wtp_ip_g1_pct' &&&&&& \\"_n

		}
		if inlist(`g',3,5,7,9,11,14){
			`fwt' " & `p_s_`g'' &&&&&&&& `p_ip_`g'' &&&&&& \\"_n
		}
		
		if inlist(`g',1,3,5,7,9,11){
			`fwt' " & \\"_n
		}
	} 
	
	`fwt' " &  \\ \hline" _n

	`fwt' "\end{tabularx}"_n
	`fwt'"\begin{tablenotes}"_n
	#delimit ; 
	`fwt' "\item[] \footnotesize
			\emph{Notes:} Table presents moments from the distribution of individually estimated WTP measures, dropping observations without enough 
			variation to estimate preferences (n=`nnovar') as well as those at or above the 95th percentile of WTP SEs (n=`nSEwinsor'), and below the 3rd and above the 97th percentile of WTP estimantes (n=`nwinsor').
			Bayesian shrinkage algorithm applied to estimates. 
			Willingness-to-Pay measured in dollars per year. 
			Standard error of mean estimate in parentheses. 
			P-value for equality of means between demographic groups displayed below each set of means.
			Willingness-to-Pay as a percent of average net costs (\textdollar`m_cost_format') displayed in brackets.
			Grads refers to student expected cohort by graduation year.
			Within means column  *, **, *** denote estimates are statistically significant at the 10\%, 5\%, and 1\% levels, respectively. " ;
	#delimit cr
	`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
	file close table

	

}

if `tableA4'{
	use ${data}surveydata.dta, clear
	bys ResponseId: gen _tag = 1 if _n == 1
	keep if _tag == 1
	tempfile indiv_data
	cap drop _merge
	save `indiv_data', replace
		
	use ${data}indiv_WTP_with_shrink.dta, clear

	gen W_wtp_social = wtp_social_s * 1000
	gen W_wtp_inperson = wtp_in_person_s * 1000

	*merge back in all observables
	keep *wtp_social *wtp_inperson ResponseId
	merge 1:1 ResponseId using `indiv_data'
	keep if _merge == 3
	
	gen nhonors = 1 - honors
	
	
	
 	cap program drop covtab
	program covtab
		syntax [, all demos(string) name(string) title(string) econ(string) noecon(string)]
		
		local lab_female1 "Female"
		local lab_low_inc1 "Lower- inc."
		local lab_nonwhite1 "Non- white"
		local lab_first_gen1 "First- gen."
		local lab_nhonors1 "Non- honors"
		local lab_fresh1 "Fresh."
		
		local lab_female0 "Male"
		local lab_low_inc0 "Higher- inc."
		local lab_nonwhite0 "White"
		local lab_first_gen0 "Second- gen."
		local lab_nhonors0 "Honors"
		local lab_fresh0 "$\geq$ Soph."
		
		foreach var in `econ' `noecon'{
			local lab_`var': var lab `var' 
		}
		
		tokenize "`demos'"
		
		local ndemo = wordcount("`demos'")
	
		forval d = 1/`ndemo'{
			foreach var in `econ' `noecon'{
				forval i = 0/1{
					clear results
					cap quietly reg `var' if ``d'' == `i'
					if abs(_b[_cons])>1{
						local fmt %4.1f
					}
					else{
						local fmt %4.2f
					}
					cap local `var'_d`d'_i`i': display `fmt' _b[_cons]
				}
				quietly ttest `var', by(``d'')
				local `var'_d`d'_p: display %3.2f r(p)
				local `var'_d`d'_p = "(``var'_d`d'_p')"
			}
		}
		
		if "`all'" != ""{
			foreach var in `econ' `noecon'{
				clear results
				cap quietly reg `var' 
				if abs(_b[_cons])>1{
					local fmt %4.1f
				}
				else{
					local fmt %4.2f
				}
				cap local `var'_d0: display `fmt' _b[_cons]
			}
		}
		
		if "`all'" != ""{
			local cols "@{\hskip 0.5em}Y"
		}
		else{
			local cols 
		}
		
		forval d = 1/`ndemo'{
			local cols = "`cols'@{\hskip 0.5em}YYc"
		}
		
		cap file close table
		local fwt "file write table"
		file open table using "${tables}`name'.tex", write replace
		`fwt' "\begin{table}[htp]\caption{`title'}\label{t:`name'} \footnotesize \setlength\tabcolsep{1pt} \centering \begin{threeparttable} \begin{tabularx}{1.05\linewidth}{l`cols'} \\" _n
		
		if "`all'" != ""{
			`fwt' " & All"
		}
		forval d = 1/`ndemo'{
			`fwt' "& `lab_``d''1' & `lab_``d''0' & "
		}
		`fwt' "\\ "_n
		`fwt' " & "
		forval d = 1/`ndemo'{
			`fwt' "&   &   & p-val"
		}
		`fwt' "\\ \hline"_n
		`fwt' "\\ "_n
		
		`fwt' " \textbf{Economic}& \\"_n
		foreach var in `econ'{
			if "`all'" != ""{
				`fwt' " `lab_`var'' & ``var'_d0'"
			}
			else{
				`fwt' " `lab_`var''"
			}
			forval d = 1/`ndemo'{
				`fwt' "& ``var'_d`d'_i1'  & ``var'_d`d'_i0' & ``var'_d`d'_p'"
			}
			`fwt' "\\[0.5em]"_n
		}
		`fwt' "\\ "_n
		`fwt' " \textbf{Non-Economic}& \\"_n
		foreach var in `noecon'{
			if "`all'" != ""{
				`fwt' " `lab_`var'' & ``var'_d0'"
			}
			else{
				`fwt' " `lab_`var'' & ``var'_d0'"
			}
			forval d = 1/`ndemo'{
				`fwt' "& ``var'_d`d'_i1'  & ``var'_d`d'_i0' & ``var'_d`d'_p'"
			}
			`fwt' "\\[0.5em]"_n
		}
	
		
		`fwt' " &  \\ \hline" _n

		`fwt' "\end{tabularx}"_n
		`fwt'"\begin{tablenotes}"_n
		#delimit ; 
		`fwt'"\item[] \small
				\emph{Notes:} Table displays mean covariate levels for various demographic divisions. P-values for equality of means in parentheses." ;
		#delimit cr
		`fwt' "\end{tablenotes} \end{threeparttable} \end{table}" _n
		file close table
	end
		
	local demos "low_inc first_gen nonwhite female nhonors fresh"
	local econ_contr "work work_hrs20 total_pay_an grants loans cost_above_instate "
	local noecon_contr "live_campus_c0 social_c0 study_hrs_c0 act_eq class_oh_an2"
	
	covtab , econ(`econ_contr') noecon(`noecon_contr') all demos("`demos'") name("tableA4") title("Covariate Summary Statistics")
	
}

************************
*	appendix figures
************************

if `figureA1a'{
	use ${data}surveydata.dta, clear
	keep if !missing(cv) & !missing(cost_all)
	drop if hypo_num == 2
	
	cap drop _*
	gen _x = _n if _n<=7
	replace _x = 5-_x if _x<=4
	gen _q = _n if _n<=7
	gen _y = .
	gen _sd = .
	gen _med = .
	gen _p25 = .
	gen _p75 = .
	gen _SE = . 
	forval q = 1/7{
		sum cv if cost_`q'==1 & hypo_num==1, d
		replace _y = r(mean) if _q == `q'
		replace _sd = r(sd) if _q == `q'
		replace _med = r(p50) if _q == `q'
		replace _p25 = r(p25) if _q == `q'
		replace _p75 = r(p75) if _q == `q'
		
		reg cv if cost_`q'==1 & hypo_num==1
		replace _SE = _se[_cons] if _q == `q'
	}
	gen _ub = _y + 1.96 * _SE
	gen _lb = _y - 1.96 * _SE
	
	gen _ubmed = _med
	gen _lbmed = _med - _sd
	
	format _y _med %2.0f
	
	#delimit ;
	twoway (bar _y _x, fc(gs10) barw(0.9) lw(none))
		(scatter _ub _x, msym(none) mlab(_y) mlabp(12) mlabc(gs4))
		//(rcap _p75 _p25 _x, lcolor(gs4))
		(rcap _ub _lb _x, lcolor(gs4))
		,
		ylabel(0(20)100, gmax)
		yline(0, lcolor(black))
		graphregion(color(white))
		xscale(noline)
		legend(off)
		xlabel(1 "+30%" 2 "+20%" 3 "+10%" 4 "+0%" 5 "-10%" 6 "-20%" 7 "-30%", notick)
		xtitle("Cost of Attendance") 
		ytitle("Likelihood Return") 
		title("")
		name("figureA1a", replace)
		;
	#delimit cr
	
	foreach g in "figureA1a"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
	

}

if `figureA1b'{

	use ${data}surveydata.dta, clear
	replace cost2019 = cost2019*1000
	format %7.0fc cost2019
	#delimit ;
	twoway 	(hist cost2019 if _indiv_tag==1, color(navy))
		,
		graphregion(color(white))
		name(figureA1b, replace)
		ylabel(,nolab notick)
		xlabel(-20000(20000)40000)
		xtitle("Net Cost of Attendance ($/year)")
		;
	#delimit cr
	
	foreach g in "figureA1b"{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
	

}

if `figureA2'{
	
	*confidential ASU data
	import delimit "${restricted_data}enrollment.csv", clear rowrange(4:) varn(4)
	gen name = first_name + " " + last_name
	ren age_today age
	gen enroll_f20 = 1 
	gen enroll_f20_ip = campus_term_start != "ONLINE"
	keep emplid asurite enroll_f20_ip enroll_f20
	ren emplid idn
	
	tempfile enroll 
	save `enroll'
	
	use ${restricted_data}surveydata_restrict.dta, clear
	destring Q135, gen(idn)
	gen asurite = substr(email,1,strpos(email,"@")-1)
	
	merge m:1 asurite using `enroll'
	drop if _merge == 2
	drop _merge
	ren enroll_f20 _enroll_f20_1
	ren enroll_f20_ip _enroll_f20_ip_1
	
	merge m:1 idn using `enroll'
	drop if _merge == 2
	drop _merge
	ren enroll_f20 _enroll_f20_2
	ren enroll_f20_ip _enroll_f20_ip_2
	
	gen enroll_f20 = (_enroll_f20_1 == 1) | (_enroll_f20_2 == 1)
	gen enroll_f20_ip = (_enroll_f20_ip_1 == 1) | (_enroll_f20_ip_2 == 1)
	
	gen cv_100 = cv/100

	keep if hypo_num == 2 & cost_1 == 1
	
	merge 1:1 ResponseId using ${data}indiv_WTP_with_shrink.dta
	cap drop _merge
	
	
	preserve 
		keep if hypo_num == 2 & cost_1 == 1
		
		gen _tpr = . 
		gen _fpr = . 
		gen _T = (_n-1)/300 if _n<=305
		replace _T = round(_T, 0.0001)
		
		levelsof _T, local(ts)
		
		forval i = 1/305{
			local t = _T[`i']
		
			gen _cutDec = cv_100>=`t' if !missing(cv_100)
			
			count if enroll_f20_ip==1
			local P=r(N)
			count if _cutDec == 1 & enroll_f20_ip==1
			local TP=r(N)
			replace _tpr = `TP'/`P' if _n==`i'
			
			count if enroll_f20_ip==0
			local N=r(N)
			count if _cutDec == 1 & enroll_f20_ip==0
			local FN=r(N)
			replace _fpr = `FN'/`N' if _n==`i'
			
			cap drop _cutDec
		}
		
		roctab enroll_f20_ip cv_100
		local auc: display %4.3f r(area)
		
		gen _xline = _n-1 if _n<=2
		gen _yline = _xline
		
		#delimit ;
		twoway (line _tpr _fpr, lwidth(*2))
			(line _yline _xline, lcolor(black) lpattern(dash))
			,
			graphregion(color(white))
			plotregion(margin(zero))
			xtitle("False Positive Rate")
			ytitle("True Positive Rate")
			legend(off)
			name(figureA2, replace)
			note("AUC: `auc'", pos(5) ring(0))
			;
		#delimit cr
		
	restore

	*save graph
	foreach g in figureA2{
		graph export ${figures}`g'.png, name(`g') replace width(3000)
	}
		

}

if `figureA3'{
	use ${data}indiv_WTP_with_shrink.dta, clear
	
	gen W_wtp_social = wtp_social_s  * 1000
	gen W_wtp_inperson = wtp_in_person_s  * 1000
	
	foreach f in social inperson{
		replace W_wtp_`f' = -4000 if W_wtp_`f'<-4000 & !missing(W_wtp_`f')
		replace W_wtp_`f' = 8000 if W_wtp_`f'>8000 & !missing(W_wtp_`f')
	}
	
	pwcorr W_wtp_social W_wtp_inperson, sig
	matrix sig = r(sig)
	local p = sig[1,2]
	local cor: display %4.3f r(rho)
	if `p'<0.1{
		local cor = "`cor'*"
	}
	if `p'<0.05{
		local cor = "`cor'*"
	}
	if `p'<0.01{
		local cor = "`cor'*"
	}

	#delimit ;
	twoway (scatter W_wtp_social W_wtp_inperson, msize(*0.4) mcolor(gs8%35))
		,
		graphregion(color(white))
		xtitle("WTP In-Person ($/year)")
		ytitle("WTP Social ($/year)")
		xline(0, lcolor(black))
		yline(0, lcolor(black))
		xlabel(-4000(2000)8000)
		ylabel(-4000(2000)8000)
		note("corr: `cor'   ", ring(0) pos(4))
		name("figureA3", replace)
		;
	#delimit cr
	
	*save graphs
	foreach g in "figureA3"{
		graph export ${figures}`g'.png, name(`g') replace width(3000) 
	}
}
