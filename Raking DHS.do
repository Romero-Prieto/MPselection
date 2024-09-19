clear
clear      matrix
clear      mata
tempfile   temp
tempfile   temp2
local      lISt           = "AO71 BJ71 BU71 KH81 CM71 ET71 GA71 GM81 HT71 MW7A ML7A MR71 NG7B NP7H PG71 RW81 SL7A SN7Z TZ7B TL71 ZM71 ZW72 LB7A UG7B ZA71" /*list of 25 surveys*/

local      pATh           = "/Users/lshjr3/Documents/DHS_MPselection/"           /*adjust path*/

/*run this code one time and then delete the las two caracters "* /" of this line. Then, all data preparation will be deactivated (green)*/
local      set            = "caseid v001 v002 v003 v005 v006 v007 v008 v011 v012 v016 v022 v023 v024 v025 v106 v169a v501" /*list of variables for identification*/
local      mm             = "mmidx_ mm1_ mm2_ mm3_ mm4_ mm5_ mm6_ mm7_ mm8_ mm9_ mm10_ mm11_ mm12_ mm13_ mm14_ mm15_ mm16_"
set        maxvar         6500

generate   survey         = ""
save     "`pATh'RakingData.dta", replace                                     /*this is the file and location of all prepared data*/

foreach survey of local lISt {                                               /*Loop for data preparation*/
	set        seed 1
	local      A              = substr("`survey'",1,2)                       /*country prefix*/
	local      B              = substr("`survey'",3,2)                       /*survey code*/
	use      "`pATh'IR/`A'IR`B'FL.DTA", clear                                /*opening of each individual recode, IR is the name of the folder*/
	keep      `set' mm*_*                                                    /*preserving relevant variables for identification and for adult mortality estimation*/
	
	generate   r              = runiform()                                   /*a–uniform–random variable to adjudicate an exact DOB to the respondent, given a month and year of birth*/
	generate   DOB            = mdy(v011 - floor((v011 - 1)/12)*12,1,floor((v011 - 1)/12) + 1900)*r + mdy(v011 +1 - floor(v011/12)*12,1,floor(v011/12) + 1900)*(1 - r)
	generate   W              = v005/1000000                                 /*weights*/
	sum        W 
	replace    W              = W/r(mean)		
	generate   respondent     = v003
	generate   cluster        = v001
	generate   household      = v002
	generate   strata         = v023 
	generate   Region         = v024
	generate   UR             = v025
	generate   Education      = v106
	recode     Education   (0 = 1) (8 = 1) (. = 1)                           /*to make a group of unknown, primary or less than primary education*/
	generate   Marital        = v501
	recode     Marital     (2 = 1) (3 = 0) (4 = 0) (5 = 0)
	generate   age            = v012
	generate   ageG           = 5*floor(age/5)
	generate   mobile         = v169a
	recode     mobile      (. = 0)
	generate   interview      = mdy(v006,v016,v007)                          	
	local      vAr            = "caseid interview respondent UR Region DOB W cluster household strata Region UR Education Marital age ageG mobile"
	save      `temp', replace
	
	contract   Region        
	generate   R              = _n
	keep       Region R
	save      `temp2', replace
	use       `temp', clear
	merge m:1  Region using `temp2', nogenerate noreport
	replace    Region         = R
	drop       R
	save      `temp', replace
		
	keep      `vAr'
	save      `temp2', replace
	use       `temp', clear
	
	reshape    long `mm', i(caseid) j(sibling) string                        /*wide to long, respondents without siblings are eliminated*/
	keep       caseid `mm' sibling interview
	rename     *_ *
	drop if    mmidx         == .
	
	egen       dateO          = min(interview)                               /*date of the first interview defined as the end of the Lexis diagram, but it can be adjusted if necessary*/	
	generate   dateA          = mdy(month(dateO),day(dateO),year(dateO) - 7) /*to define a 7-year Lexis diagram, but it can be adjusted if necessary*/
	generate   r              = runiform()                                   /*another–uniform–random variable to adjudicate an exact DOB to the siblings, given a month and year of birth*/
	generate   DOB_s          = mdy(mm4 - floor((mm4 - 1)/12)*12,1,floor((mm4 - 1)/12) + 1900)*r + mdy(mm4 +1 - floor(mm4/12)*12,1,floor(mm4/12) + 1900)*(1 - r) 
	replace    r              = runiform()                                   /*another–uniform–random variable to adjudicate an exact DOD to the siblings*/
	generate   DOD_s          = mdy(mm8 - floor((mm8 - 1)/12)*12,1,floor((mm8 - 1)/12) + 1900)*r + mdy(mm8 +1 - floor(mm8/12)*12,1,floor(mm8/12) + 1900)*(1 - r)
	replace    DOD_s          = max(DOB_s,DOD_s) if DOD_s != .               /*to prevent negative duration for siblings dying within the first month of life*/
	format     %tdDD/NN/CCYY interview DOB* DOD dateA dateO
	generate   age_D          = 1/365.25*(DOD_s - DOB_s)
	generate   exposure       = 0
	
	forvalues j = 15(5)55 {                                                  /*to calculate exposure and count events for all ages and for each 5-year age interval*/
		generate   exposure`j'  = 1/365.25*max(min(min(dateO,mdy(month(DOB_s),month(DOB_s),year(DOB_s) + `j' + 5)),DOD_s) - max(dateA,mdy(month(DOB_s),month(DOB_s),year(DOB_s) + `j')),0)
		generate   events`j'    = 1   if  (DOD_s >= dateA & DOD_s < dateO) & (age_D >= `j' & age_D < `j' + 5)
		replace    exposure     = exposure + exposure`j'
		recode     events`j'   (. = 0)
		}
		
	drop if    exposure      == 0                                            /*to delete sibling not exposed to the risk of dying*/
	drop       exposure
	generate   sexS           = mm1	
	drop if    mm2           == 8                                            /*to delete siblings with unknown survival status, assuming they were missed at random–not mortality selected*/
	drop if    mm1           == .                                            /*to delete siblings with unknown sex, assuming they were missed at random–not mortality selected*/
	keep       caseid DOB_s DOD_s sexS exposure* events*
	order      caseid DOB_s DOD_s sexS exposure* events*
	label drop _all
	bysort     caseid: generate   k = _n         
	bysort     caseid: generate   K = _N
	
	merge m:1  caseid using `temp2', nogenerate
	recode     k K         (. = 0)
	save      `temp', replace                                                /*to save the IR data*/
		
	use      "`pATh'PR/`A'PR`B'FL.DTA", clear                                /*opening of each household member recode, PR is the name of the folder*/
	keep       hhid hvidx hv001 hv002 hv003 hv005 hv023 hv101 hv102 hv103 hv104 hv105 hv106 hv107 hv112 hv114 hv201 hv206 hv215
	
	generate   cluster        = hv001
	generate   household      = hv002
	generate   respondent     = hvidx
	
	generate   strata         = hv023
	generate   sex            = hv104
	generate   Education      = hv106
	recode     Education   (0 = 1) (8 = 1) (. = 0)
	generate   age            = hv105
	recode     age         (. = 50)
	
	generate   Electricity    = hv206
	generate   Roofing        = 0
	recode     Roofing     (0 = 1)    if hv215 == 31 | hv215 == 33 | hv215 == 34 | hv215 == 35 | hv215 == 36 /*good material excluding wood*/
	generate   Water          = 0
	recode     Water       (0 = 1)    if hv201 == 11 | hv201 == 12 | hv201 == 13 | hv201 == 14 | hv201 == 21 | hv201 == 31 | hv201 == 41 | hv201 == 51 | hv201 == 62 | hv201 == 71 
	
	generate   usual_resident = 1     if hv102 == 1
	bysort     hhid: egen   HH_size = sum(usual_resident)
	generate   householdS     = 1     if HH_si  < 5
	replace    householdS     = 2     if HH_si  < 9  & householdS == .
	recode     householdS  (. = 3)
	
	generate   Wh             = hv005/1000000/HH_size
	generate   head           = 1     if hv101 == 1
	recode     head        (. = 0)
	recode     hv112 hv114 (. = 0)
	generate   with_parents   = max(min(1,hv112),min(1,hv114))
	keep       hhid hvidx cluster household respondent strata sex Education age head Electricity Roofing Water usual_resident with_parents householdS Wh
	save      `temp2', replace
	
	keep if    head          == 1
	keep       hhid age sex Education
	rename    (age sex Education) =_head

	merge 1:m  hhid using `temp2', nogenerate
	keep       hvidx cluster household respondent strata Electricity Roofing Water Wh head with_parents householdS usual_resident *_head hhid
	
	pca        Electricity Roofing Water [aw = Wh]                           /*to estimate the amenity index*/
	predict    P1                                                            /*prediction of the index*/
	sort       P1                                                            /*sort the population to calculate quintiles*/
	generate   w              = sum(Wh)                                      
	replace    w              = w/w[_N]
	local      Q              = 2
	generate   Q              = ceil(`Q'*w)
	drop       w
	tabulate   Q [aw = Wh]
	keep       cluster household respondent Q householdS usual_resident with_parents head *_head hhid Electricity Roofing Water
	clear      matrix
	save      `temp2', replace                                               /*to save the PR data*/
		
	use       `temp.dta', clear
	merge m:1  cluster household respondent using `temp2', nogenerate keep(master match)
	generate   survey         = "`survey'"
	merge m:1  survey using "`pATh'survey_names.dta", nogenerate noreport keep(master match)
	order      survey cluster household respondent hhid
	
	egen       GO             = cut(age), at(15,20,30,40,50,65) icodes
	replace    GO             = GO + 1
	replace    Electricity    = Electricity + 1

	
	append     using "`pATh'RakingData.dta"
	save     "`pATh'RakingData.dta", replace
	}
	*/

use      "`pATh'RakingData.dta", clear
keep       survey caseid K k
save     "`pATh'RakingOuTpUT.dta", replace

foreach survey of local lISt {                                               /*Loop for post-stratification*/
	dis      "`survey'" 
	use      "`pATh'RakingData.dta", clear
	keep if    survey        == "`survey'"
	save      `temp', replace
	keep if    k             <= 1                                            /*to contract the data to include only respondents*/
		
	generate   Total          = 1
	tabulate   Total [aw = W], matcell(Total)
	local      listA          = "GO UR Region householdS Education Electricity"
	keep       caseid W mobile age ageG `listA'
	
	local      variable       = ""
	local      totals         = "_cons = 10000"
	foreach var of local listA {
		local            variable       = "`variable'" + " i.`var'"
		tabulate        `var' [aw = W], matcell(`var')
		forvalues i = 1(1)`= rowsof(`var')' {
			local            number         = `var'[`i',1]/Total[1,1]*10000
			local            totals         = "`totals'" + " `i'.`var' = `number'"
			}	
		}

	sort             W
	local            up             = 10.0
	local            lw             = 00.1
	dis "`totals'"
	svycal     rake `variable' if mobile == 1 [pw = W], force generate(WR) totals(`totals') ul(`up') ll(`lw')
	sum        WR
	replace    WR             = WR/r(mean)
	keep       caseid WR
	save      `temp2', replace
	
	use       `temp', clear
	merge m:1  caseid using `temp2', nogenerate noreport
	merge 1:1  survey caseid K k using "`pATh'RakingOuTpUT.dta", nogenerate noreport keep(match using)
	save     "`pATh'RakingOuTpUT.dta", replace	
	}
export delimited using "`pATh'RakingOuTpUT.csv", replace

foreach survey of local lISt {
	use      "`pATh'RakingOuTpUT.dta", clear  
	keep if    survey        == "`survey'"
	dis      "`survey'"
	svyset     cluster [pw = W], strata(strata) || caseid	
	svy: logit mobile `variable'
	}
