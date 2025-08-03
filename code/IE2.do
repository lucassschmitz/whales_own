
cd ..
************************************************************* *
**# Bookmark #1
**** one possibility is that bigger/smaller ships produce different products, eg.. big ones spetialize on sperm and small ones on bones. Separate ships by their size and look at their shares. 
************************************************************* *
use "Data/temp", clear

summarize tonn_avg, d 
local p50 = r(p50)
gen big = 1 if tonn_avg >= `p50'
replace big = 0 if tonn_avg < `p50'
bysort big: summ tonn_avg, d

*normalize bone to 1. 
gen rat_oil = oil/bone 
gen rat_sperm = sperm/bone 


collapse rat_oil rat_sperm, by(yearout big) 

replace rat_oil = .  if rat_oil > 10 // outlier, not plausible  
replace rat_sperm = .  if rat_sperm > 6 // outlier, not plausible  

twoway (lpoly rat_oil yearout if big == 1 )  /// 
		(lpoly rat_oil yearout if big == 0 ) /// 
		, legend(order(1 2) label(1 "Big") label(2 "Small")) /// 
		title("Production ratio by ship size") xtitle("Year") ///
		ytitle("Ratio: oil/bone")

twoway (lpoly rat_sperm yearout if big == 1 )  /// 
		(lpoly rat_sperm yearout if big == 0 ) /// 
		, legend(order(1 2) label(1 "Big") label(2 "Small")) /// 
		title("Production ratio by ship size") xtitle("Year") ///
		ytitle("Ratio: sperm/bone")

		
************************************************************* *
**# Bookmark #2
/* Plot 
1. share and number of  ships that are schooners and Full rig a
2. avg and total tonnage per year. 
*/ 
*************************************************************
use "Data/temp", clear

gen full_rig = strpos(rig,"hip")>0
gen schooner = strpos(rig,"chr")>0

keep vesselid tonn_avg builtdate enddate life full_rig schooner 
bys vesselid (builtdate) : keep if _n==1    // 1 row per vessel

expand life +1                         // one row per ship-year
bys vesselid : gen year = built + _n - 1 // the calendar year

twoway (lpolyci full_rig year) (lpolyci schooner year), /// 
	title("Types of ships") xtitle("Year") ytitle("Share") /// 
	legend(order(2 3) label(2 "Full rig") label(3 "Schooner"))
graph export "Figures/IE2_sharetypeship_year.png", replace

drop if missing(enddate) | missing(builtdate) | missing(tonn_avg)

twoway lpolyci tonn_avg year,  title("Average tonnage") ///
		xtitle("Year") ytitle("Tonns") legend(off)
graph export "Figures/IE2_avgtonnage_year.png", replace

collapse (sum) tonn_avg full_rig schooner, by(year) 

twoway lpolyci tonn_avg year, title("Total tonnage" ) xtitle("Year") /// 
	ytitle("Tonns") legend(off)
graph export "Figures/IE2_totaltonnage_year.png", replace

twoway (lpolyci full_rig year) (lpolyci schooner year), /// 
	title("Types of ships") xtitle("Year") ytitle("Number") /// 
	legend(order(2 3) label(2 "Full rig") label(3 "Schooner"))
graph export "Figures/IE2_totaltypeship_year.png", replace
 

************************************************************* *
**# Bookmark #3
*/Tower (1907) mentions that at some point it was not profitable to use the ships, examine this possibility by plotting the use ratio for each year */ 
************************************************************* *
use "Data/temp", clear
drop if vesselid == "AS0000"

drop if missing(yearin) 
drop life duration 
by vesselid: gen life = enddate - builtdate + 1
gen duration = yearin - yearout + 1

expand duration 
bysort vesselid yearout yearin: gen year_count = yearout + _n - 1
bysort year_count: gen N_sea = _N

preserve 
by year_count: keep if _n == 1

keep year_count N_sea
tempfile at_sea
save "`at_sea'", replace
restore 

bysort vesselid: keep if _n == 1
expand life
drop year_count
bysort vesselid: gen year_count = builtdate + _n -1 
bysort year_count: gen N_available = _N 
by year_count: keep if _n == 1
keep year_count N_available
*existing boats 
mmerge year_count using "`at_sea'", type(1:1) 

gen use_ratio = N_sea/ N_available 


twoway lpolyci use_ratio year_count,  title("Yearly use ratio") ///
		xtitle("Year") ytitle("Ships at sea/Ships available") legend(off)
graph export "Figures/IE2_useratio_year.png", replace



************************************************************* *
**# Bookmark #4
/*At the firm level years with lot of entry also have lot of exit. Examine this possibility in the case of ship entry-exit. 
-> it is not the case. 
*/  
************************************************************* *
use "Data/temp", clear

drop if missing(yearin) 
bysort vesselid: keep if _n == 1
bysort builtdate: gen entry = _N  
bysort enddate: gen exit = _N  


keep builtdate enddate entry exit

preserve
tempfile ent
collapse   entry, by(builtdate)     // entrants per year
rename builtdate year
save `ent'
restore

collapse   exit,  by(enddate)     // exits per year
rename enddate year
merge 1:1 year using `ent', nogenerate
replace entry = 0 if missing(entry) // fill gaps

corr entry exit 
twoway (lpolyci entry year) (lpolyci exit year),  title("Yearly Entry/Exit") ///
		xtitle("Year") ytitle("Number of ships") legend(order(2 3) label(2 "Entry") label(3 "Exit"))
graph export "Figures/IE2_entryexit_year.png", replace

************************************************************* *
**# Bookmark #5
/*Distance by type of ship each year. 
*/  
************************************************************* *
 
use "Data/temp", clear

keep vesselid yearout yearin rig max_distance
gen full_rig = strpos(rig,"hip")>0
gen schooner = strpos(rig,"chr")>0
collapse max_distance, by(full_rig schooner yearout)

twoway (lpolyci max_distance yearout if full_rig == 1) /// 
		(lpolyci max_distance yearout if schooner == 1), /// 
		legend(order(2 3) label(2 "Full rig") label(3 "Schooner")) /// 
		title("Distance by type of ship") ytitle("Voyage distance (1000km)") xtitle("Years") 
graph export "Figures/IE2_distancebytypeship_year.png", replace
		

************************************************************* *
**# Bookmark #6
/* Ratio of product over time 
-> striking that bones do not increase more, since they got much more expensive.  
-> motivated by the Tower (1907, p. 72) mentions that before the increase in bone prices, bones were thrown away
*/  
**************************************************************

use "Data/temp.dta", clear

gen ratio1 = bone/oil 
gen ratio2 = bone/sperm 
gen ratio3 = log(bone/oil) 
gen ratio4 = log(bone/sperm) 


collapse (mean) ratio*, by(yearout)

twoway (lpolyci ratio1 yearout) (lpolyci ratio2 yearout), /// 
		legend(order(2 3) label(2 "Bones/oil") label(3 "Bones/sperm")) /// 
		title("Evolution of products") ytitle("Ratio of products") xtitle("Years") 
graph export "Figures/IE2_ratioproducts_year.png", replace
		

twoway (lpolyci ratio3 yearout) (lpolyci ratio4 yearout), /// 
		legend(order(2 3) label(2 "log(Bones/oil)") label(3 "log(Bones/sperm)")) /// 
		title("Evolution of products") ytitle("Ratio of products") xtitle("Years") 
graph export "Figures/IE2_logratioproducts_year.png", replace
		

		
	
************************************************************* *
**# Bookmark #7
/*Tower 1907 (p.73) writes "in 1871 the entire Arctic fleet of thirty-four vessels was completely destroyed by pack ice". this is not visible in our data 1871 and 1872 have only 13 and 10 ships with a nonmissing returncode. Which is far from 34 even if we count returncodes like "So" that is sold in foreing port. 
*/  
************************************************************* *

use "Data/temp", clear

bysort vesselid: keep if _n == 1 
tab returncode, missing 
keep if inlist(returncode, "C", "L", "So") // So: sold foreign port 

contract yearout returncode

twoway 	(lpoly _freq yearout if return == "C" ) /// 
		(lpoly _freq yearout if return == "L" ) /// 
		 , ytitle("Fate") xtitle("Year of departure") ///
         title("Voyage outcomes ") ///
         legend(order(1 "Condemned" 2 "Lost/sank/wrecked" )) 
		 
graph export "Figures/IE2_vesselfate_year.png", replace
