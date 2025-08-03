cd ..
////////////////////////////////
**# Bookmark #1 dist of # of voyages and weight/revenue by voyage number. 
////////////////////////////////
use Data/temp, clear

egen mastercode = group(masterid)
sort mastercode
drop if inlist(masterid, "AM0000")

bysort masterid: gen Total = _N 
bysort masterid (yearout): gen experience = _n -1


keep mastercode masterid product yearin yearout Total experience tonn_avg
tab Total if experience == 1 

histogram Total if experience == 1, ytitle("Density") /// 
	xtitle("Total # of trips")	title("Trips by masters")

graph export "Figures/M1_master_voyages.png", replace


* hazard_ratio
preserve
bysort experience: gen N = _N
bysort experience: keep if _n == 1
keep experience N
sort experience
gen N_next = N[_n + 1]
gen hazard_ratio = (N - N_next) / N
list experience N N_next hazard_ratio
twoway lpoly hazard_ratio experience, title("Captain hazard rate") /// 
	ytitle("Share of last trip") xtitle("Trip number")
graph export "Figures/M1_hazardratio.png", replace
restore


drop if experience > 9
twoway (lpoly tonn_avg experience, yaxis(1)) ///
       (lpoly product experience, yaxis(2)), ///
       ytitle("Ship Weigth (Tonnes)", axis(1)) ///
       ytitle("Voyage Revenue", axis(2)) ///
	   xtitle("# of prior voyages") /// 
       legend(label(1 "Weight") label(2 "Revenue")) /// 
	   note("Sample: up to the 10th voyage of every captain")
graph export "Figures/M1_rvenue&weight_by_experience.png", replace

drop if experience > 4 
twoway (lpoly tonn_avg experience, yaxis(1)) ///
       (lpoly product experience, yaxis(2)), ///
       ytitle("Ship Weigth (Tonnes)", axis(1)) ///
       ytitle("Voyage Revenue", axis(2)) ///
	   xtitle("# of prior voyages") /// 
       legend(label(1 "Weight") label(2 "Revenue")) /// 
	   note("Sample: up to the 5th voyage of every captain")
graph export "Figures/M1_rvenue&weight_by_experience(2).png", replace



////////////////////////////////
**# Bookmark #2 revenue/weight depending on past performance. 
////////////////////////////////
 
use Data/temp, clear

egen mastercode = group(masterid)
sort mastercode
drop if inlist(masterid, "AM0000")
bysort masterid: gen Total = _N 
bysort masterid (yearout): gen experience = _n 
keep mastercode  product yearin yearout Total experience tonn_avg
keep if Total > 1 & experience < 5

// for each trip define whether the production is below or above the mean 
bysort experience: egen avg_pro = mean(product) 
gen high = (product > 1.2*avg_pro) if  !missing(product)
replace  high = 0 if (product <= .8*avg_pro &  !missing(product)) 
bysort mastercode (experience): gen next_prod = product[_n+1] 
bysort mastercode (experience): gen next_tonn = tonn_avg[_n+1]


twoway (lpoly next_tonn experience if high == 1) /// 
		(lpoly next_tonn experience if high == 0) ///
		, ytitle("Ship Weigth (Tonnes)") title("Weight of next match") /// 
		xtitle("# of prior voyages") legend(label(1 "High") label(2 "Low")) /// 
	   note("Sample: up to the 5th voyage of every captain." ///
     "High(low): dummy indicating voyage produced more(less) than 1.2(0.8)x the average.")
graph export "Figures/M1_path_dependence.png", replace

		 /////

		
*gen diff = (next_prod-product)/product
*twoway (lpoly diff experience if high == 1) (lpoly diff experience if high == 0)


////////////////////////////////
**# Bookmark #3  captain's age.
////////////////////////////////
use Data/masters, clear
bysort voyageid (lay): keep if _n == 1 
keep voyageid lay age
isid voyageid
tempfile temporary
save "`temporary'", replace

use Data/temp, clear
drop if voyageid == "AV19340" & voyagename == "Unknown : 1775-"
mmerge voyageid using "`temporary'", type(1:1) 
keep if _merge == 3

egen mastercode = group(masterid)
sort mastercode
drop if inlist(masterid, "AM0000")
bysort masterid: gen Total = _N 
bysort masterid (yearout): gen experience = _n 

/// descriptives. 
destring age, replace
histogram age  , ytitle("Density") xtitle("Age") title("Captain's age")
graph export "Figures/M1_master_agedistribution.png", replace

gen exp_aux = experience 
replace exp_aux = 6 if experience > 6 
twoway (kdensity age if exp_aux == 1) (kdensity age if exp_aux == 5)

levelsof exp_aux, local(xvals)
local plotcmds
local i = 1
foreach x of local xvals {
    local plotcmds `plotcmds' (kdensity age if exp_aux == `x', ///
        lwidth(medthin)   ///
        legend(label(`i' "Voyage: `x'")))
    local ++i
}
twoway `plotcmds', ///
    xtitle("Age") ytitle("Density") title("Age Dist. by Experience Level")
graph export "Figures/M1_agedist_by_experience.png", replace
	
	
summ age, d
keep if age >= r(p5) & age <= r(p95)
  

twoway (lpolyci tonn_avg age, yaxis(1) color(blue%20)) /// 
		(lpolyci product age, yaxis(2) color(red%20)) /// 
		, ytitle("Ship Weigth (Tonnes)", axis(1)) ///
       ytitle("Voyage Revenue", axis(2)) /// 
		xtitle("Captain's age") ///	
       legend(order(1 3) label(1 "Weight(tonns)") label(3 "Revenue"))  

	   graph export "Figures/M1_outcomes_age.png", replace
		
			
/// keep mastercode age Total experience product tonn_avg voyagename