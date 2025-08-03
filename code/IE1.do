cd ..
************************************************************* *
**** Longer trips (duration-distance) produce more 
************************************************************* *
use "Data/temp", clear

drop if vesselid == "AS0000"

corr max_distance years
reghdfe product years, absorb(vesselid yearout)
reghdfe product max_distance, absorb(vesselid yearout)

reghdfe bone years, absorb(vesselid yearout)
reghdfe bone max_distance, absorb(vesselid yearout)
reghdfe oil years, absorb(vesselid yearout)
reghdfe oil max_distance, absorb(vesselid yearout)
reghdfe sperm years, absorb(vesselid yearout)
reghdfe sperm max_distance, absorb(vesselid yearout)

egen dist_bin = cut(max_dist), at(0(1)50)

preserve 
collapse (mean) bone sperm oil, by(years)
twoway ///
    (scatter bone years,      mcolor(blue) yaxis(2)) ///
	 (lfit bone years,      lcolor(blue) yaxis(2) ) ///
	(scatter oil years,      mcolor(orange) ) ///
	(lfit oil years,      lcolor(orange)  ) ///
	(scatter sperm years,      mcolor(green) ) ///
	(lfit sperm years,      lcolor(green)  ) ///
    , ///
	ytitle("Oil / Sperm",  axis(1)) ///
    ytitle("Bone",         axis(2)) ///
    xtitle("Years") ///
    legend( order(1 3 5)  label(1 "Bone")   label(3 "Oil") label(5 "Sperm")) /// 
	title("Production vs. distance")
graph export "Figures/IE_production_duration_voyages.png", replace
restore

preserve
collapse (mean) bone sperm oil, by(dist_bin)
twoway ///
    (scatter bone dist_bin, mcolor(blue)   yaxis(2)) ///
    (lfit    bone dist_bin, lcolor(blue)   yaxis(2)  legend(off)) ///
    (scatter oil  dist_bin, mcolor(orange)) ///
    (lfit    oil  dist_bin, lcolor(orange) legend(off)) ///
    (scatter sperm dist_bin, mcolor(green)) ///
    (lfit    sperm dist_bin, lcolor(green) legend(off)) ///
    , ///
    ytitle("Oil / Sperm", axis(1)) ///
    ytitle("Bone",        axis(2)) ///
    xtitle("Distance") ///
    legend(order(1 3 5) label(1 "Bone") label(3 "Oil") label(5 "Sperm")) ///
	title("Production vs. distance")
graph export "Figures/IE_production_distance_voyages.png", replace

restore


************************************************************* *
**** Ships optimize the ratio of products (bone/oil/sperm) 
************************************************************* *
use "Data/temp.dta", clear
 
drop if vesselid == "AS0000"

gen ratio1 = log(bone/oil) 
gen ratio2 = log(sperm/oil) 
gen ratio1p = log(price_bone_real/price_oil_real)
gen ratio2p = log( price_sperm_real/price_oil_real)

summ ratio1, d 
summ ratio2, d 

reghdfe ratio1 ratio1p, absorb(vesselid) 
reghdfe ratio2 ratio2p , absorb(vesselid)


************************************************************* *
**** Considerable variation in service life. Ships are idle low share of time
/* We create new built/scraped var because if there are missing voyages we would be underestimating the use-ratio. 
*/ 
************************************************************* *
use "Data/temp.dta", clear

drop if vesselid == "AS0000"

histogram life, xtitle("Years") title("Distribution of ship's service life") ///
	note("Service life: number of years between its launch and its final recorded voyage.")
graph export "Figures/IE_hist_servicelife.png", replace

bys vesselid : egen aux_built = min(yearout) 
by vesselid: egen aux_end = max(yearin)
by vesselid: gen aux_life = aux_end - aux_built if _n == 1 
replace duration = 1 if yearin == yearout // to consider short trips 

by vesselid: egen use = sum(duration) 
bys vesselid life: gen use_ratio = use/aux_life if _n == 1
replace use_ratio = 1 if use_ratio > 1 & !missing(use_ratio)

histogram use_ratio, xtitle("Years at sea/Service life") title("Distribution of use ratio") 
graph export "Figures/IE_hist_useratio.png", replace


************************************************************* *
**** Bigger ships do longer (in time and distance) voyages
************************************************************* *
use "Data/temp.dta", clear

egen ton_bin = cut(tonn_avg), at(0(50)850) icodes   
gen ton_mid = 25 + 50*ton_bin                      

collapse (mean) years max_dist, by(ton_mid)

twoway  (scatter years ton_mid, msymbol(O) msize(med) mcolor(blue)) ///
       (lfit    years ton_mid, lcolor(red)  lwidth(medthick)), ///
       ytitle("Years")  xtitle("Tonnage") ///
       title("Trip Duration vs. Vessel Size") ///
       legend(off)	   
graph export "Figures/IE_time_tonnage.png", replace

twoway  (scatter max_dist ton_mid, msymbol(O) msize(med) mcolor(blue)) ///
       (lfit    max_dist ton_mid, lcolor(red)  lwidth(medthick))  , ///
       ytitle("Trip Distance") xtitle("Tonnage") ///
       title("Trip Max. Distance vs. Vessel Size") ///
       legend(off)
graph export "Figures/IE_distance_tonnage.png", replace
	   

************************************************************* *
**** Product has no relation with lay
************************************************************* *
use ../Whales/Data/all_master_lays, clear
drop if missing(voyageid)
duplicates drop voyageid, force
tempfile clean_lays
save "`clean_lays'", replace

use "Data/temp.dta", clear
drop if missing(voyageid)
duplicates drop voyageid, force
mmerge voyageid using "`clean_lays'", type(1:1) 

keep if product < 7000
reghdfe product lay, absorb(vesselid)


************************************************************* *
**** Ships are being scrapped later over time, and they are scraped later if prices are high 
************************************************************* *
use "Data/temp.dta", clear

*price index: revenue of ship with avg. production 
summarize oil
scalar w_oil = r(mean)
quietly summarize bone
scalar w_bone = r(mean)
quietly summarize sperm
scalar w_sperm = r(mean)

gen price_index = w_oil * price_oil_real  + w_bone * price_bone_real + w_sperm * price_sperm_real
 
sort vesselid yearout
by vesselid: gen trips = _N  

by vesselid: egen use = sum(duration) 
by vesselid: gen use_ratio = use/life if _n == 1

collapse (mean) life use_ratio price_index  price_oil_real price_bone_real price_sperm_real bone sperm oil , by(enddate)

* when prices are high the ships are sold as scrap are older. 
twoway  (scatter life price_index, msymbol(O) msize(med) mcolor(blue)) ///
       (lfit    life price_index, lcolor(red)  lwidth(medthick))  , ///
       ytitle("Service life") xtitle("Price Index") legend(off) ///
	   note("Price index: real yearly prices multiplied by the average production")
graph export "Figures/IE_life_priceindex.png", replace

reg life price_index

keep if enddate >= 1810

twoway ///
    (scatter life enddate,  mcolor(blue)  yaxis(1)) ///
    (lowess  price_index enddate, lcolor(red) lwidth(medthick) yaxis(2)) ///
    , ///
    xtitle("Scrap year") ///
    ytitle("Average service life (years)", axis(1)) ///
    ytitle("Price index",                   axis(2)) ///
    legend(order(1 "Service life" 2 "Price index")) ///
    title("Service life and price index over time")

graph export "Figures/IE_life&price_scrapyear.png", replace

// ships scrapped each year get older over time
twoway  (scatter life enddate , msymbol(O) msize(med) mcolor(blue)) ///
       (lfit    life enddate, lcolor(red)  lwidth(medthick))  , ///
       ytitle("Service life") xtitle("Scrap year") ///
       legend(off)

graph export "Figures/IE_life_scrapyear.png", replace


************************************************************* *
**** REvenue share of bones is increasing over time and reaches almost 1
************************************************************* *
use "Data/temp.dta", clear

gen rev_oil = 42*w_oil * price_oil_real 
gen rev_bone = w_bone * price_bone_real 
gen rev_sperm = 42*w_sperm * price_sperm_real

gen rev_total = rev_oil + rev_bone + rev_sperm

gen share_oil = rev_oil / rev_total 
gen share_sperm = rev_sperm / rev_total 
gen share_bone = rev_bone / rev_total 

collapse (mean) share_* rev_* price_oil_real price_bone_real price_sperm_real bone sperm oil, by(yearout)


* when prices are high the ships sold as scrap are older. 
twoway  (scatter share_oil yearout)  (scatter share_bone yearout) /// 
		(scatter share_sperm yearout), /// 
		xtitle("Year departure") ytitle("Share") /// 
		legend(label(1 "Share Oil") label(2 "Share Bone") label(3 "Share Sperm")) ///
		title("Share of revenue by product") 
		
graph export "Figures/IE_productsshare_year.png", replace

		
		   