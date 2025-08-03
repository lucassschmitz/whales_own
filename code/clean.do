cd ".."

//////////////////////////////////////////////////////
**# Bookmark #1 Clean voyages 
//////////////////////////////////////////////////////
clear all
set more off

insheet using ../Whales/Data/AmericanOffshoreWhalingVoyages/AmericanOffshoreWhalingVoyages/voyages_20211020.txt, tab clear
destring sperm, force replace
destring bone, force replace
destring oil, force replace
egen voyagecode = group(voyageid)
egen boatcode = group(vesselid)

drop if voyagerank > 1


mmerge voyageid using "../Whales/Data/voyage_distances.dta", type(n:1) unmatched(master)

duplicates drop
gen years = yearin - yearout
drop if years < 0
keep if yearout >= 1790 & yearout <= 1910

* up to this point the cleaning follows aggregates.dta

gen year = yearout //for merging 
merge m:1 year using "../Whales/Data/prices.dta", gen(_merge2) // A: relevant prices for decisions are departure prices since captains observe them upon departure. 
drop year
drop if _merge2 == 2 

gen nonmissing = !missing(oil) + !missing(sperm) + !missing(bone) // A: if all are missing data problem, if only some are missing then they are 0 
foreach var in oil sperm bone {
    replace `var' = 0 if missing(`var') & nonmissing > 0
}
drop nonmissing

gen product = bone * price_bone_real + 42* oil * price_oil_real + 42* sperm * price_sperm_real //oil, sperm quantities are in barrels and prices in gallons  

gen product_nom = bone * price_bone + 42* oil * price_oil + 42* sperm * price_sperm //oil, sperm quantities are in barrels and prices in gallons 

/// clean tonnage var 
	replace tonnage = trim(tonnage)
	replace tonnage = subinstr(tonnage, " or ", "/", .)
	replace tonnage = subinstr(tonnage, " ?", "", .)
	replace tonnage = subinstr(tonnage, "?", "", .)
	split tonnage, parse("/") gen(ton_)
	destring  ton_*, replace

	foreach v of varlist ton_1-ton_5 {
		replace `v' = . if `v' > 1000
	}
	egen tonn_avg = rowmean(ton_1 ton_2 ton_3 ton_4 ton_5)
	drop ton_*


/// generate the entry and exit dates 
	gen aux1 = min(yearin, yearout) 
	gen aux2 = max(yearin, yearout)
	 
	bysort vesselid: egen built = min(aux1) // earliest year there are records 
	bysort vesselid: egen scrap = max(aux2) // latest year with records. robust to having trips with no 'yearin' recorded 
	drop aux* 
	corr built builtdate
	replace builtdate = built if missing(builtdate) |  (built < builtdate) // if there is a trip later than builtdate, there was an error 
	scatter built builtdate

	gen str cleaned = ustrregexra(end,"[^0-9]"," ")
	split cleaned, parse(" ") destring gen(aux)
	egen enddate = rowmax(aux*)
	replace enddate = . if enddate < 1700
	
	replace enddate = scrap if missing(enddate) | (scrap > enddate) // 
	scatter enddate scrap // !! enddates much greater than their last trip. 
	
	by vesselid: gen life = enddate - builtdate if _n == 1
	gen duration = yearin - yearout 
	drop cleaned aux* scrap built
	

//// remove outliers in production 
	histogram bone 
	summarize bone, d // reasonable distribution 

	histogram sperm 
	summarize sperm, d // reasonable distribution 

	histogram oil 
	summarize oil, d // p99 == 4186 and max() == 1.7M 
	drop if (oil > 1.5* r(p99)) & !missing(oil)


/// clean returncode 
/* returncode indicates the fate of the vessel, but had mainly missings, we use 'end' to fill missing values. 

The data has inconsistencies, for example end = 'wrecked off' should be returncode = 'L' (vesselid AS0656) but for AS0735 is A. Here we will state the inconsistencies and for each inconsistency we provide a vesselid -6character code-as an example of the problem we are mentioning. 

Inconsistencies: 
- 'foundered'  sometimes is L (AS0624) and sometimes A (AS0644) 'withdr' is classified as So (AS0902) or L (AS0833) 
- 'broken up' sometimes is classified as C (AS0273, AS2759) or So (AS0312) we choose C
- obs with 'wreck' are classifed as C or L (see  AS0065 and AS0013)  


Moreover, some observations could be assigned different categories. We try to follow the original coding, we use the following criteria:  
- obs with 'cond' and 'sold' are classificed as C (see AS1367) 
- obs with 'cond' and 'wrecked' are classified as C (see  AS0462 or AS0631) 


We did not classify the following values for 'end': stone fleet (??), withdr. 
 
Since the coding is inconsistent we will try to be consistent but we are arbitrary, for example we classify 'foundered' as L and not A 

Code
A	Abandoned
B	Burned
C	Condemned
F	Chartered for freighting
L	Lost, sank, wrecked
M	Missing
S	Seized
So	Sold, in a foreign port
?, prob	Uncertain

*/

	* rc expands returncode to all obs. within vesselid
	gen rc = returncode    
	sort vesselid rc
	by  vesselid: replace rc = rc[_N]    if missing(rc)             
	bysort vesselid (end): assert end == end[1] // end doesnt change within vesselid 

	gen  aux = lower(end)

	gen f_burn   = regexm(aux,"burn|fire")
	gen f_sold   = regexm(aux,"sold")
	gen f_aband  = regexm(aux,"aband")
	gen f_lost   = regexm(aux,"wrecked|crushed|lost|sunk|sank|capsize|blown up|ice|foundered") // 'wreck' classifies as lost (see AS0013)
	gen f_cond   = regexm(aux,"cond|broken up")

	gen rc2 = "So" if  f_sold  
	replace rc2 = "B" if f_burn
	replace rc2 = "C"  if  f_cond  
	replace rc2 = "A"  if  f_aband 
	replace rc2 = "L"  if  f_lost  

	replace rc2 = rc if missing(rc2)
	replace returncode = rc2
	drop f_* aux rc*

	save "Data/temp.dta", replace
	
	
//////////////////////////////////////////////////////
**# Bookmark #2 Clean masters 
//////////////////////////////////////////////////////


clear all
set more off

insheet using ../Whales/Data/AmericanOffshoreWhalingCrewlists/AmericanOffshoreWhalingCrewlists/crewentries_20200302.csv, comma clear

* clean lay
split lay, p("-")
replace lay2 = lay3 if !missing(lay3) 
drop lay3 
destring lay1, force replace
destring lay2, force replace
drop lay
gen lay = lay1/lay2
drop lay1 lay2 

replace citizenship = lower(citizenship)
replace rank = lower(rank) 
keep if rank == "master"

bysort voyage name_first name_last: keep if _n == 1  

save "Data/masters.dta", replace

//////////////////////////////////////////////////////
**# Bookmark #3 merge data 
//////////////////////////////////////////////////////

/// merge masters with their voyages. 
use Data/masters, clear
bysort voyageid (lay): keep if _n == 1 
keep voyageid lay age
isid voyageid
tempfile temporary
save "`temporary'", replace

use Data/temp, clear
egen mastercode = group(masterid)
sort mastercode
drop if inlist(masterid, "AM0000")
drop if voyageid == "AV19340" & voyagename == "Unknown : 1775-" // duplicated voyage 
mmerge voyageid using "`temporary'", type(1:1) 
keep if _merge == 3

bysort masterid: gen Total = _N 
bysort masterid (yearout): gen experience = _n 

replace duration = 1 if (duration == 0) & (yearin == yearout)
gen re_prod = product/duration
bysort mastercode (yearout): gen n_voyage = _n 
bysort mastercode (n_voyage): gen l_re_prod = re_prod[_n-1]

bysort mastercode (n_voyage): gen change_boat = 1 if (boatcode != boatcode[_n-1]) & _n > 1
bysort mastercode (n_voyage): gen stay_boat = 1 if (boatcode == boatcode[_n-1]) & _n > 1 


save Data/masters_voyages_merged, replace
export delimited using "Data/masters_voyages_merged.csv", replace


