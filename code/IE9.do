cd ..
use Data/masters_voyages_merged, clear

/* Our ML algorithm requires 
- captainID 
- X_mat: boat characteristics 
- catch of each product 
- positive production dummies 
- duration 
- voyage identifier. 
*/ 


keep  yearout bone sperm oil  duration  tonn_avg mastercode 
replace duration = 1 if duration == 0 
bys mastercode (yearout): gen n_voy = _n
drop yearout 
sort mastercode n_voy

generate voyageID = _n

correlate bone sperm oil

* rename vars 
rename bone    Y1
rename sperm   Y2
rename oil     Y3
rename tonn_avg X1
rename duration Tau

* 4. Reshape from wide (one row per voyage) to long (three rows per voyage)
reshape long Y, i(voyageID) j(productID)
replace Y = 0 if missing(Y)
gen d = 1 if Y > 0 
replace d = 0 if Y == 0 


/// fill the missing 
	foreach v of varlist _all {
		count if missing(`v')
		display "`v': " r(N) " missing"
	}

	bysort mastercode: egen mean_X1  = mean(X1)
	bysort mastercode: egen mean_Tau = mean(Tau)
	replace X1  = mean_X1  if missing(X1)
	replace Tau = mean_Tau if missing(Tau)
	drop mean_X1 mean_Tau

	*still X1 missing, for captains that have no X1 observation. We can not identify the captain random coefficient without any info about weight, we drop this observations 
	drop if missing(X1)
	drop if missing(Tau)

egen captainID = group(mastercode) 

export excel using "Data/clean_ML_estimation.xlsx", firstrow(variables) replace



// Moments to target in the simulation 
summ Tau, d // mean 2.5, s.d 1.5, min 1 

summ X1 // mean 272, sd 103, min 45 

bysort productID: egen share_d_pos = mean(d>0) // 1 is .41, 2 is .81 and 3 is .69
drop share_d_pos

bysort productID: egen avg_Y = mean(Y) if Y> 0 // 1 is 17k, 2 is 700 and 3 is 1239



