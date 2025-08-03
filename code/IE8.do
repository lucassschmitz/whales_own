cd ..
use Data/masters_voyages_merged, clear

**# Bookmark #1 Clean data 
	
	keep yearout product duration  tonn_avg mastercode boatcode product_nom re_prod 

	// 5 and 10 year periods 
	*generate year10 = int(yearout / 10) * 10
	*generate year5 = int(yearout / 5) * 5

	replace duration = 1 if duration == 0 

	// gen production vars 
	*gen log_p = log(product) 
	*gen log_rp = log(re_prod)

	bys mastercode (yearout): gen n_voy = _n

	summ product, d
	summ re_prod, d 


	bysort mastercode (yearout): gen lag_prod = re_prod[_n-1]

	gen n_voy_top = n_voy 
	replace n_voy_top = 5 if n_voy_top > 5 // experience top-coded to max of 5 voyages 
	
	
	bys mastercode (n_voy): gen master_last = (_N == _n) 
	
	bys mastercode (yearout): gen exp_years = sum(duration) - duration
	drop if product == . 
	by mastercode: gen Total = _N
label var re_prod "Product"
 label var n_voy "Exp. voy"
label var lag_prod "Production lag"
label var exp_years "Exp. years"
label var master_last "Last voy."
label var tonn_avg "Weight"

/// Table 1 
order mastercode n_voy
xtset mastercode n_voy

gen re_prod2 = re_prod/1000
label var re_prod2 "Product" 

save Data/temps/temp, replace

est clear

eststo m1: logit master_last L.re_prod2

eststo m2: logit master_last re_prod2

eststo m3: logit master_last re_prod2 L.re_prod2

twoway (scatter master_last re_prod2) (lfit master_last re_prod2)


summ Total if n_voy ==1 

eststo m2_2: logit master_last re_prod2 if Total <= 3
eststo m2_3: logit master_last re_prod2 if Total == 4
eststo m2_4: logit master_last re_prod2 if Total > 4 


esttab  m1 m2  m3 m2_2 m2_3 m2_4   using "Writeup/Tables/IE8do_temp.tex", replace ///
    title("Logit ")    /// 
	mtitles("All" "All" "All" "<3" "4 voys" "<4")   nocon /// keep(lag_prod)
    label  b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    s(r2_p N )  /// 
	nogaps noomit  addnotes(" ")
	
//////////////// Coefficient over time //// /
use Data/temps/temp, clear
summ Total if n_voy == 1, d

// 1. Run regressions by n_voy, storing coefficient and its standard error
statsby coef=_b[re_prod2] se=_se[re_prod2], by(n_voy) nodots clear: ///
    regress master_last re_prod2

// 2. Compute 95% confidence‐interval bounds
gen lb = coef - invnormal(0.975)*se
gen ub = coef + invnormal(0.975)*se

// 3. Plot coefficient with vertical "whiskers" for the CI, restricting to n_voy<10
twoway ///
    (rcap ub lb n_voy if n_voy < 6, sort) ///  // vertical bars from lb to ub
    (line coef n_voy if n_voy < 6, sort lcolor(blue) lwidth(medium)), ///
    ylabel(, angle(horizontal)) ///
    xtitle("Number of Voyages (n\_voy)") ///
    ytitle("Coefficient on re\_prod2") ///
    title("β̂ and 95% CI for master\_last on re\_prod2 by n\_voy") note("93% of captains do less than 6 voyages.")

 graph export "Figures/meeting_250610/IE8_real_coeff_time.png", replace

 
 //////////////// Coefficient over time with interaction //// /
use Data/temps/temp, clear

bysort mastercode (n_voy): gen lag_prod2 = re_prod2[_n-1]

gen interaction = re_prod2 * lag_prod2

preserve 
statsby coef=_b[interaction] se=_se[interaction], by(n_voy) nodots clear: regress master_last interaction

gen lb = coef - invnormal(0.975)*se
gen ub = coef + invnormal(0.975)*se

twoway ///
    (rcap ub lb n_voy if n_voy < 6, sort) ///  // vertical bars from lb to ub
    (line coef n_voy if n_voy < 6, sort lcolor(blue) lwidth(medium)), ///
    ylabel(, angle(horizontal)) ///
    xtitle("Number of Voyages (n\_voy)") ///
    ytitle("Coefficient on interaction") ///
    title("β̂ and 95% CI for master\_last on re\_prod2 by n\_voy") note("No controls")

 graph export "Figures/meeting_250624/IE8_interaction_no_controls.png", replace
restore

preserve 
statsby coef=_b[interaction] se=_se[interaction], by(n_voy) nodots clear: ///
	regress master_last interaction lag_prod2 re_prod2

gen lb = coef - invnormal(0.975)*se
gen ub = coef + invnormal(0.975)*se

twoway ///
    (rcap ub lb n_voy if n_voy < 6, sort) ///  // vertical bars from lb to ub
    (line coef n_voy if n_voy < 6, sort lcolor(blue) lwidth(medium)), ///
    ylabel(, angle(horizontal)) ///
    xtitle("Number of Voyages (n\_voy)") ///
    ytitle("Coefficient on interaction") ///
    title("β̂ and 95% CI for master\_last on re\_prod2 by n\_voy") note("Controls for current and lag production")

 graph export "Figures/meeting_250624/IE8_interaction_controls.png", replace
restore
 
 
//////////
 
 use Data/temps/temp, clear
	
bysort mastercode (n_voy): egen y_sumcum = sum(re_prod2)
gen prev_y = y_sumcum - re_prod2

gen avg_y_prev = prev_y / (n_voy - 1) 
replace avg_y_prev = .  if n_voy == 1 	

egen FE_tercile = xtile(avg_y_prev) , by(n_voy) nq(2) 

// 1. Run regressions by n_voy, storing coefficient and its standard error
statsby coef=_b[re_prod2] se=_se[re_prod2], by(FE_tercile n_voy) nodots clear: ///
    regress master_last re_prod2 if n_voy < 7

	
twoway ///
    (line coef n_voy if FE_tercile == 1 , sort lcolor(blue)    lpattern(solid)) ///
    (line coef n_voy if FE_tercile == 2, sort lcolor(orange)  lpattern(longdash)), ///
    legend(order(1 "Lower skill" 2 "Higher skill" )) ///
    ylabel(, angle(horizontal)) ///
    xtitle("Time") ///
    ytitle("β̂_{t,q} (coef on y)") ///
    title("Coefficient on y by Time and shown quality")
 

graph export "Figures/meeting_250610/IE8_real_coeff_time_byquality.png", replace

  