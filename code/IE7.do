cd ..
use Data/masters_voyages_merged, clear

**# Bookmark #1 Clean data 
	
	keep yearout product duration  tonn_avg mastercode boatcode product_nom

	// 5 and 10 year periods 
	generate year10 = int(yearout / 10) * 10
	generate year5 = int(yearout / 5) * 5


	replace duration = 1 if duration == 0 

	// gen production vars 
	gen re_prod = product/duration 
	gen re_prod_nom = product_nom/duration
	gen log_p = log(product) 
	gen log_rp = log(re_prod)

	bys mastercode (yearout): gen n_voy = _n

	summ product, d
	summ re_prod, d 


	bysort mastercode (yearout): gen lag_prod = re_prod[_n-1]

	gen n_voy_top = n_voy 
	replace n_voy_top = 5 if n_voy_top > 5 // experience top-coded to max of 5 voyages 
	
	
	by mastercode: egen max_year = max(yearout)
	gen master_last = (yearout == max_year)
	drop max_year
	
	bys mastercode (yearout): gen exp_years = sum(duration) - duration
	drop if product == . 
label var re_prod "Product"
 label var n_voy "Exp. voy"
label var lag_prod "Production lag"
label var exp_years "Exp. years"
label var master_last "Last voy."
label var tonn_avg "Weight"

/// Table 1 
* Was not possible to include master fixed effects. Included only weight. 


eststo m1: reghdfe re_prod n_voy,      absorb(mastercode boatcode yearout)
eststo m2: reghdfe re_prod exp_years,       absorb(mastercode boatcode yearout)
eststo m3: reghdfe re_prod n_voy lag_prod,      absorb(mastercode boatcode yearout)
eststo m4: reghdfe re_prod exp_years lag_prod,       absorb(mastercode boatcode yearout)

eststo m5: reghdfe re_prod n_voy lag_prod master_last,      absorb(mastercode boatcode yearout)
eststo m6: reghdfe re_prod exp_years lag_prod master_last,       absorb(mastercode boatcode yearout)
 
esttab m1 m2 m3 m4 m5 m6 using "Writeup/Tables/IE7do_tab1.tex", replace ///
    title("Effect of Experience on Production")        /// mtitles("(1)" "(2)" "(3)" "(4)" "(5)") 
    label keep(n_voy exp_years lag_prod master_last)     /// 
    b(3) se(3)                                                 ///
    star(* 0.10 ** 0.05 *** 0.01)                   ///
    nogaps noomit                                       /// nolabel
    addnotes("All specifications include ship, captain and year fixed effects")
	

	
	
// Table 2: Arellano-Bond 

		
xtset mastercode n_voy

eststo m7: xtabond2 re_prod L.re_prod n_voy, gmm(L.re_prod, lag(2 .) collapse)  iv(n_voy) twostep robust   
  
eststo m7b: xtabond2 re_prod L.re_prod exp_years, gmm(L.re_prod, lag(2 .) collapse)  iv(exp_years) twostep robust   

eststo m7c: xtabond2 re_prod L.re_prod n_voy, gmm(L.re_prod, eq(d))  iv(n_voy) twostep robust 
eststo m7d: xtabond2 re_prod L.re_prod exp_years, gmm(L.re_prod, eq(d))  iv(exp_years) twostep robust   

	
eststo m8: xtabond2 re_prod L.re_prod n_voy tonn_avg, gmm(L.re_prod, lag(2 .) collapse)  iv(n_voy tonn_avg) twostep robust 

eststo m8b: xtabond2 re_prod L.re_prod exp_years tonn_avg, gmm(L.re_prod, lag(2 .) collapse)  iv(exp_years tonn_avg) twostep robust 

eststo m9: xtabond2 re_prod L.re_prod n_voy tonn_avg i.boatcode, gmm(L.re_prod, lag(2 .) collapse)  iv(n_voy tonn_avg i.boatcode)  robust //twostep 

eststo m9b: xtabond2 re_prod L.re_prod exp_years tonn_avg i.boatcode, gmm(L.re_prod, lag(2 .) collapse)  iv(exp_years tonn_avg i.boatcode)  robust // twostep
 
eststo m9c: xtabond2 re_prod L.re_prod n_voy tonn_avg i.boatcode, gmm(L.re_prod, eq(d))  iv(n_voy tonn_avg i.boatcode) twostep robust 

gen not_change = (boatcode == L.boatcode) 

eststo m10: xtabond2 re_prod L.re_prod n_voy tonn_avg if not_change == 1, gmm(L.re_prod, lag(2 .) collapse)  iv(n_voy tonn_avg) twostep robust 

eststo m10b: xtabond2 re_prod L.re_prod exp_years tonn_avg if not_change == 1, gmm(L.re_prod, lag(2 .) collapse)  iv(exp_years tonn_avg) twostep robust 

eststo m11: xtabond2 re_prod L.re_prod n_voy tonn_avg master_last, gmm(L.re_prod, lag(2 .) collapse)  iv(n_voy tonn_avg) twostep robust 

eststo m11b: xtabond2 re_prod L.re_prod exp_years tonn_avg master_last, gmm(L.re_prod, lag(2 .) collapse)  iv(exp_years tonn_avg) twostep robust 


esttab  m8 m8b m9 m9b m10 m10b m11 m11b  using "Writeup/Tables/IE7do_tab2.tex", replace ///
    title("Effect of Experience on Production")    /// mtitles("(1)" "(2)" "(3)" "(4)" "(5)")
    label keep( L.re_prod n_voy exp_years tonn_avg master_last)     ///   
    b(2) se(2)  width("2")                                       ///
    star(* 0.10 ** 0.05 *** 0.01)                              ///
    nogaps noomit                                       /// nolabel
    addnotes(" ")

	
/// Table 3 
est clear
bysort mastercode:  gen Total = _N


eststo m20: reghdfe re_prod lag_prod
eststo m21: reghdfe re_prod lag_prod, absorb(mastercode )
estadd local captain "Yes"

eststo m22: reghdfe re_prod lag_prod, absorb(boatcode )
estadd local ship "Yes"

eststo m23: reghdfe re_prod lag_prod, absorb(mastercode boatcode yearout)
estadd local captain "Yes"
estadd local ship "Yes"
estadd local time "Yes"


eststo m24: reghdfe re_prod lag_prod if Total == 3
eststo m25: reghdfe re_prod lag_prod if Total == 4  
eststo m26: reghdfe re_prod lag_prod if Total >= 5  


esttab  m20 m21 m22 m23  m24 m25 m26 using "Writeup/Tables/IE7do_tab3.tex", replace ///
    title("Effect of Experience on Production")    /// 
	mtitles("Product" "" "" "" "3 voy" "4 voy" "5+ voy") ///
    label keep(lag_prod) b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    s(N captain ship time, label("N" "Ship FE" "Captain FE" "Time FE"))  /// 
	nogaps noomit  addnotes(" ")

/// Table 4


eststo m31: reghdfe re_prod lag_prod if Total == 3, absorb(mastercode )
estadd local captain "Yes"
eststo m32: reghdfe re_prod lag_prod if Total == 4, absorb(mastercode )
estadd local captain "Yes"
eststo m33: reghdfe re_prod lag_prod if Total  >=  5, absorb(mastercode )
estadd local captain "Yes"

eststo m34: reghdfe re_prod lag_prod if Total == 3, absorb(boatcode mastercode yearout)
estadd local ship "Yes"
estadd local captain "Yes"
eststo m35: reghdfe re_prod lag_prod if Total == 4, absorb(boatcode mastercode yearout )
estadd local ship "Yes"
estadd local captain "Yes"
eststo m36: reghdfe re_prod lag_prod if Total  >=  5, absorb(boatcode mastercode yearout )
estadd local ship "Yes"
estadd local captain "Yes"

esttab  m31 m32 m33 m34  m35 m36 using "Writeup/Tables/IE7do_tab4.tex", replace ///
    title("Effect of Experience on Production")    /// 
	mtitles("3 voy" "4 voy"  "5+  voy" "3 voy" "4 voy"  "5+  voy") ///
    label keep(lag_prod) b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    s(N captain ship time, label("N"  "Captain FE" "Ship FE" "Time FE"))  /// 
	nogaps noomit  addnotes(" ")