cd ..
use Data/masters_voyages_merged, clear

keep yearout product duration  tonn_avg mastercode boatcode product_nom
generate year10 = int(yearout / 10) * 10
generate year5 = int(yearout / 5) * 5

replace duration = 1 if duration == 0 
gen re_prod = product/duration 
gen re_prod_nom = product_nom/duration
gen log_p = log(product) 
gen log_rp = log(re_prod)

bys mastercode (yearout): gen n_voy = _n

summ product, d
summ re_prod, d 


bysort mastercode (yearout): gen lag_prod = re_prod[_n-1]

gen n_voy_top = n_voy 
replace n_voy_top = 5 if n_voy_top > 5

//////////////////// 
eststo m1: reghdfe re_prod n_voy,                absorb(mastercode boatcode yearout)
eststo m2: reghdfe re_prod n_voy lag_prod
eststo m3: reghdfe re_prod n_voy lag_prod,       absorb(mastercode)
eststo m4: reghdfe re_prod n_voy lag_prod,       absorb(boatcode)
eststo m5: reghdfe re_prod n_voy lag_prod,       absorb(mastercode boatcode yearout)

esttab m1 m2 m3 m4 m5 using "Writeup/Tables/do_reg.tex", replace ///
    title("Table 1: Effect of Experience on Production")        ///
    mtitles("(1)" "(2)" "(3)" "(4)" "(5)")                   ///
    label                                                      ///
    keep(n_voy lag_prod)                                      ///
    b(3) se(3)                                                 ///
    star(* 0.10 ** 0.05 *** 0.01)                              ///
    ///
    /// suppress default stats, add your own if desired
    nolabel nogaps noomit                                       ///
    ///
    /// add a note about SEs
    addnotes("Standard errors in parentheses")  


/// ANOVA test /// 

reghdfe re_prod n_voy lag_prod, absorb(mastercode yearout)
scalar sse_r = e(rss)
scalar df_r  = e(df_r)
	
reghdfe re_prod n_voy lag_prod, absorb(mastercode boatcode yearout)
scalar sse_u = e(rss)
scalar df_u  = e(df_r)
	
scalar df1 = df_r - df_u
scalar df2 = df_u
scalar F   = ((sse_r - sse_u)/df1) / (sse_u/df2)
scalar p   = Ftail(df1, df2, F)

display "F(" df1 "," df2 ") = " %6.2f F ",  p = " %6.3f p

////////////////// 
drop if n_voy > 5

 sort mastercode n_voy_top
twoway  (lpolyci  tonn_avg n_voy_top ,    lwidth(medium)), legend(off) ///
  xtitle("Experience") ytitle("Weight(tons)")  title("Weight Tonnage vs Experience") 
graph export "Figures/meeting_250527/IE6_2_weight_experience.png", replace

reghdfe tonn_avg , absorb(mastercode) residuals(resid_tonn)
  
twoway  (lpolyci  resid_tonn n_voy_top ,    lwidth(medium)), ///
  xtitle("Experience") ytitle("Residualized Weight") ///
 title("") note("Weight is residualized by captain FE") legend(off)
graph export "Figures/meeting_250527/IE6_2_residualizedweight_experience.png", replace
  
 
reghdfe tonn_avg , absorb(mastercode yearout) residuals(resid_tonn2)
  
twoway  (lpolyci  resid_tonn2 n_voy_top ,    lwidth(medium)), ///
  xtitle("Experience") ytitle("Residualized Weight") ///
 title("") note("Weight is residualized by captain FE") legend(off)
graph export "Figures/meeting_250527/IE6_2_residualizedweight_experience2.png", replace
  

  
  