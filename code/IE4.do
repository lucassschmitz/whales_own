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

summ product, d
summ re_prod, d 

histogram re_prod

sort re_prod tonn_avg
twoway (scatter re_prod tonn_avg) (lfit re_prod tonn_avg)
graph export "Figures/IE4_production_tonnage.png", replace

/// plots in levels 
	twoway (lpolyci re_prod yearout, yaxis(1) color(blue%20)) /// 
			(lpolyci product yearout, yaxis(2) color(red%20)) /// 		
			, ytitle("Product/Duration", axis(1)) ///
		   ytitle("Product", axis(2)) /// 
			xtitle("Year") ///	
		   legend(order(1 3) label(1 "Time adjusted") label(3 "Gross Revenue"))  
	graph export "Figures/IE4_prod&reprod_year.png", replace

	twoway (lpolyci re_prod_nom yearout, yaxis(1) color(blue%20)) /// 
			(lpolyci product_nom yearout, yaxis(2) color(red%20)) /// 		
			, ytitle("Product/Duration", axis(1)) ///
		   ytitle("Product", axis(2)) /// 
			xtitle("Year") title("Nominal Production") ///	
		   legend(order(1 3) label(1 "Time adjusted") label(3 "Gross Revenue"))  
	graph export "Figures/IE4_prod&reprod_year_nom.png", replace

	
	forvalues y = 1810(10)1905 {
		label define year10lbl `y' "`y'", add
	}
	label values year10 year10lbl


	regress product i.year10
	coefplot, drop(_cons)  	vertical    baselevels  xtitle("Year") /// 
		ytitle("Coefficient") title("Gross Production")
	graph export "Figures/IE4_prodFE.png", replace

	regress re_prod i.year10
	coefplot, drop(_cons) 	vertical    baselevels  xtitle("Year") ///
		ytitle("Coefficient") title("Production/Duration")
	graph export "Figures/IE4_reprodFE.png", replace


 ////// 
reghdfe re_prod i.year5, absorb(mastercode boatcode)
coefplot, drop(_cons) vertical  baselevels  xtitle("Year") /// 
	ytitle("Coefficient") title("Gross Production") 


reghdfe product i.year5, absorb(mastercode boatcode)
coefplot, drop(_cons)  	vertical    baselevels  xtitle("Year") /// 
	ytitle("Coefficient") title("Gross Production") 

