clear all
cd ..

set seed 123

// Parameters
local N = 30000  // individuals
local T = 14   // periods

// Create panel structure
set obs `N'
gen id = _n
expand `T'
bysort id: gen time = _n

scalar mu_c = 10 
scalar sd_c = 2
scalar sd_e = 6
scalar standard = 9 // captains with a lower standard are not hired.

// Generate individual fixed effects (mu)
bysort id: gen masterFE = rnormal(mu_c , sd_c ) if _n==1
bysort id: replace masterFE = masterFE[1]

// Generate dependent variable y
gen error = rnormal(0, sd_e)
gen y = masterFE + error 

// sum up to the point 
bysort id (time): gen y_cumsum = sum(y)

gen s2_post = ((1/sd_c^2) +  (time/sd_e^2) ) ^(-1) 
gen mu_post = s2_post *((mu_c/sd_c^2)+ (y_cumsum/sd_e^2)) 

gen fired = (standard >= mu_post) 
* gen hired = (standard < mu_post) 

bysort id (time): gen not_hired  =  (sum(fired) > fired)
gen observed = 1 - not_hired

bys id: egen Total = total(observed)

summarize Total if time == 1

bysort id (time): gen master_last = (observed == 1) & ///
                              (observed[_n+1] == 0 | _n == _N)

//////////// Initial regressions  analysis /////////////
xtset id time 

eststo m1: reg master_last y // < 0 
eststo m2: reg master_last L.y // < 0 but smaller, requries big sample to get power (more than 20k) 
*eststo m3: reg master_last L.y if Total < T-1  // still negative. 
eststo m4: reg master_last y L.y  // b_y < b_lag(y) < 0 

// the coefficient is an order of magnitude bigger for captains who leave early. 
eststo m5: reg master_last y if Total < 6
eststo m6: reg master_last y if Total > 6 

// the coefficient is an order of magnitude bigger for bad captaisn 
eststo m7: reg master_last y if masterFE < 9
eststo m8: reg master_last y if masterFE > 11 

esttab  m1 m2 m4 m5 m6 m7 m8  using "Writeup/Tables/IE8do_simulations.tex", replace ///
    title("Determinants of captain exit")    /// 
	mtitles("All" "All" "All" "Total<6" "Total>6" "Low type" "High type"  )   nocon /// keep(lag_prod)
    label  b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    s(r2 N )  /// 
	nogaps noomit  addnotes(" ")
	
	



reg Total masterFE if time == 1 // > 0 
binscatter Total masterFE if time==1, nquantiles(20)
graph export "Figures/meeting_250610/IE8_sim_Totalvoyages_quality.png", replace

est clear 
save Data/temp3, replace

//////// Trajectory beliefs over time /////// 

use Data/temp3, clear

preserve
egen FE_group = cut(masterFE), group(3) // low, medium, high skill groups
collapse (mean) mu_post, by(time FE_group)
xtline mu_post, overlay i(FE_group) t(time) xtitle("Voyage number")  ///
	legend(label(1 "Low Skill") label(2 "Medium Skill") label(3 "High Skill")) ///
	ytitle("Posterior") note("Captains grouped into skill levels the first period. The plot includes production not realized because captains were fired.")
graph export "Figures/meeting_250610/IE8_sim_posterior_time.png", replace
restore

preserve
egen FE_group = cut(masterFE), group(3) // low, medium, high skill groups
keep if observed == 1
collapse (mean) mu_post, by(time FE_group)
xtline mu_post, overlay i(FE_group) t(time) xtitle("Voyage number")  ///
	legend(label(1 "Low Skill") label(2 "Medium Skill") label(3 "High Skill")) ///
	ytitle("Posterior") note("Captains grouped into skill levels the first period. The plot includes production realized not including fired captains.")
graph export "Figures/meeting_250610/IE8_sim_posterior_time(2).png", replace
restore

/////////////// Cumulative firing rates. 
use Data/temp3, clear

xtile FE_quartile = masterFE if time==1, nq(4)

// Extend quartile groups to all periods for each id
bys id (time): replace FE_quartile = FE_quartile[1]

// 2. Generate cumulative firing indicator by time for each captain
bysort id (time): gen ever_fired = sum(fired) > 0

// 3. Compute cumulative fired rates over time by quartile
bysort FE_quartile time: egen cum_fire_rate = mean(ever_fired)

// 4. Plot cumulative firing rates over time, separately by quartiles
twoway ///
    (line cum_fire_rate time if FE_quartile==1, sort) ///
    (line cum_fire_rate time if FE_quartile==2, sort) ///
    (line cum_fire_rate time if FE_quartile==3, sort) ///
    (line cum_fire_rate time if FE_quartile==4, sort), ///
    xlabel(1(1)11) ylabel(0(.1)1) ///
    legend(label(1 "Q1 (Lowest skill)") ///
           label(2 "Q2") ///
           label(3 "Q3") ///
           label(4 "Q4 (Highest skill)")) ///
    title("Cumulative firing rate by initial skill quartile")
graph export "Figures/meeting_250610/IE8_sim_cumulativefiringrate_time.png", replace

///////// Belief upadting by CatpainFE quartiles. ///// 
use Data/temp3, clear

bysort id (time): gen update = mu_post - mu_post[_n-1] if observed == 1 

// Define initial skill quartiles
xtile FE_quartile = masterFE if time==1, nq(4)
bys id (time): replace FE_quartile = FE_quartile[1]

// Compute average update per time and quartile
bysort FE_quartile time: egen avg_update = mean(update)

// Plot average updates over time by quartile
twoway ///
    (line avg_update time if FE_quartile==1, sort) ///
    (line avg_update time if FE_quartile==2, sort) ///
    (line avg_update time if FE_quartile==3, sort) ///
    (line avg_update time if FE_quartile==4, sort), ///
    xlabel(1(1)11) ///
    legend(label(1 "Q1 (Lowest skill)") ///
           label(2 "Q2") ///
           label(3 "Q3") ///
           label(4 "Q4 (Highest skill)")) ///
    title("Average Update in Beliefs by Initial Skill Quartile") ///
    ytitle("Average Update in Posterior") xtitle("Time")

graph export "Figures/meeting_250610/IE8_sim_update_time.png", replace

///////
use Data/temp3, clear
keep if observed == 1

xtile y_quint = y  , nq(5)
collapse (mean) prob_last = master_last, by(y_quint)


twoway bar prob_last y_quint, ///
    xlabel(1 "Q1 (lowest L.y)" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5 (highest L.y)") ///
    ytitle("P(master_last = 1)") xtitle("L.y Quintile") ///
    note("Only considers observed voyages.") ///
    title("Probability of Last Voyage by L.y Quintile")
graph export "Figures/meeting_250610/IE8_sim_probexit_production.png", replace

  
//////// 
	use Data/temp3, clear
	// Create initial skill groups (quartiles)
xtile FE_quartile = masterFE if time==1, nq(4)
bys id (time): replace FE_quartile = FE_quartile[1]

// Compute probability of master_last per quartile and time
bysort FE_quartile time: egen prob_last = mean(master_last)

keep if observed == 1 
// Plot probability of last voyage over time by skill quartile
twoway (line prob_last time if FE_quartile==1 & time < T -1 , sort) ///
       (line prob_last time if FE_quartile==2 & time < T -1, sort) ///
       (line prob_last time if FE_quartile==3 & time < T -1, sort) ///
       (line prob_last time if FE_quartile==4 & time < T -1, sort), ///
       legend(label(1 "Q1 (Lowest skill)") ///
              label(2 "Q2") ///
              label(3 "Q3") ///
              label(4 "Q4 (Highest skill)")) ///
       xlabel(1(1)11)  /// ///
       title("Probability of Last Voyage by Skill Quartile") ///
       xtitle("Time") ytitle("Probability of Last Voyage")
graph export "Figures/meeting_250610/IE8_sim_exit_time_production.png", replace

////////// ////////////////////////////////////////

use Data/temp3, clear

statsby coef=_b[y], by(time) nodots clear: regress master_last y

// (3) You can list or plot these coefficients:
list time coef

twoway line coef time if time < 13 , sort ///
    ylabel(, angle(horizontal)) ///
    xtitle("Time") ///
    ytitle("β_{time} (coef on y)") ///
    title("Slope of master_last on y by Period")
 
 graph export "Figures/meeting_250610/IE8_sim_coeff_time.png", replace

 
 
 
use Data/temp3, clear

by id (time): gen lag_y     = y[_n-1]
gen      y_lag_int = y * lag_y

preserve
statsby coef=_b[y_lag_int] se  =_se[y_lag_int], by(time) nodots clear: ///
    regress master_last y lag_y y_lag_int

gen lb = coef - invnormal(0.975)*se
gen ub = coef + invnormal(0.975)*se

// 4. plot for time < 13
twoway (rcap ub lb time if time < 13, sort) /// 
    (line coef time if time < 13, sort lcolor(blue) lwidth(medium)), ///
    ylabel(, angle(horizontal)) ///
    xtitle("Time") ///
    ytitle("Coef on y×lag_y") ///
    title("β̂ and 95% CI for master_last on y×lag_y by Time")
graph export "Figures/meeting_250624/IE8sim_interaction_controls.png", replace

restore
	

preserve
statsby coef=_b[y_lag_int] se  =_se[y_lag_int], by(time) nodots clear: ///
    regress master_last  y_lag_int

gen lb = coef - invnormal(0.975)*se
gen ub = coef + invnormal(0.975)*se

// 4. plot for time < 13
twoway (rcap ub lb time if time < 13, sort) /// 
    (line coef time if time < 13, sort lcolor(blue) lwidth(medium)), ///
    ylabel(, angle(horizontal)) ///
    xtitle("Time") ///
    ytitle("Coef on y×lag_y") ///
    title("β̂ and 95% CI for master_last on y×lag_y by Time")
graph export "Figures/meeting_250624/IE8sim_interaction_no_controls.png", replace
restore	
	
////////////////////////////////////////////////////////////
 

 use Data/temp3, clear

gen prev_y = y_cumsum - y 

gen avg_y_prev = prev_y / (time - 1) 
replace avg_y_prev = .  if time == 1 	

egen FE_quartile = xtile(avg_y_prev) if observed == 1, by(time) nq(4) 

statsby coef = _b[y], by(FE_quartile time) nodots clear: regress master_last y
keep if time < 13 

twoway ///
    (line coef time if FE_quartile == 1 , sort lcolor(blue)    lpattern(solid)) ///
    (line coef time if FE_quartile == 2 , sort lcolor(red)     lpattern(dash))  ///
    (line coef time if FE_quartile == 3 , sort lcolor(green)   lpattern(dot))   ///
    (line coef time if FE_quartile == 4, sort lcolor(orange)  lpattern(longdash)), ///
    legend(order(1 "Q1 (Lowest skill)" ///
                 2 "Q2"                  ///
                 3 "Q3"                  ///
                 4 "Q4 (Highest skill)")) ///
    ylabel(, angle(horizontal)) ///
    xtitle("Time") ///
    ytitle("β̂_{t,q} (coef on y)") ///
    title("Coefficient on y by Time and shown quality Quartile")
graph export "Figures/meeting_250610/IE8_sim_coeffbyquartile_time.png", replace
