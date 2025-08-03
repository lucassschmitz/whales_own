clear all
set seed 123

// Parameters
local N = 30000  // individuals
local T = 11   // periods

// Create panel structure
set obs `N'
gen id = _n
expand `T'
bysort id: gen time = _n

// Generate individual fixed effects (mu)
bysort id: gen masterFE = rnormal(0,1) if _n==1
bysort id: replace masterFE = masterFE[1]

// Generate dependent variable y
gen error = rnormal(0,1)
gen y = masterFE + error + .1 * time 
gen y2 = masterFE + error 

xtset id time
bysort id (time): gen y_mean_prev = sum(L.y) / (time- 1) if time > 1

reg y L.y // positive lag coefficient
reg y L.y time
reg y2 L.y2 time

reghdfe y L.y, absorb(id) // negative
reghdfe y L.y time, absorb(id)


reghdfe y2 L.y2, absorb(id) // lag < 0 
reghdfe y2 L.y2 time, absorb(id) // lag < 0 time =0 

// wrong coefficients
xtabond2 y L.y, gmm(L.y, lag(2 .) collapse)  twostep robust   // coeff == 1 
xtabond2 y L.y , gmm(L.y, eq(d)) iv() twostep robust		// coeff < 0 

// correct coefficients once controls for experience. 
xtabond2 y L.y time, gmm(L.y, lag(2 .) collapse) iv(time) twostep robust // lag = 0, time = .1 
xtabond2 y L.y time, gmm(L.y, eq(d)) iv(time) twostep robust // lag  = 0, time = .1  

// if experience has no effect then it is not necessary to control for it. but does not change coefficient if included. 
xtabond2 y2 L.y2, gmm(L.y, lag(2 .) collapse)  twostep robust   // coeff == 0 
xtabond2 y2 L.y2 , gmm(L.y2, eq(d)) iv() twostep robust		// coeff = 0 

xtabond2 y2 L.y2 time, gmm(L.y, lag(2 .) collapse)  iv(time) twostep robust   // coeff == 0 
xtabond2 y2 L.y2 time , gmm(L.y2, eq(d)) iv(time) twostep robust		// coeff = 0 