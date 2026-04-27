*   Change your working directory to the example folder, then run this script:
*   cd /path/to/survAH_stata/example
*   do example-stah.do
 
discard

* Load PBC data
clear all
version 13
set more off
capture log close

local base_dir = ".."

log using "`base_dir'/example/example-stah.log", replace text

* Add the local ado folder
adopath ++ "`base_dir'/ado"

* Import the dataset
import delimited "`base_dir'/example/pbc_data.csv", clear

* rescale time from days to years
gen years = time/365.25

* Create event variable (death = status 2)
gen byte event = (status == 2)

* Convert string treatment to numeric
gen trt_num = .
replace trt_num = 1 if trt == "1"
replace trt_num = 2 if trt == "2"

* Create bilirubin stratification variable
gen bili_strata = (bili >= 3.0) if !missing(bili)

* Keep only complete cases
drop if missing(trt_num) | missing(bili)

* Set survival data
stset years, failure(event)

* 1. Single-arm analysis
stah, tau(5)

* 2. Two-sample analysis (treatment comparison)
stah trt_num, tau(5) reference(2)

* 3. Stratified analysis (by bilirubin)
stah trt_num, strata(bili_strata) tau(5) reference(2)

* 4. Stratified analysis with custom weights (emphasize high bilirubin, 1:2)
stah trt_num, strata(bili_strata) tau(5) reference(2) weights(1 2)

* Sensitivity analysis across tau values (stratified by bilirubin, default proportional weights)
di _newline "=================================================="
di "Sensitivity analysis: Stratified AH and Cox HR across tau = 1 to 7 years"
foreach tau in 1 2 3 4 5 6 7 {
    di _newline "--------------------------------------------------"
    di "tau = `tau' years"

    qui stset years, failure(event)

    * Stratified Average Hazard (RAH)
    stah trt_num, strata(bili_strata) tau(`tau') reference(2)

    * Cox HR (censored at tau)
    di _newline "Cox regression (censored at `tau' years):"
    qui stset years, failure(event) exit(time `tau')
    stcox ib2.trt_num
}

* Restore original stset
qui stset years, failure(event)

log close
