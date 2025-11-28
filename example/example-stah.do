* Load PBC data
clear all
import delimited "pbc_data.csv", clear

* rescale time from day to year
gen year = time/365.25

* Create event variable (use "2" only)
gen event = (status == 2)

* Convert string treatment to numeric
gen trt_num = .
replace trt_num = 1 if trt == "1"
replace trt_num = 2 if trt == "2"

* Create bilirubin stratification
gen bili_strata = (bili >= 3.0) if !missing(bili)

* Keep only complete cases
drop if missing(trt_num) | missing(bili)

// Set survival data
stset year, failure(event)
sts graph, by(trt_num)

// AH
* 1. Single-arm analysis
stah, tau(5)

* 2. Two-sample analysis (treatment comparison)
stah trt_num, tau(5) reference(2)

* 3. Stratified analysis (by bilirubin)
stah trt_num, strata(bili_strata) tau(5) reference(2)


