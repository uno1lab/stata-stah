/* stah.ado - Complete Average Hazard (AH) analysis
   Version 1.0

   Computes the Average Hazard with Survival Weight (AHSW), a summary measure that 
   quantifies the rate of event occurrence over a restricted time period [0, tau].
   Hereafter referred to as Average Hazard (AH).

   Usage:
   - Single-arm:   stah, tau(21)
   - Two-sample:   stah armvar, tau(21)
   - Stratified:   stah armvar, tau(21) strata(stratavar)
*/

program stah, rclass byable(recall)
    version 13
    // Must be stset
    if "`_dta[_dta]'"!="st" {
        di as err "data not st"
        exit 119
    }

    // detect if first token is an option (comma) -> single arm
    gettoken first_token : 0, parse(" ,")
    if substr("`first_token'",1,1)=="," {
        syntax [if] [in], [tau(numlist > 0 max=1 miss) level(cilevel)]
        marksample touse
        local completevars _t _d
        markout `touse' `completevars'
        quietly count if `touse'
        di " "
        di as input "Number of observations for analysis = " r(N)
        stah_single_arm, touse(`touse') tau(`tau') level(`level')
        return add
    }
    else {
        syntax varlist(max=1 numeric) [if] [in], ///
            [tau(numlist > 0 max=1 miss) level(cilevel) reference(numlist >= 0 max=1) strata(varlist max=1 numeric)]
        marksample touse
        local completevars _t _d `varlist' `strata'
        markout `touse' `completevars'
        quietly count if `touse'
        di " "
        di as input "Number of observations for analysis = " r(N)

        if "`strata'"!="" {
            stah_two_sample_stratified `varlist', touse(`touse') tau(`tau') level(`level') ///
                reference(`reference') strata(`strata')
            return add
        }
        else {
            stah_two_sample `varlist', touse(`touse') tau(`tau') level(`level') ///
                reference(`reference')
            return add
        }
    }
end

// --------------------------------------------------
// Single arm
// --------------------------------------------------
program stah_single_arm, rclass
    version 13
    syntax, touse(string) [tau(numlist > 0 max=1 miss) level(cilevel)]

    // globals used inside ah_calculate
    global level = cond("`level'"=="", 95, `level')
    global tau   = "`tau'"
    global time  "_t"
    global fail  "_d"

    // default tau if not provided
    if "$tau"=="" {
        quietly count if $fail==1 & `touse'
        if r(N)<10 {
            di as err "Insufficient events for analysis (need at least 10)"
            exit 198
        }
        preserve
            quietly keep if `touse'
            capture {
                quietly sts gen nrisk = n
                quietly count if nrisk>=10 & $fail==1
                if r(N)>0 {
                    quietly su $time if nrisk>=10 & $fail==1
                    local defaulttau = r(max)
                }
                else {
                    quietly su $time if $fail==1, detail
                    local defaulttau = r(p75)
                }
            }
        restore
        if `defaulttau'==. {
            quietly su $time if $fail==1 & `touse', detail
            local defaulttau = r(p75)
        }
        global tau = `defaulttau'
        di as input "The truncation time: tau = " %9.3f $tau " was specified (default)."
    }
    else {
        di as input "The truncation time: tau = $tau was specified."
    }

    // validate tau in used subset
    quietly su $time if `touse'
    if $tau>r(max) {
        di as err "tau ($tau) exceeds maximum observed time (" %9.3f r(max) ")"
        exit 198
    }

    // compute AH
    preserve
        quietly keep if `touse'
        ah_calculate
        local ah_est     = r(ah_est)
        local ah_se      = r(ah_se)
        local ah_var     = r(ah_var)
        local rmst_val   = r(rmst)
        local cuminc_val = r(cum_inc_tau)
        matrix ah_results = r(ah_output)
    restore

    // counts
    quietly count if `touse'
    local total_n = r(N)
    quietly count if $fail==1 & $time <= $tau & `touse'
    local events_by_tau = r(N)
    quietly count if $time > $tau & `touse'
    local at_risk_at_tau = r(N)
    quietly count if $fail==0 & $time <= $tau & `touse'
    local censored_by_tau = r(N)

    // display
    di _n as gr "Single-arm Average Hazard (AH)"
    di as gr "Truncation time (tau): " %7.3f $tau
    di _n as gr "Number of observations:"
    di as gr "{hline 12}{c TT}{hline 62}"
    di as gr "        " _col(13) "{c |}" _col(16) "Total N" _col(26) "Event by tau" _col(40) "Censor by tau" _col(54) "At risk at tau"
    di as gr "{hline 12}{c +}{hline 62}"
    di as gr %8s "Single arm" _col(12) " {c |}" as ye ///
        _col(18) %5.0f `total_n' _col(30) %5.0f `events_by_tau' _col(44) %5.0f `censored_by_tau' _col(58) %5.0f `at_risk_at_tau'
    di as gr "{hline 12}{c BT}{hline 62}"

    di _n as gr "Average Hazard:"
    di as gr "{hline 70}"
    di as gr "AH(tau)" _col(15) "{c |}" _col(20) "Estimate" _col(32) "Std. Err." _col(44) "[$level% Conf. Interval]"
    di as gr "{hline 14}{c +}{hline 55}"
    di as gr %12s "AH($tau)" _col(15) "{c |}" as ye ///
        _col(20) %8.3f el(ah_results,1,1) _col(32) %8.3f el(ah_results,1,2) ///
        _col(44) %8.3f el(ah_results,1,3) _col(55) %8.3f el(ah_results,1,4)
    di as gr "{hline 70}"

    // returns
    return matrix results   = ah_results
    return scalar ah        = `ah_est'
    return scalar ah_se     = `ah_se'
    return scalar ah_var    = `ah_var'
    return scalar rmst      = `rmst_val'
    return scalar cum_inc   = `cuminc_val'
    return scalar surv_tau  = 1-`cuminc_val'
    return scalar tau       = $tau
    return scalar n_total   = `total_n'
    return scalar n_events  = `events_by_tau'
    return scalar n_censored= `censored_by_tau'
    return scalar n_atrisk  = `at_risk_at_tau'
    return local  analysis_type "single_arm"
end

// --------------------------------------------------
// Two-sample (unstratified)
// --------------------------------------------------
program stah_two_sample, rclass
    version 13
    syntax varlist(max=1 numeric), touse(string) [tau(numlist > 0 max=1 miss) level(cilevel) reference(numlist >= 0 max=1)]

    global arm  "`varlist'"
    global level = cond("`level'"=="", 95, `level')
    global tau  = "`tau'"
    global time "_t"
    global fail "_d"
    global reference = "`reference'"

    quietly tab $arm if `touse', matrow(matname)
    if rowsof(matname)!=2 {
        di as err "error: exactly 2 treatment arms supported"
        exit 198
    }

    if "$reference"=="" {
        global reference = el(matname,1,1)
    }
    // find treatment (the other value)
    forvalues i=1/`=rowsof(matname)' {
        if el(matname,`i',1)!=$reference global treatment = el(matname,`i',1)
    }

    // reference
    preserve
        quietly keep if `touse' & $arm==$reference
        ah_calculate
        local ah_ref_est    = r(ah_est)
        local ah_ref_var    = r(ah_var)
        local ah_ref_logvar = r(ah_logvar)
        matrix ah_ref = r(ah_output)
    restore

    // treatment
    preserve
        quietly keep if `touse' & $arm==$treatment
        ah_calculate
        local ah_trt_est    = r(ah_est)
        local ah_trt_var    = r(ah_var)
        local ah_trt_logvar = r(ah_logvar)
        matrix ah_trt = r(ah_output)
    restore

    // contrasts
    local dah_est = `ah_trt_est' - `ah_ref_est'
    local dah_var = `ah_trt_var' + `ah_ref_var'
    local dah_se  = sqrt(`dah_var')
    local z = invnorm(1 - (100-$level)/2/100)
    local dah_low = `dah_est' - `z'*`dah_se'
    local dah_upp = `dah_est' + `z'*`dah_se'
    local dah_p   = cond(`dah_se'>0, 2*(1-normprob(abs(`dah_est'/`dah_se'))), .)

    local rah_est = .
    local rah_low = .
    local rah_upp = .
    local rah_p   = .
    if `ah_ref_est'>0 & `ah_trt_est'>0 {
        local l   = ln(`ah_trt_est') - ln(`ah_ref_est')
        local lv  = `ah_trt_logvar' + `ah_ref_logvar'
        if `lv'>0 {
            local lse = sqrt(`lv')
            local rah_est = exp(`l')
            local rah_low = exp(`l' - `z'*`lse')
            local rah_upp = exp(`l' + `z'*`lse')
            local rah_p   = 2*(1 - normprob(abs(`l'/`lse')))
        }
    }

    matrix unadj_results = (`dah_est', `dah_low', `dah_upp', `dah_p' \ `rah_est', `rah_low', `rah_upp', `rah_p')
    matrix rownames unadj_results = DAH RAH
    matrix colnames unadj_results = Estimate Lower$level% Upper$level% P

    // counts for display
    quietly count if `touse' & $arm==$reference
    local total_n_ref = r(N)
    quietly count if `touse' & $arm==$treatment
    local total_n_trt = r(N)
    quietly count if $fail==1 & $time <= $tau & `touse' & $arm==$reference
    local events_ref = r(N)
    quietly count if $fail==1 & $time <= $tau & `touse' & $arm==$treatment
    local events_trt = r(N)
    quietly count if $fail==0 & $time <= $tau & `touse' & $arm==$reference
    local censored_ref = r(N)
    quietly count if $fail==0 & $time <= $tau & `touse' & $arm==$treatment
    local censored_trt = r(N)
    quietly count if $time > $tau & `touse' & $arm==$reference
    local atrisk_ref = r(N)
    quietly count if $time > $tau & `touse' & $arm==$treatment
    local atrisk_trt = r(N)

    // display
    di _n as gr "Number of observations:"
    di as gr "{hline 12}{c TT}{hline 62}"
    di as gr "        " _col(13) "{c |}" _col(16) "Total N" _col(26) "Event by tau" _col(40) "Censor by tau" _col(54) "At risk at tau"
    di as gr "{hline 12}{c +}{hline 62}"
    di as gr %8s "arm$reference" _col(12) " {c |}" as ye ///
        _col(18) %5.0f `total_n_ref' _col(30) %5.0f `events_ref' _col(44) %5.0f `censored_ref' _col(58) %5.0f `atrisk_ref'
    di as gr %8s "arm$treatment" _col(12) " {c |}" as ye ///
        _col(18) %5.0f `total_n_trt' _col(30) %5.0f `events_trt' _col(44) %5.0f `censored_trt' _col(58) %5.0f `atrisk_trt'
    di as gr "{hline 12}{c BT}{hline 62}"

    di _n as gr "Average Hazard (AH) by arm:"
    di as gr "{hline 12}{c TT}{hline 40}"
    di as gr "        " _col(13) "{c |}" _col(16) "Est." _col(26) "Lower $level%" _col(38) "Upper $level%"
    di as gr "{hline 12}{c +}{hline 40}"
    di as gr %8s "AH (arm $reference)" _col(12) " {c |}" as ye ///
        _col(16) %8.3f el(ah_ref,1,1) _col(26) %8.3f el(ah_ref,1,3) _col(38) %8.3f el(ah_ref,1,4)
    di as gr %8s "AH (arm $treatment)" _col(12) " {c |}" as ye ///
        _col(16) %8.3f el(ah_trt,1,1) _col(26) %8.3f el(ah_trt,1,3) _col(38) %8.3f el(ah_trt,1,4)
    di as gr "{hline 12}{c BT}{hline 40}"

    di _n as gr "Between-group contrast:"
    di as gr "{hline 20}{c TT}{hline 44}"
    di as gr "     Contrast" _col(21) "{c |}" _col(24) "Est." _col(34) "Lower $level%" _col(46) "Upper $level%" _col(58) "P>|z|"
    di as gr "{hline 20}{c +}{hline 44}"
    di as gr %18s "DAH ($treatment - $reference)" _col(20) " {c |}" as ye ///
        _col(24) %8.3f el(unadj_results,1,1) _col(34) %8.3f el(unadj_results,1,2) _col(46) %8.3f el(unadj_results,1,3) _col(58) %8.3f el(unadj_results,1,4)
    di as gr %18s "RAH ($treatment / $reference)" _col(20) " {c |}" as ye ///
        _col(24) %8.3f el(unadj_results,2,1) _col(34) %8.3f el(unadj_results,2,2) _col(46) %8.3f el(unadj_results,2,3) _col(58) %8.3f el(unadj_results,2,4)
    di as gr "{hline 20}{c BT}{hline 44}"

    // returns
    return matrix ah_reference     = ah_ref
    return matrix ah_treatment     = ah_trt
    return matrix unadjusted_results = unadj_results
    return scalar tau        = $tau
    return scalar reference  = $reference
    return scalar treatment  = $treatment
    return local  analysis_type "two_sample"
end

// --------------------------------------------------
// Two-sample (stratified)
// --------------------------------------------------

program stah_two_sample_stratified, rclass
    version 13
    syntax varlist(max=1 numeric), touse(string) [tau(numlist > 0 max=1 miss) level(cilevel) ///
                                   reference(numlist >= 0 max=1)] strata(varlist max=1 numeric)

    global arm = "`varlist'"
    global level = `level'
    global tau = "`tau'"
    global time "_t"
    global fail "_d"
    global reference = "`reference'"
    global strata = "`strata'"

    quietly tab $arm if `touse', matrow(__arms)
    if rowsof(__arms) != 2 {
        di as error "Error: stratified analysis currently supports exactly 2 treatment arms"
        exit 198
    }

    if "$reference"=="" {
        global reference = el(__arms,1,1)
    }
    forvalues i=1/`=rowsof(__arms)' {
        if el(__arms,`i',1) != $reference {
            global treatment = el(__arms,`i',1)
        }
    }

    // unique strata values
    quietly tab $strata if `touse', matrow(__SV)
    local n_strata = rowsof(__SV)

    // storage matrices
    tempname F0 R0 F1 R1 N0 N1 NN
    matrix `F0' = J(`n_strata',1,.)
    matrix `R0' = J(`n_strata',1,.)
    matrix `F1' = J(`n_strata',1,.)
    matrix `R1' = J(`n_strata',1,.)
    matrix `N0' = J(`n_strata',1,0)
    matrix `N1' = J(`n_strata',1,0)
    matrix `NN' = J(`n_strata',1,0)

    matrix strata_ah_ref = J(`n_strata',4,.)
    matrix strata_ah_trt = J(`n_strata',4,.)

    // conventional IVW accumulators
    local pooled_dah_num      = 0
    local pooled_dah_denom    = 0
    local pooled_rah_lognum   = 0
    local pooled_rah_logdenom = 0
    local valid_strata        = 0

    // Get total sample sizes and display summary
    quietly count if `touse' & $arm==$reference
    local total_n0 = r(N)
    quietly count if `touse' & $arm==$treatment  
    local total_n1 = r(N)

    // Display overall sample counts by strata
    di _n as gr "Number of observations:"
    di as gr "{hline 12}{c TT}{hline 30}"
    di as gr "        " _col(13) "{c |}" _col(16) "total" _col(24) "arm$reference" _col(32) "arm$treatment"
    di as gr "{hline 12}{c +}{hline 30}"
    forvalues s = 1/`n_strata' {
        local g = el(__SV,`s',1)
        quietly count if `touse' & $strata==`g' & $arm==$reference
        local n0s = r(N)
        quietly count if `touse' & $strata==`g' & $arm==$treatment
        local n1s = r(N)
        local total_s = `n0s' + `n1s'
        di as gr %8s "strata`g'" _col(12) " {c |}" as ye ///
            _col(18) %3.0f `total_s' _col(26) %3.0f `n0s' _col(34) %3.0f `n1s'
    }
    di as gr %8s "total" _col(12) " {c |}" as ye ///
        _col(18) %3.0f `=`total_n0'+`total_n1'' _col(26) %3.0f `total_n0' _col(34) %3.0f `total_n1'
    di as gr "{hline 12}{c BT}{hline 30}"

    // Display event counts by arm
    quietly count if $fail==1 & $time <= $tau & `touse' & $arm==$reference
    local events_ref = r(N)
    quietly count if $fail==1 & $time <= $tau & `touse' & $arm==$treatment
    local events_trt = r(N)
    quietly count if $fail==0 & $time <= $tau & `touse' & $arm==$reference
    local censored_ref = r(N)
    quietly count if $fail==0 & $time <= $tau & `touse' & $arm==$treatment
    local censored_trt = r(N)
    quietly count if $time > $tau & `touse' & $arm==$reference
    local atrisk_ref = r(N)
    quietly count if $time > $tau & `touse' & $arm==$treatment
    local atrisk_trt = r(N)

    di _n as gr "        Total N Event by tau Censor by tau At risk at tau"
    di as gr "arm$reference" %8.0f `total_n0' %11.0f `events_ref' %12.0f `censored_ref' %13.0f `atrisk_ref'
    di as gr "arm$treatment" %8.0f `total_n1' %11.0f `events_trt' %12.0f `censored_trt' %13.0f `atrisk_trt'

    // ---------- PASS 1: per-stratum components ----------
    forvalues s = 1/`n_strata' {
        local g = el(__SV,`s',1)

        // reference arm in stratum g
        preserve
            quietly keep if `touse' & $arm==$reference & $strata==`g'
            quietly count
            matrix `N0'[`s',1] = r(N)
            if r(N)>0 {
                quietly ah_calculate
                matrix __tmp_ref = r(ah_output)
                matrix strata_ah_ref[`s',1] = el(__tmp_ref,1,1)
                matrix strata_ah_ref[`s',2] = el(__tmp_ref,1,2)
                matrix strata_ah_ref[`s',3] = el(__tmp_ref,1,3)
                matrix strata_ah_ref[`s',4] = el(__tmp_ref,1,4)
                matrix `F0'[`s',1] = r(cum_inc_tau)
                matrix `R0'[`s',1] = r(rmst)
                local ref_est    = r(ah_est)
                local ref_var    = r(ah_var)
                local ref_logvar = r(ah_logvar)
            }
            else {
                local ref_est    = .
                local ref_var    = .
                local ref_logvar = .
            }
        restore

        // treatment arm in stratum g
        preserve
            quietly keep if `touse' & $arm==$treatment & $strata==`g'
            quietly count
            matrix `N1'[`s',1] = r(N)
            if r(N)>0 {
                quietly ah_calculate
                matrix __tmp_trt = r(ah_output)
                matrix strata_ah_trt[`s',1] = el(__tmp_trt,1,1)
                matrix strata_ah_trt[`s',2] = el(__tmp_trt,1,2)
                matrix strata_ah_trt[`s',3] = el(__tmp_trt,1,3)
                matrix strata_ah_trt[`s',4] = el(__tmp_trt,1,4)
                matrix `F1'[`s',1] = r(cum_inc_tau)
                matrix `R1'[`s',1] = r(rmst)
                local trt_est    = r(ah_est)
                local trt_var    = r(ah_var)
                local trt_logvar = r(ah_logvar)
            }
            else {
                local trt_est    = .
                local trt_var    = .
                local trt_logvar = .
            }
        restore

        // conventional IVW updates
        if `ref_est'<. & `trt_est'<. & `ref_var'>0 & `trt_var'>0 {
            local d = `trt_est' - `ref_est'
            local v = `trt_var' + `ref_var'
            if `v'>0 {
                local w = 1/`v'
                local pooled_dah_num      = `pooled_dah_num'   + `w'*`d'
                local pooled_dah_denom    = `pooled_dah_denom' + `w'

                if `ref_logvar'>0 & `trt_logvar'>0 & `ref_est'>0 & `trt_est'>0 {
                    local l  = ln(`trt_est') - ln(`ref_est')
                    local vl = `trt_logvar' + `ref_logvar'
                    if `vl'>0 {
                        local wr = 1/`vl'
                        local pooled_rah_lognum   = `pooled_rah_lognum'   + `wr'*`l'
                        local pooled_rah_logdenom = `pooled_rah_logdenom' + `wr'
                    }
                }
                local valid_strata = `valid_strata' + 1
            }
        }
    }

    // Compute unadjusted AH estimates (pooling across strata, not stratified)
preserve
    quietly keep if `touse' & $arm==$reference
    quietly ah_calculate
    local unadj_ah0_est = r(ah_est)
    local unadj_ah0_se  = r(ah_se)
    local unadj_ah0_low = el(r(ah_output),1,3)
    local unadj_ah0_upp = el(r(ah_output),1,4)
    // log-based CI
    local unadj_ah0_logvar = r(ah_logvar)
    local unadj_ah0_low_log = exp(ln(`unadj_ah0_est') - 1.96*sqrt(`unadj_ah0_logvar'))
    local unadj_ah0_upp_log = exp(ln(`unadj_ah0_est') + 1.96*sqrt(`unadj_ah0_logvar'))
restore

    preserve
        quietly keep if `touse' & $arm==$treatment
        quietly ah_calculate
        local unadj_ah1_est = r(ah_est)
        local unadj_ah1_se  = r(ah_se)
        local unadj_ah1_low = el(r(ah_output),1,3)
        local unadj_ah1_upp = el(r(ah_output),1,4)
        // log-based CI
        local unadj_ah1_logvar = r(ah_logvar)
        local unadj_ah1_low_log = exp(ln(`unadj_ah1_est') - 1.96*sqrt(`unadj_ah1_logvar'))
        local unadj_ah1_upp_log = exp(ln(`unadj_ah1_est') + 1.96*sqrt(`unadj_ah1_logvar'))
    restore

    
    // sample-size weights for direct standardization
    scalar __Ntot = 0
    forvalues s = 1/`n_strata' {
        matrix `NN'[`s',1] = el(`N0',`s',1) + el(`N1',`s',1)
        scalar __Ntot = __Ntot + el(`NN',`s',1)
    }
    if __Ntot==0 {
        di as error "No observations after subsetting."
        exit 2000
    }

    // pooled F,R via direct standardization
    scalar __F0bar = 0
    scalar __R0bar = 0
    scalar __F1bar = 0
    scalar __R1bar = 0
    forvalues s = 1/`n_strata' {
        scalar __wt = el(`NN',`s',1) / __Ntot
        scalar __F0bar = __F0bar + __wt * el(`F0',`s',1)
        scalar __R0bar = __R0bar + __wt * el(`R0',`s',1)
        scalar __F1bar = __F1bar + __wt * el(`F1',`s',1)
        scalar __R1bar = __R1bar + __wt * el(`R1',`s',1)
    }
    scalar __eta0_bar = cond( __R0bar>0, __F0bar/__R0bar, .)
    scalar __eta1_bar = cond( __R1bar>0, __F1bar/__R1bar, .)
    scalar __dah_std  = __eta1_bar - __eta0_bar
    scalar __rah_std  = cond( __eta0_bar>0, __eta1_bar/__eta0_bar, .)
    
    // Initialize variance accumulators
    scalar __eta0_var     = 0
    scalar __eta1_var     = 0
    scalar __log_eta0_var = 0
    scalar __log_eta1_var = 0

    forvalues s = 1/`n_strata' {
        local g   = el(__SV,`s',1)
        scalar __wt = el(`NN',`s',1) / __Ntot

        local n0s = el(`N0',`s',1)
        local n1s = el(`N1',`s',1)

        // arm 0: use pooled F0bar/R0bar but the stratum's event grid & RMST path
        preserve
            quietly keep if `touse' & $arm==$reference & $strata==`g'
            quietly count
            if r(N)>0 & `n0s'>0 {
                quietly __ah_var_given_FRbar, fbar(`=__F0bar') rbar(`=__R0bar') nsub(`n0s')
                scalar __eta0_var     = __eta0_var     + r(var_ah)    * ((__wt^2) / `n0s')
                scalar __log_eta0_var = __log_eta0_var + r(var_logah) * ((__wt^2) / `n0s')
            }
        restore
        
        // arm 1
        preserve
            quietly keep if `touse' & $arm==$treatment & $strata==`g'
            quietly count
            if r(N)>0 & `n1s'>0 {
                quietly __ah_var_given_FRbar, fbar(`=__F1bar') rbar(`=__R1bar') nsub(`n1s')
                scalar __eta1_var     = __eta1_var     + r(var_ah)    * ((__wt^2) / `n1s')
                scalar __log_eta1_var = __log_eta1_var + r(var_logah) * ((__wt^2) / `n1s')
            }
        restore
    }


    local __z = invnorm(1 - (100-$level)/2/100)
    local __z_precise = 1.959963984540    
    local __z_use = `__z_precise'

    // DS contrasts with variance and precise z-value
    scalar __dah_se = sqrt(__eta0_var + __eta1_var)
    scalar __dah_L  = __dah_std - `__z_use'*__dah_se
    scalar __dah_U  = __dah_std + `__z_use'*__dah_se
    scalar __dah_P  = cond(__dah_se>0, 2*(1 - normprob(abs(__dah_std/__dah_se))), .)

    scalar __rah_L = .
    scalar __rah_U = .
    scalar __rah_P = .
    scalar __vlog = __log_eta0_var + __log_eta1_var
    if (__rah_std>0 & __vlog>0) {
        scalar __lhat = ln(__eta1_bar) - ln(__eta0_bar)
        scalar __lse  = sqrt(__vlog)
        scalar __rah_L = exp(__lhat - `__z_use'*__lse)
        scalar __rah_U = exp(__lhat + `__z_use'*__lse)
        scalar __rah_P = 2*(1 - normprob(abs(__lhat/__lse)))
    }

    matrix stratified_results_ds = ( __dah_std, __dah_L, __dah_U, __dah_P \ __rah_std, __rah_L, __rah_U, __rah_P )
    matrix rownames stratified_results_ds = DAH RAH
    matrix colnames stratified_results_ds = Estimate Lower$level% Upper$level% P

    // ---------- IVW contrasts ----------
    capture matrix drop stratified_results_ivw
    matrix stratified_results_ivw = J(2,4,.)
    matrix rownames stratified_results_ivw = DAH RAH
    matrix colnames stratified_results_ivw = Estimate Lower$level% Upper$level% P

    if `valid_strata' >= 1 & `pooled_dah_denom' > 0 {
        local zcrit = invnormal((100 + $level)/200)

        local Dhat = `pooled_dah_num' / `pooled_dah_denom'
        local Dse  = sqrt(1/`pooled_dah_denom')
        local Dlo  = `Dhat' - `zcrit' * `Dse'
        local Dup  = `Dhat' + `zcrit' * `Dse'
        local Dp   = cond(`Dse' > 0, 2*(1 - normal(abs(`Dhat'/`Dse'))), .)

        matrix stratified_results_ivw[1,1] = `Dhat'
        matrix stratified_results_ivw[1,2] = `Dlo'
        matrix stratified_results_ivw[1,3] = `Dup'
        matrix stratified_results_ivw[1,4] = `Dp'

        if `pooled_rah_logdenom' > 0 {
            local Lhat = `pooled_rah_lognum' / `pooled_rah_logdenom'
            local Lse  = sqrt(1/`pooled_rah_logdenom')

            local Rest = exp(`Lhat')
            local Rlo  = exp(`Lhat' - `zcrit' * `Lse')
            local Rup  = exp(`Lhat' + `zcrit' * `Lse')
            local Rp   = cond(`Lse' > 0, 2*(1 - normal(abs(`Lhat'/`Lse'))), .)

            matrix stratified_results_ivw[2,1] = `Rest'
            matrix stratified_results_ivw[2,2] = `Rlo'
            matrix stratified_results_ivw[2,3] = `Rup'
            matrix stratified_results_ivw[2,4] = `Rp'
        }
    }
    // Display standardized (adjusted) AH by arm 
    scalar __se0 = sqrt(__eta0_var)
    scalar __se1 = sqrt(__eta1_var)

    scalar __L0 = .
    scalar __U0 = .
    scalar __L1 = .
    scalar __U1 = .

    // log-scale CI if available; otherwise Wald on AH scale
    if (__eta0_bar>0 & __log_eta0_var>0) {
        scalar __L0 = exp( ln(__eta0_bar) - `__z_use'*sqrt(__log_eta0_var) )
        scalar __U0 = exp( ln(__eta0_bar) + `__z_use'*sqrt(__log_eta0_var) )
    }
    else {
        scalar __L0 = __eta0_bar - `__z_use'*__se0
        scalar __U0 = __eta0_bar + `__z_use'*__se0
    }

    if (__eta1_bar>0 & __log_eta1_var>0) {
        scalar __L1 = exp( ln(__eta1_bar) - `__z_use'*sqrt(__log_eta1_var) )
        scalar __U1 = exp( ln(__eta1_bar) + `__z_use'*sqrt(__log_eta1_var) )
    }
    else {
        scalar __L1 = __eta1_bar - `__z_use'*__se1
        scalar __U1 = __eta1_bar + `__z_use'*__se1
    }

    // Put in a matrix for return & display
    matrix stratified_ahsw = ( __eta0_bar, __se0, __L0, __U0 \ __eta1_bar, __se1, __L1, __U1 )
    matrix rownames stratified_ahsw = AH_arm$reference AH_arm$treatment
    matrix colnames stratified_ahsw = Estimate StdErr Lower$level% Upper$level%

    di _n as gr "<Adjusted analysis> Average Hazard (AH) by arm:"
    di as gr "{hline 12}{c TT}{hline 52}"
    di as gr "        " _col(13) "{c |}" _col(16) "Est." _col(26) "Std. Err." _col(38) "Lower $level%" _col(52) "Upper $level%"
    di as gr "{hline 12}{c +}{hline 52}"
    di as gr %8s "AH (arm$reference)"  _col(12) " {c |}" as ye ///
        _col(16) %8.3f el(stratified_ahsw,1,1) _col(26) %8.3f el(stratified_ahsw,1,2) ///
        _col(38) %8.3f el(stratified_ahsw,1,3) _col(52) %8.3f el(stratified_ahsw,1,4)
    di as gr %8s "AH (arm$treatment)" _col(12) " {c |}" as ye ///
        _col(16) %8.3f el(stratified_ahsw,2,1) _col(26) %8.3f el(stratified_ahsw,2,2) ///
        _col(38) %8.3f el(stratified_ahsw,2,3) _col(52) %8.3f el(stratified_ahsw,2,4)
    di as gr "{hline 12}{c BT}{hline 52}"


    // Display results (only standardized version)
    di _n in gr "Between-group contrast:"
    di in smcl in gr "{hline 20}{c TT}{hline 44}"
    di in smcl in gr " Contrast" _col(21) "{c |}" ///
                     _col(24) "Est." _col(34) "Lower $level%" _col(46) "Upper $level%" _col(58) "P>|z|"
    di in smcl in gr "{hline 20}{c +}{hline 44}"
    di in smcl in gr %18s "DAH ($treatment - $reference)" _col(20) " {c |}" in ye ///
         _col(24) %8.3f el(stratified_results_ds,1,1) ///
         _col(34) %8.3f el(stratified_results_ds,1,2) ///
         _col(46) %8.3f el(stratified_results_ds,1,3) ///
         _col(58) %8.3f el(stratified_results_ds,1,4)
    di in smcl in gr %18s "RAH ($treatment / $reference)" _col(20) " {c |}" in ye ///
         _col(24) %8.3f el(stratified_results_ds,2,1) ///
         _col(34) %8.3f el(stratified_results_ds,2,2) ///
         _col(46) %8.3f el(stratified_results_ds,2,3) ///
         _col(58) %8.3f el(stratified_results_ds,2,4)
    di in smcl in gr "{hline 20}{c BT}{hline 44}"

    // returns
    return matrix stratified_results    = stratified_results_ivw
    return matrix stratified_ahsw = stratified_ahsw
    return matrix stratified_results_ds = stratified_results_ds
    return matrix strata_ah_reference   = strata_ah_ref
    return matrix strata_ah_treatment   = strata_ah_trt
    return scalar tau        = $tau
    return scalar reference  = $reference
    return scalar treatment  = $treatment
    return scalar n_strata   = `n_strata'
    return local analysis_type "two_sample_stratified"
    return local strata_var "$strata"
end

// --------------------------------------------------
// Variance helper function 
// --------------------------------------------------
// Variance helper using pooled Fbar/Rbar (R-style)
// var = sum_t [ dH(t) * weight / G(t) ], with dH=D/N, G=N/n_sub
// => sum_t [ (D/N^2) * weight * n_sub ]
// --------------------------------------------------
program __ah_var_given_FRbar, rclass
    version 13
    syntax , Fbar(real) Rbar(real) Nsub(real)

    preserve
        capture {
            quietly sts gen __S = s __N = n __D = d
        }
        if _rc {
            quietly sts list, saving(__tmp_sts2, replace)
            use __tmp_sts2, clear
            rename time $time
            rename survivor __S
            rename n_risk   __N
            rename n_event  __D
        }

        keep $time __S __N __D
        quietly duplicates drop
        keep if $time <= $tau
        sort $time

        capture drop __had_tau
        gen byte __had_tau = ($time == $tau)
        quietly count if __had_tau
        if r(N) == 0 {
            quietly count
            local newobs = r(N) + 1
            quietly set obs `newobs'
            replace $time = $tau in `newobs'
            replace __S   = __S[_n-1] if $time == $tau & _n > 1
            replace __S   = 1          if $time == $tau & _n == 1
            replace __N   = 0          if $time == $tau
            replace __D   = 0          if $time == $tau
        }
        drop __had_tau

        // RMST path for this stratum (left-continuous)
        gen double __dt    = $time - $time[_n-1] if _n > 1
        replace   __dt    = $time if _n == 1
        gen double __Sleft = 1
        replace   __Sleft = __S[_n-1] if _n > 1
        gen double __Rpath = sum(__dt * __Sleft)

        tempvar dV wU wX vU vX
        gen double `dV' = .
        replace `dV' = __D / (__N^2) if __D > 0 & __N > 0     // dH/G = (D/N^2) * n_sub; we'll multiply by n_sub at the end

        gen double `wU' = ( 1/`rbar' - ( `fbar' * __Rpath) / (`rbar'^2) )^2
        gen double `wX' = .
        replace `wX'    = ( 1/`fbar' - ( __Rpath / `rbar') )^2 if `fbar' > 0

        gen double `vU' = `dV' * `wU' if `dV' < . & `wU' < .
        gen double `vX' = `dV' * `wX' if `dV' < . & `wX' < .

        quietly summarize `vU', meanonly
        scalar __VU = r(sum)  * `nsub'     // ← multiply by n_sub = 1/G scaling
        quietly summarize `vX', meanonly
        scalar __VX = r(sum)  * `nsub'

        return scalar var_ah    = __VU
        return scalar var_logah = __VX
    restore
end


// --------------------------------------------------
// Core calculation function with variance and CIs
// --------------------------------------------------
program ah_calculate, rclass
    version 13
    syntax [, debug]

    quietly {
        count
        if r(N)==0 {
            di as error "No observations in subset"
            exit 2000
        }

        capture {
            quietly sts gen surv = s nrisk = n nevent = d
        }
        if _rc != 0 {
            preserve
                quietly sts list, saving(__tmp_sts, replace)
                use __tmp_sts, clear
                rename time     $time
                rename survivor surv
                rename n_risk   nrisk
                rename n_event  nevent
                keep $time surv nrisk nevent
                quietly duplicates drop
                keep if $time <= $tau
                sort $time
                capture drop __had_tau
                gen byte __had_tau = ($time == $tau)
                quietly count if __had_tau
                if r(N) == 0 {
                    quietly count
                    local newobs = r(N) + 1
                    quietly set obs `newobs'
                    replace $time = $tau in `newobs'
                }
                drop __had_tau
                recast double surv nrisk nevent, force

                scalar __S_tau = surv[_N]
                if missing(__S_tau) {
                    quietly count if surv < .
                    local lastnm = r(N)
                    scalar __S_tau = cond(`lastnm'>0, surv[`lastnm'], 1)
                }
                replace surv = __S_tau if $time == $tau
                scalar __F_tau = 1 - __S_tau

                gen double dt = .
                replace dt = $time in 1
                replace dt = $time - $time[_n-1] if _n > 1
                gen double S_left = 1
                replace S_left = surv[_n-1] if _n > 1
                gen double area = dt * S_left
                gen double rmst_path = sum(area)
                scalar __R_tau = rmst_path[_N]

                scalar __AH = cond( __R_tau>0, __F_tau/__R_tau, 0 )

                tempvar dV wU wX vU vX
                gen double `dV' = .
                replace `dV' = nevent / (nrisk^2) if nevent > 0 & nrisk > 0
                gen double `wU' = ( 1/__R_tau - (__F_tau * rmst_path)/(__R_tau^2) )^2
                gen double `wX' = .
                replace `wX' = ( 1/__F_tau - (rmst_path / __R_tau) )^2 if __F_tau > 0
                gen double `vU' = `dV' * `wU' if `dV' < . & `wU' < .
                gen double `vX' = `dV' * `wX' if `dV' < . & `wX' < .
                quietly summarize `vU', meanonly
                scalar __V_AH = r(sum)
                scalar __SE    = sqrt(max(0, __V_AH))
                quietly summarize `vX', meanonly
                scalar __V_LOG = r(sum)

                local z = invnorm(1 - (100-$level)/2/100)
                scalar __L = .
                scalar __U = .
                if (__AH>0 & __V_LOG>0) {
                    scalar __L = exp( ln(__AH) - `z' * sqrt(__V_LOG) )
                    scalar __U = exp( ln(__AH) + `z' * sqrt(__V_LOG) )
                }
                else {
                    scalar __L = __AH - `z' * __SE
                    scalar __U = __AH + `z' * __SE
                }

                matrix __M = (__AH, __SE, __L, __U)
                matrix rownames __M = AH
                matrix colnames __M = Estimate StdErr Lower$level% Upper$level%

                return matrix ah_output   = __M
                return scalar ah_est      = __AH
                return scalar ah_se       = __SE
                return scalar ah_var      = __V_AH
                return scalar ah_logvar   = __V_LOG
                return scalar rmst        = __R_tau
                return scalar cum_inc_tau = __F_tau
            restore
            exit
        }

        // If sts gen worked (no fallback file roundtrip):
        keep $time surv nrisk nevent
        quietly duplicates drop
        keep if $time <= $tau
        sort $time
        capture drop __had_tau
        gen byte __had_tau = ($time == $tau)
        quietly count if __had_tau
        if r(N) == 0 {
            quietly count
            local newobs = r(N) + 1
            quietly set obs `newobs'
            replace $time = $tau in `newobs'
        }
        drop __had_tau
        recast double surv nrisk nevent, force

        scalar __S_tau = surv[_N]
        if missing(__S_tau) {
            quietly count if surv < .
            local lastnm = r(N)
            scalar __S_tau = cond(`lastnm'>0, surv[`lastnm'], 1)
        }
        replace surv = __S_tau if $time == $tau
        scalar __F_tau = 1 - __S_tau

        gen double dt = .
        replace dt = $time in 1
        replace dt = $time - $time[_n-1] if _n > 1
        gen double S_left = 1
        replace S_left = surv[_n-1] if _n > 1
        gen double area = dt * S_left
        gen double rmst_path = sum(area)
        scalar __R_tau = rmst_path[_N]

        scalar __AH = cond( __R_tau>0, __F_tau/__R_tau, 0 )

        tempvar dV wU wX vU vX
        gen double `dV' = .
        replace `dV' = nevent / (nrisk^2) if nevent > 0 & nrisk > 0
        gen double `wU' = ( 1/__R_tau - (__F_tau * rmst_path)/(__R_tau^2) )^2
        gen double `wX' = .
        replace `wX' = ( 1/__F_tau - (rmst_path / __R_tau) )^2 if __F_tau > 0
        gen double `vU' = `dV' * `wU' if `dV' < . & `wU' < .
        gen double `vX' = `dV' * `wX' if `dV' < . & `wX' < .
        quietly summarize `vU', meanonly
        scalar __V_AH = r(sum)
        scalar __SE    = sqrt(max(0, __V_AH))
        quietly summarize `vX', meanonly
        scalar __V_LOG = r(sum)

        local z = invnorm(1 - (100-$level)/2/100)
        scalar __L = .
        scalar __U = .
        if (__AH>0 & __V_LOG>0) {
            scalar __L = exp( ln(__AH) - `z' * sqrt(__V_LOG) )
            scalar __U = exp( ln(__AH) + `z' * sqrt(__V_LOG) )
        }
        else {
            scalar __L = __AH - `z' * __SE
            scalar __U = __AH + `z' * __SE
        }

        matrix __M = (__AH, __SE, __L, __U)
        matrix rownames __M = AH
        matrix colnames __M = Estimate StdErr Lower$level% Upper$level%

        return matrix ah_output   = __M
        return scalar ah_est      = __AH
        return scalar ah_se       = __SE
        return scalar ah_var      = __V_AH
        return scalar ah_logvar   = __V_LOG
        return scalar rmst        = __R_tau
        return scalar cum_inc_tau = __F_tau
    }
end