*! version 1.8  15April2026


capture program drop csestudy
program define csestudy, eclass
    version 14
    syntax varlist [if], EVENTstartdate(string) ///
        FIRSTPREeventdate(string) LASTPREeventdate(string) ///
        [gls npc(real 100) coefsonly WOODbury]

    _xt, trequired
    local panelvar = r(ivar) 
    local timevar = r(tvar)

    // Evaluate eventstartdate, firstpreeventdate, and lastpreeventdate
    local eventstartdate = `eventstartdate'
    local firstpreeventdate = `firstpreeventdate'
    local lastpreeventdate = `lastpreeventdate'

    local n_pre_event_days = `lastpreeventdate' - `firstpreeventdate' + 1

    // Check for valid event and pre-event dates
    capture assert `eventstartdate' > `lastpreeventdate'
    if _rc {
        di as error "Event start date must be after last pre-event date"
        exit 198
    }
    capture assert `lastpreeventdate' > `firstpreeventdate'
    if _rc {
        di as error "Last pre-event date must be after first pre-event date"
        exit 199
    }

    if !mi("`woodbury'") & mi("`gls'") {
        di as error "woodbury option requires gls"
        exit 198
    }

    if !mi("`gls'") {
        capture assert `npc' < = `n_pre_event_days'
        if _rc {
            di as error "Number of principal components must be less than or equal to the number of pre-event days"
            exit 200
        }
        if mi("`coefsonly'") {
            qui sum `timevar'
            capture assert r(min) <= `firstpreeventdate' - (`eventstartdate'-`firstpreeventdate')
            if _rc {
                local required_window = `eventstartdate'-`firstpreeventdate'
                di as error "Time variable must have `required_window' observations before the first pre-event date"
                exit 201
            }
        }
    }


    local cmdline "csestudy `0'"
    
    tokenize `varlist'
    local lhsvar `1'
    macro shift
    local rhsvars `*'
    scalar Dim = wordcount("`rhsvars' constant")

    // Mark all valid event and pre-event observations and create touse
    marksample touse_all

    // Set Data View
    tempvar data_window
    if mi("`gls'") {
        gen byte `data_window' = inrange(`timevar',`firstpreeventdate',`eventstartdate')
    }
    else {
        gen byte `data_window' = inrange(`timevar',`firstpreeventdate'-(`eventstartdate'-`firstpreeventdate'),`eventstartdate')
    }

    mata st_view(A = ., ., (               ///
        st_varindex("`panelvar'"),         ///
        st_varindex("`timevar'"),          ///
        st_varindex( "`touse_all'"),       ///
        st_varindex("`lhsvar'"),           ///
        st_varindex(tokens("`rhsvars'"))), ///
            st_varindex("`data_window'"))
    mata long_data = get_data_views(A)
    mata full_index = get_data_indexes(long_data,"`gls'")


    ****************************************************************************
    *                                 Run Test                                 *
    ****************************************************************************


    tempname b nobs

    //pre-allocate matrices
    local colnames `rhsvars' :_cons
    local ncols: word count `colnames'
    matrix `b' = J(1,`ncols',.)
    matrix rownames `b' = y1




    // Get event period coefficients
    mata current_index = get_current_indexes(full_index, ///
        `eventstartdate', `lastpreeventdate', `firstpreeventdate')

    // Check for number of valid observations in event-period (with or without balancing)
    mata st_local("valid_obs", strofreal(rows(current_index.touse_index)))
    if "`gls'" == "gls" {
        if `valid_obs' == 0 {
            di as error "Error: After balancing there are NO valid observations at event date `eventstartdate'."
            di as error "Event and pre-period windows have no panels with sequential observations."
            di as error "Check data structure to make sure trading date scheme is continguous."
            exit 202
        }

        if `valid_obs' < `npc' + 10 {
            di as error "Error: After balancing there are only `valid_obs' valid observations at event date `eventstartdate'."
            di as error "Event and pre-event panels lack a sufficient number"
            di as error "of sequential observations to perform PCA with `npc' components."
            di as error "Check for an unusually large number of gaps."
            exit 202
        }
    }
    else {
        if `valid_obs' == 0 {
            di as error "Error: There are no valid observations at date `eventstartdate'."
            di as error "Check data structure to make sure trading date scheme is continguous."
            exit 203
        }
    }


    mata _get_coefficients(long_data, current_index, "`b'", "`nobs'", "`gls'", `npc', "`woodbury'")


    // Label beta matrix
    local colnames `rhsvars' :_cons
    matrix rownames `b' = y1
    matrix colnames `b' = `colnames'
  

    if mi("`coefsonly'") {
        tempname all_betas all_nobs
        // Allocate all_betas matrix and store 
        // event period coefficient at start of matrix
        matrix `all_betas' = J(`n_pre_event_days'+1,`ncols',.)    
        
        local time_format: format `timevar'
        if substr("`time_format'",1,3) == "%tb" {
            local label_format `time_format'
        }

        local date_label: di `time_format' `eventstartdate'
        local all_betas_colnames `date_label'
        forval i = `lastpreeventdate' (-1) `firstpreeventdate' {
            local date_label: di `time_format' `i'
            local all_betas_colnames `all_betas_colnames' `date_label'
        }
        matrix colnames `all_betas' = `colnames'
        matrix rownames `all_betas' = `all_betas_colnames'

        matrix `all_betas'[1,1] = `b'

        // Store event period number of obs at start of all_nobs
        matrix `all_nobs' = J(`n_pre_event_days'+1,1,.)
        matrix colnames `all_nobs' = "N"
        matrix rownames `all_nobs' = `all_betas_colnames'

        matrix `all_nobs'[1,1] = `nobs'

        // Set reporting completion marker
        local percent_complete_last = 0

        if "`gls'" == "gls" {
            // GLS needs to test data back to n_pre_events prior
            local all_data_start_date = `firstpreeventdate' - (`eventstartdate'-`firstpreeventdate')
        }
        else {
            local all_data_start_date = `firstpreeventdate'            
        }

        // Loop through each pre-event date and run regression
        if "`gls'" == "gls" {
            di "Percent complete = 0% " _continue
        }

        tempname pre_event_b pre_event_nobs
        forval noevent_date = `lastpreeventdate'(-1)`firstpreeventdate' {

            if "`gls'" == "gls" {
                local noevent_lastpreeventdate = `noevent_date' - (`eventstartdate' - `lastpreeventdate')
                local noevent_firstpreeventdate = `noevent_date' - (`eventstartdate' - `firstpreeventdate')

            }
            else {
                // Under OLS, the pre-non-event window is unused
                local noevent_lastpreeventdate = .
                local noevent_firstpreeventdate = .
            }
            mata current_index = get_current_indexes(full_index, ///
                `noevent_date', `noevent_lastpreeventdate', `noevent_firstpreeventdate') 

            // Check for number of valid observations in event-period (with or without balancing)
            mata st_local("valid_obs", strofreal(rows(current_index.touse_index)))
            if "`gls'" == "gls" {
                if `valid_obs' == 0 {
                    di as error "Error: After balancing there are NO valid observations at pseudo-event date `eventstartdate'."
                    di as error "Psuedo-event and pre-period windows have no panels with sequential observations."
                    di as error "Check data structure to make sure trading date scheme is continguous."
                    exit 202
                }

                if `valid_obs' < `npc' + 10 {
                    di as error "Error: After balancing there are only `valid_obs' valid observations at event date `eventstartdate'."
                    di as error "Event and pre-event panels lack a sufficient number"
                    di as error "of sequential observations to perform PCA with `npc' components."
                    di as error "Check for an unusually large number of gaps."
                    exit 202
                }

            }
            else {
                if `valid_obs' == 0 {
                    di as error "Error: There are no valid observations at date `eventstartdate'."
                    di as error "Check data structure to make sure trading date scheme is continguous."
                    exit 203
                }
            }


        
            mata _get_coefficients(long_data, current_index, "`pre_event_b'", "`pre_event_nobs'", "`gls'", `npc', "`woodbury'")

            local j =  `lastpreeventdate' - `noevent_date' + 2
            matrix `all_betas'[`j',1] = `pre_event_b'
            matrix `all_nobs'[`j',1] = `pre_event_nobs'

            // Report percentage of gls regressions finished
            local pe_regno =  `lastpreeventdate' - `noevent_date' + 1
            local pe_total = `lastpreeventdate' - `firstpreeventdate' + 1
            local percent_complete = `pe_regno'/`pe_total' * 100
            if `percent_complete' - `percent_complete_last' >= 10 & "`gls'" == "gls" {
                di "... " %-2.0f `percent_complete' "% " _continue
                local percent_complete_last = `percent_complete'
            }            
        }
        
        tempname pcdf ts_z
        mata _get_significance_stats("`all_betas'", "`pcdf'", "`ts_z'")    
        matrix rownames `pcdf' = y1 
        matrix colnames `pcdf' = `rhsvars' :_cons

        matrix rownames `ts_z' = y1
        matrix colnames `ts_z' = `rhsvars' :_cons

        local event_start_date: display `label_format' `eventstartdate'
        local pre_event_start: display `label_format' `firstpreeventdate'
        local pre_event_end: display `label_format' `lastpreeventdate'


        di _n
        if !mi("`gls'") {
            if !mi("`woodbury'") {
                di as text "GLS Estimates with Time Series Corrected Errors (Woodbury)"
            }
            else {
                di as text "GLS Estimates with Time Series Corrected Errors"
            }
        }
        else {
            di as text "OLS Estimates with Time Series Corrected Errors"
        }

        di as text "{hline 61}"
        di as text "Event start date: "     _col(25) as result "`event_start_date'" 
        di as text "Pre-event window: "  _col(25) as result "`pre_event_start'"  ///
            as text " to " as result "`pre_event_end'"
        di as text "{hline 61}" _n
        
        di _col(36) as text "Number of obs  = " as result %9.0fc `nobs'
        di _col(24) as text "Number of pre-period dates = " as result %9.0fc `n_pre_event_days'

        di as text "{hline 13}{c TT}{hline 47}"
        di as text %12s abbrev("`lhsvar'",12)  " {c |}  Coefficient" _col(29) %~12s  "CDF p-val" _col(41)  %~12s  "Parametric p-val" 
        di as text "{hline 13}{c +}{hline 47}"
        foreach colnm in `rhsvars' _cons {
            di as text %12s abbrev("`colnm'",12) " {c |}"  ///
            _col(17) as result %9.0g `b'[1, colnumb(`b',"`colnm'") ]  ///
            _col(29) %9.3f `pcdf'[1, colnumb(`pcdf',"`colnm'") ]  ///
            _col(41) as result %9.3f `ts_z'[1, colnumb(`ts_z',"`colnm'") ] 
        }
        di as text "{hline 13}{c BT}{hline 47}" _n
    }

    ereturn post `b' , depname("`lhsvar'") esample(`touse')
    ereturn scalar N = `nobs'
    if mi("`coefsonly'") {
        ereturn matrix betas = `all_betas'
        ereturn matrix N_all_dates = `all_nobs'
        ereturn matrix p = `pcdf'
        ereturn matrix z = `ts_z'
        ereturn local event_start_date = "`event_output'"
        // Check whether there are an unusually small number of observations for some pre-event windows
        mata N_all_dates = st_matrix("e(N_all_dates)")
        mata st_local("event_nobs", strofreal(N_all_dates[1]))
        mata st_local("min_obs", strofreal(min(N_all_dates)))
        local min_frac = `min_obs'/ `event_nobs'
        if `min_frac' < 0.75 {
            di as error "Warning: Some pre-event windows have fewer than 75%"
            di as error "of the expected number of observations."
            di as error "Number of observations in event window = " _continue
            di as result "`event_nobs'"
            di as error "Minimum number of observations in a pre-event window = " _continue
            di as result "`min_obs'"
            di as error "Check matrix e(N_all_dates) for more details"
        }

    } 

end


findfile "csestudy.mata"
include "`r(fn)'"
