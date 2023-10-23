*! version 1.2  21Ocotber2023
* updated to new GLS routine - checks for missing preevent vars
* GLS uses return residuals
* correct CDF reporting 
* Assume sequential date use, either setup via bcal or manually

capture program drop csestudy
program define csestudy, eclass sortpreserve
    syntax varlist [if], EVENTdate(string) ///
    [ EVENTENDdate(string)   ///
    NPREeventdays(integer -99) STARTpreeventdate(string) ENDpreeventdate(string) ///
    noBALance gls npc(integer 100) PRESAMPLEmarker(name) noCHECKbalance timeit ///
    CUMRETurns(name) PRECALCulated newvar(name)]  


    _xt
    local panelvar = r(ivar) 
    local timevar = r(tvar)
    // Get format of timevar if bcal is enabled
    local datefmt : format `timevar'
    if substr("`datefmt'",1,3) != "%tb" {
        di _n "{hline 55}"
        display as text "Date is not formatted using bcal. While not necessary,"
        display as text "doing so is encouraged for ease of use and reporting."
        display as text "Non-sequential calendar date conventions which do not "
        display as text "account for weekends and holidays will produce errors."
        display as text "See help file for more details."
        di "{hline 55}" _n

    }


    // Validate nonmissing date inputs which allow string equations
    // like bofd("mycal",mdy(#,#,####))
    foreach dateinput in eventdate startpreeventdate endpreeventdate ///
        eventenddate {
        * Evaluate input macro
        if !mi("``dateinput''") {
            capture local `dateinput' = ``dateinput''
            capture confirm number ``dateinput''
            if _rc != 0 {
                di as error "Invalid value for `dateinput': ``dateinput''"
                exit _rc
            }
        }
    }

    // If eventenddate is missing, assume = eventdate (a single day event)
    if mi("`eventenddate'") {
        local eventenddate = `eventdate'
    }

    local event_window_length = `eventenddate' - `eventdate' + 1
    
    capture assert  `event_window_length' > 0
    if _rc != 0 {
        di as error "Event end date is before event start date"
        exit _rc
    }

    // If startpreeventdate is not specified, make it one event_window before eventdate
    if mi("`endpreeventdate'") {
        local endpreeventdate = `eventdate'  - `event_window_length'
    }

    // Make NPREeventdays = 200 the default if both it and ENDpreeventdate are missing
    // Note: npreeventdays defaults to -99 so it can be optional
    if `npreeventdays' == -99 & mi("`startpreeventdate'") {
        local npreeventdays = 200
        local startpreeventdate = `endpreeventdate' - `npreeventdays' + 1
    }
    // If npreeventdays is missing and endpreeventdate is specified,
    // calculate number of days
    else if `npreeventdays' == -99 & !mi("`startpreeventdate'") {
        local npreeventdays = `endpreeventdate' - `startpreeventdate' - 1
    }
    // Declare an error if both options are specified
    else if `npreeventdays' != -99 & !mi("`startpreeventdate'") {
        di as error "Cannot specify both -npreeventdays- and -startpreeventdate-"
        exit
    }
    else if `npreeventdays' < 1 {
        di as error "-npreeventdays- must be positive "
        exit
    } 
    else if `npreeventdays' > 1  {
        local startpreeventdate = `endpreeventdate' - `npreeventdays' + 1
    }
    else {
        di as error "Something has gone terribly wrong with this switch! FAIL!"
        exit
    }


    // Notify users who use the if option that it only applies
    // to event window observations
    if !mi("`if'") {
        di _n as text "NOTE: the -if- expression only applies to observations in the event window."
        di    as text "The pre-event window will exclude any panel ids excluded by -if- in the "
        di    as text "event window, but if you wish to exclude other pre-event window observations"
        di    as text "based on certain criteria, you must do so manually." _n
    }


    if !mi("`balance'") & !mi("`gls'") {
        di as error "GLS requires a balanced panel and does not allow nobalance."
        exit
    }

    local cmdline "csestudy `0'"

    time_section, clear `timeit'
    time_section, label(Setup) `timeit'

    capture confirm variable `presamplemarker'
    if _rc ==0 {
        di as error "Pre-sample marker variable name (`presamplemarker') already exists."
        exit
    }



    qui desc, short varlist
    local sortlist = r(sortlist)
    tokenize "`sortlist'"
    local sortvar1 `1'
    local sortvar2 `2'


    if `startpreeventdate'  >= `endpreeventdate' {
        di as error "Pre-event start date must come before end date."
        exit
    }

    if `eventdate'  <= `endpreeventdate' {
        di as error "Event Date come after pre-event end date."
        exit
    }

    di "{hline 50}"
    di as text "Event Window Length:    " as result %-4.0f `event_window_length'
    di as text "Event Window:           " as result `datefmt' `eventdate' " to " `datefmt'  `eventenddate'
    di as text "# of pre-event periods: " as result %-4.0f `npreeventdays'
    di as text "Pre-Event Window:       " as result `datefmt' `startpreeventdate' " to " `datefmt'  `endpreeventdate'
    di as text "{hline 50}"

    tokenize `varlist'
    local lhsvar `1'
    macro shift
    local rhsvars `*'
    scalar Dim = wordcount("`rhsvars' constant")


    time_section, label(Panel Setup) `timeit'

    * Create event variable to track event period
    tempvar event
    qui count if `timevar' == `eventdate'
    if r(N) == 0 {
        di as error "Event date does not occur in data. Check to make sure your format is correct."
        exit
    }
    qui gen byte `event' = `timevar' == `eventdate'

    marksample touse
    qui replace `touse' = 0 if !`event'
    
    local varlist_delim : subinstr local varlist " " ",", all 
    tempvar touse_pre_event
    qui gen byte `touse_pre_event' = inrange(`timevar',`startpreeventdate' , `endpreeventdate') ///
        & !mi(`varlist_delim')


    
    if "`sortvar1' `sortvar2'" != "`panelvar' `timevar'" {
        display "NOTE: Your data is not xt sorted. Sorting now..."
        sort  `panelvar' `timevar' 
    }

    // If event window length is greater than 1 and the precalculated option has not been specified
    // calculate the rolling window of returns
    local reset_lhsvar = 0
    if `event_window_length' > 1 & mi("`precalculated'") {
        time_section, label(Calc multi-day window) `timeit'
        local reset_lhsvar = 1

        if mi("`newvar'") {
            tempvar newvar
        }


        // Generate rolling lead window equation        
        local ret_eq `lhsvar'
        local j = `event_window_length' - 1
        forval i = 1/`j' {
            local ret_eq `ret_eq' + f`i'.`lhsvar'
        }
        
        qui gen `newvar' = `ret_eq' if inrange(`timevar',  `startpreeventdate' , `endpreeventdate')  | `timevar' == `eventdate'  

        // Preserve original lhsvar
        tempvar temp_lhsvar
        rename `lhsvar' `temp_lhsvar' 
        
        gen `lhsvar' = `newvar'        


    }


    // Try to balance the panel by excluding all panel ids which do not have
    // the maximum number of observed days
    if mi("`balance'") {
        time_section, label(Balancing Routine) `timeit'

        // Get maximum number of total trading dates for all stocks
        tempvar tcount
        qui by `panelvar' (`timevar'): gen long `tcount' /// 
            = sum(`touse'|`touse_pre_event') if `touse'|`touse_pre_event'
        // Copy non-empty count to the bottom of each panelvar
        qui by `panelvar' (`timevar'): replace `tcount' = `tcount'[_n - 1] /// 
            if mi(`tcount') ///
            | (!mi(`tcount'[_n - 1]) & `tcount'[_n - 1] > `tcount')        
        // Generate totaldays as last value of tcount, which has been copied to the bottom
        tempvar totaldays
        qui by `panelvar' (`timevar'): gen long `totaldays' = `tcount'[_N] if `touse'|`touse_pre_event'
        
        sum `totaldays', meanonly
        local maxdays = r(max)
        // Zero out the pre-events and post-events insample indicators where there are not the max obs
        qui replace `touse' = 0 if `totaldays' != `maxdays'
        qui replace `touse_pre_event' = 0 if `totaldays' != `maxdays'

        // Balance check if every panel-id exists for every date in touse_pre_event
        // This operation is somewhat expensive and it is very unlikely that the data 
        // passes the max-trading day operation but is somehow unbalanced, so it
        // can be supressed
        if mi("`checkbalance'") {
            time_section, label(Check Balance) `timeit'
            mata is_balanced("`panelvar'", "`timevar'","","`touse_pre_event'")
            if !r(all_id_unique) {
                di as smcl as err "Pre-period data cannot be balanced. Check for unusual gaps in data"
                exit
            }
        }
    }

    // Even if the data is not balanced, the command should still throw away all pre-period observations
    // that are screened out by -if- in the event period.
    else {
        tempvar tcount
        qui by `panelvar' (`timevar'): gen long `tcount' = `touse'
        // Copy non-empty count to the bottom of each panelvar
        qui by `panelvar' (`timevar'): replace `tcount' = `tcount'[_n - 1] if !mi(`tcount'[_n-1]) & `tcount'==0

        // Zero out touse_pre_event if last value of tcount is 0
        qui by `panelvar' (`timevar'): replace `touse_pre_event' = 0 if `tcount'[_N] == 0
    }

    qui count if `touse'
    local nobs = r(N)


    ****************************************************************************
    *                                 Run Test                                 *
    ****************************************************************************

    
    * Check that the number of pca components is less than the number of observations
    capture assert `nobs' > `npc'
    if _rc {
        di as error "After balancing pre-period, you have more pca components than test observations."
        exit
    }

    tempname b V ts_z pcdf ndates betas
    

    time_section, label(Mata Code) `timeit'

    mata _tsregress("`varlist'", "`timevar'", "`panelvar'" , "`gls'", "`balance'" , "`touse_pre_event'", "`touse'", `npc', "`b'" , "`V'", "`pcdf'", "`ts_z'", "`betas'")

    * local rhsvars w_ltdd w_me w_btm
    matrix rownames `b' = y1
    matrix colnames `b' = `rhsvars' :_cons

    matrix rownames `pcdf' = y1
    matrix colnames `pcdf' = `rhsvars' :_cons

    matrix rownames `ts_z' = y1
    matrix colnames `ts_z' = `rhsvars' :_cons

    matrix rownames `betas' = y1
    matrix colnames `betas' = `rhsvars' :_cons

    matrix rownames `V' = `rhsvars' :_cons
    matrix colnames `V' = `rhsvars' :_cons

    time_section, label(Post results to e) `timeit'
  
    if !mi("`gls'") {
        di as text "GLS Estimates with Time Series Corrected Errors"
    }
    else {
        di as text "OLS Estimates with Time Series Corrected Errors"
    }
    di _col(36) as text "Number of obs  = " as result %9.0fc `nobs'
    di _col(24) as text "Number of pre-period dates = " as result %9.0fc `ndates'

    di as text "{hline 13}{c TT}{hline 47}"
    di as text %12s abbrev("`lhsvar'",12)  " {c |}  Coefficient" _col(29) %~12s  "CDF p-val" _col(41)  %~12s  "TS Z-Score" 
    di as text "{hline 13}{c +}{hline 40}"
    foreach colnm in `rhsvars' _cons {
        di as text %12s abbrev("`colnm'",12) " {c |}"  ///
        _col(17) as result %9.0g `b'[1, colnumb(`b',"`colnm'") ]  ///
        _col(29) %9.3f `pcdf'[1, colnumb(`pcdf',"`colnm'") ]  ///
        _col(41) as result %9.3f `ts_z'[1, colnumb(`ts_z',"`colnm'") ] 
    }
    di as text "{hline 13}{c BT}{hline 40}" _n

    ereturn post `b' , depname("`lhsvar'") esample(`touse')
    ereturn local  cmd  "csestudy"
    ereturn local  cmdline  "`cmdline'"
    ereturn scalar N = `nobs'
    ereturn scalar n_pe_dates = `ndates'
    ereturn matrix z = `ts_z'
    ereturn matrix p = `pcdf'
    ereturn matrix betas = `betas'

    if !mi("`presamplemarker'") {
        qui gen byte `presamplemarker' = `touse_pre_event'
        la var `presamplemarker' "Pre-Event Sample Used"
    }

    if `reset_lhsvar' {
        drop `lhsvar'
        rename `temp_lhsvar' `lhsvar'
    }

    time_section, off `timeit'
    time_section, list `timeit'
end




capture mata mata drop _tsregress()
mata:
    void _tsregress(string scalar regvars, string scalar timevar, ///
        string scalar panelvar, string scalar gls, string scalar nobalance, ///
        string scalar pre_event_sample, string scalar event_sample, ///
        real scalar npc, ///
        string scalar bmat, string scalar Vmat, string scalar pcdf, ///
        string scalar ts_zmat, string scalar betas) {

        // Convert variable names to column indices
        real rowvector regcols
        real scalar time_column, startobs

        regcols          = st_varindex(tokens(regvars))
        time_column      = st_varindex(timevar)
        panel_column     = st_varindex(panelvar)
        event_sample_col = st_varindex(event_sample)
        pre_event_sample_col = st_varindex(pre_event_sample)


        // Preallocate regression objects
        real matrix residuals, XX, Xy, M, E, y, X, Xstack, ystack
        real scalar obs
        real vector beta, pre_event_dates, event_beta 
        // Setup Matrices for pre and post period
        st_view(E, . , (time_column,panel_column,regcols), event_sample_col)
        // M needs to be sorted, and partial sorting of the data may be required,
        // so M is data not a view
        M = st_data(. , (time_column,panel_column,regcols), pre_event_sample_col)
        M = sort(M,(1,2))

        pre_event_dates = uniqrows(M[.,1])
        firms = rows(uniqrows(M[.,2]))
        ndates = rows(pre_event_dates)
        pre_event_beta = J(ndates, 1 + cols(regcols), .)
        residuals = J(ndates, firms, .)
     
        // Allocate stacked X and y variables for fast GLS estimation
        if (gls == "gls") {
            Xstack = J(firms,ndates*cols(regcols),.)
            ystack = J(firms,ndates,.)
        }


        //calculate pre-event betas
        startobs = 1
        i = 1
        last_obs_flag = 0

        // Loop through to calculate OLS betas and, if GLS is specified, residuals.
        for (obs=1; obs <= rows(M); obs++) {

            if (obs == rows(M)) last_obs_flag = 1
            else if (M[obs,1] != M[obs+1,1]) last_obs_flag = 1
            
            if (last_obs_flag) {
                    y = M[startobs..obs, 3]
                    X = M[startobs..obs, 4..cols(M)] , J(obs - startobs + 1,1,1)
                    XX = quadcross(X,X)
                    Xy = quadcross(X,y)
                    // ------------ COMPUTE COEFFICIENTS --------------------
                    beta    = cholsolve(XX,Xy)
                    // ------------ STORE OUTPUT ----------------------------
                    curdate = M[obs,1]
                    pre_event_beta[i,.] = curdate, beta'
   
                    if (gls == "gls") {
                        // If GLS option is invoked, store residuals for PCA estimation
                        residuals[i,.] = (y-X*beta)'
                        // Stack X and y variables into matrices
                        // Note: this requires a balanced panel as per the GLS option
                        Xstack[.,1+cols(regcols)*(i-1)..cols(regcols)*i] = X
                        ystack[.,i] = y

                    }  

                    startobs = obs + 1 
                    last_obs_flag = 0
                    i = i + 1
            }
        } 


        if (gls == "gls") {    
            
            //Allocate PCA objects
            real matrix U, Vt, pca_coeff, pca_score
            real matrix sig2_e, Omega, L
            real vector s
            A = residuals :- mean(residuals)
            
            
            fullsvd(A,U,s,Vt)
            

            pca_coeff = Vt'[,1..npc]
            pca_score = A*pca_coeff
            sig2_e = variance(A - pca_score*pca_coeff')
            // Omega should be symmetric, but Mata doesn't recognize that it is
            // so we coerce it to be symmetric 
            Omega =  makesymmetric(pca_coeff * diag(variance(pca_score)) * // 
                pca_coeff' + diag(sig2_e))
            
            L = cholesky(Omega)
            

            // Solve for LX and Ly using cholesky decomposition
            Xstack = solvelower_wrapper(L,Xstack)
            ystack = solvelower_wrapper(L,ystack)
            

            // Loop though the stacked X and y tilde to get the time series GLS betas

            
            for (i = 1; i<=cols(ystack); i++) {
                X = Xstack[.,1+cols(regcols)*(i-1)..cols(regcols)*i]
                y = ystack[.,i]
                XX = quadcross(X,X)
                Xy = quadcross(X,y)
                beta = cholsolve(XX,Xy)
                pre_event_beta[i,2..cols(regcols)+1] = beta'
            }
            
        }

        st_subview(y, E, ., 3)
        st_subview(X, E, ., 4\.)
        X = X,J(rows(X),1,1)

        if (gls == "gls") {
            X = solvelower_wrapper(L,X)
            y = solvelower_wrapper(L,y)
        }
        XX = quadcross(X,X)
        Xy = quadcross(X,y)
        event_beta = E[1,1],cholsolve(XX,Xy)'


        // -----------------------------------------------------------------------------
        // Gather output and test
        // -----------------------------------------------------------------------------
        // Generate percentiles under empirical CDF
        coef_cols = 2..cols(pre_event_beta)
        mean_coefs = mean(pre_event_beta[.,coef_cols])
        event_pctile = colsum(abs(pre_event_beta[.,coef_cols]:-mean_coefs):> abs(event_beta[coef_cols]-mean_coefs)):/rows(pre_event_beta)
        
        // If percentile is 1 or 0, adjust by half the rank of N
        event_pctile = (event_pctile :== 0) * 1/(rows(pre_event_beta)*2) + ///
            (event_pctile :== 1) * (1- 1/(rows(pre_event_beta)*2)) + ///
            (event_pctile :!= 1 :& event_pctile :!= 0):* event_pctile

        // Generate TS empirical z-score
        X = (1\J(rows(pre_event_dates),1,0)),J(rows(pre_event_dates)+1,1,1)
        y = event_beta[.,coef_cols]\pre_event_beta[.,coef_cols]
        ts_beta = cholsolve(quadcross(X,X), quadcross(X,y))[1,.]
        sd = (sqrt(diagonal(quadvariance(pre_event_beta[.,coef_cols])) * ///
            (rows(pre_event_beta)+1) / rows(pre_event_beta)))'
        ts_z = ts_beta :/ sd
        VCV = diag(-abs(event_beta[coef_cols]):*invnormal(event_pctile/2):^-1):^2
        st_matrix(bmat, event_beta[.,coef_cols])
        st_matrix(Vmat, VCV)
        st_matrix(pcdf, event_pctile)
        st_matrix(ts_zmat, ts_z)
        st_matrix(betas, pre_event_beta[.,coef_cols])
        st_local("ndates", strofreal(ndates))
    }
end


// Replace solvelower function with solvelowerlapacke if 
// Stata is verison 17 or greater

if c(stata_version) >=17 {
    capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return(solvelowerlapacke(A,B))
}
else {
    capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return (solvelower(A,B))
}



capture mata mata drop is_balanced()
mata:
    void is_balanced( string scalar varlist, ///
        string scalar panelid, string scalar timeid, ///
        string scalar touse) {
        
        real matrix X, unique_ids, sample_marker
        real scalar valid_time, i

        st_view(X,.,(panelid, timeid, tokens(varlist)), touse)
        valid_time = 1
        for (i=2; i<= rows(unique_ids); i++) {
            if (unique_ids[1,2] != unique_ids[i,2]) {
                valid_time = 0
                break
            }
        }
        st_numscalar("r(all_id_unique)", valid_time)
    }
end





capture program drop time_section
program define time_section
    syntax [,label(string asis) clear list off timeit]
    
    * If the timeit option is not fed into this sub-program, it exits without doing anything
    if mi("`timeit'") {
        exit
    }

    if !mi("`list'") {
        timer off 1
        qui timer list
        local total_timer_val = r(t1)
        di as result "{hline 60}"
        forval x = 2/$timer_num {
            local tnum = `x' - 1
            local timer_val = r(t`x')
            local percentage = (`timer_val'/`total_timer_val') * 100
            local cum_percentage = `cum_percentage' + `percentage'
            di as text "`tnum': ${timer_label`x'}:" _col(30) as result %9.2f `timer_val' ///
                %9.1f `percentage' "%" %9.1f `cum_percentage' "%"
        }
        di as result "{hline 60}"
        di as text "${timer_label1}:" _col(30) as result %9.2f `total_timer_val'
        macro drop timer_label* timer_num
    }

    else if !mi("`clear'") {
        timer clear
        macro drop timer_* 
    }

    else if !mi("`off'") {
        timer off $timer_num
    }

    else if mi("${timer_num}") {
        timer on 1
        global timer_label1 "Total"
        global timer_num = 2
        global timer_label${timer_num} "`label'"
        timer on $timer_num
    }

    else {
        timer off $timer_num 
        global timer_num =${timer_num} + 1
        global timer_label${timer_num} "`label'"
        timer on $timer_num
    }
end
