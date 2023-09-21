*! version 1.1  21sep2023
* updated to new GLS routine - checks for missing preevent vars
* GLS uses return residuals
* correct CDF reporting 

capture program drop csestudy_092023
program define csestudy_092023, eclass
    syntax varlist [if] [in], EVentdate(string) STARTpreeventdate(string) ///
    ENDpreeventdate(string) ///
    [noBALance gls npc(integer 100) PRESAMPLEmarker(name) nobalcheck timeit]  

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


    _xt
    local panelvar = r(ivar) 
    local timevar = r(tvar)

    qui desc, short varlist
    local sortlist = r(sortlist)
    tokenize "`sortlist'"
    local sortvar1 `1'
    local sortvar2 `2'

    * Validate eventdate and start and end pre-event dates
    foreach dateinput in eventdate startpreeventdate endpreeventdate {
        * Evaluate input macro
        capture local `dateinput' = ``dateinput''
        capture confirm number ``dateinput''
        if _rc != 0 {
            di as error "Invalid value for `dateinput': ``dateinput''"
            exit _rc
        }
    }

    if `startpreeventdate'  >= `endpreeventdate' {
        di as error "Pre-event start date must come before end date."
        exit
    }

    if `eventdate'  <= `endpreeventdate' {
        di as error "Event Date come after pre-event end date."
        exit
    }

    
    tokenize `varlist'
    local lhsvar `1'
    macro shift
    local rhsvars `*'
    scalar Dim = wordcount("`rhsvars' constant")

    tempvar preevent
    qui gen `preevent' = inrange(`timevar',`startpreeventdate' , `endpreeventdate')

    time_section, label(Panel Setup) `timeit'

    * Create event variable to track event period
    tempvar event
    qui count if `timevar' == `eventdate'
    if r(N) == 0 {
        di as error "Event date does not occur in data. Check to make sure your format is correct."
        exit
    }
    gen `event' = `timevar' == `eventdate'

    marksample touse
    qui replace `touse' = 0 if !`event'
    
    tempvar touse_pre_event
    tempvar missing
    egen `missing' = rowmiss(`varlist')
    qui gen `touse_pre_event' = `preevent' & `missing' == 0


    * Generate max number of trading dates assuming no gaps


    if mi("`balance'") {
        time_section, label(Balancing Routine) `timeit'
        if "`sortvar1' `sortvar2'" != "`panelvar' `timevar'" {
            display "NOTE: Your data is not sorted. Sorting..."
            sort  `panelvar' `timevar' 
        }
        * Get min and max pre-dates
        tempvar tcount
        qui by `panelvar' (`timevar'): gen `tcount' = sum(`preevent')
        qui by `panelvar' (`timevar'): replace `tcount' = . if _n != _N
        sum `tcount', meanonly
        local n_pre_dates = r(max)
        
        tempvar in_sample_test full_obs_test
        
        qui by `panelvar' (`timevar'): gen `in_sample_test' = cond(sum(`touse')>0,1,0)
        qui by `panelvar' (`timevar'): gen `full_obs_test' = sum(`touse_pre_event')
        
        * Mark out pre_event sample which are not in sample during the event period
        qui by `panelvar' (`timevar'): replace `touse_pre_event' = 0 if `in_sample_test'[_N] != 1 | `full_obs_test'[_N] != `n_pre_dates'

        * Mark out the test sample if there is not a full set of obserations
        * for the pre-period
        qui by `panelvar' (`timevar'): replace `touse' = 0 if `full_obs_test'[_N] != `n_pre_dates' | `in_sample_test'[_N] != 1

        * Balance check if every panel-id exists for every date in touse_pre_event
        if mi("`balcheck'") {
            time_section, label(Check Balance) `timeit'
            mata is_balanced("`panelvar'", "`timevar'","","`touse_pre_event'")
            if !r(all_id_unique) {
                di as smcl as err "Pre-period data cannot be balanced. Check for unusual gaps in data"
                exit
            }
        }
        // I changed this sort to match unbalanced data will make _ts_regress faster?
        
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

    sort `timevar' `panelvar' 
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
    di as result "Note: The e(V) matrix is reverse-engineered from the CDF p-values"
    di as result "to be compatible with table output and should not be used literally."

    ereturn post `b' `V', depname("`lhsvar'") esample(`touse')
    ereturn local  cmd  "csestudy"
    ereturn local  cmdline  "`cmdline'"
    ereturn scalar N = `nobs'
    ereturn scalar n_pe_dates = `ndates'
    ereturn matrix ts_z = `ts_z'
    ereturn matrix pcdf = `pcdf'
    ereturn matrix betas = `betas'

    if !mi("`presamplemarker'") {
        gen byte `presamplemarker' = `touse_pre_event'
        la var `presamplemarker' "Pre-Event Sample Used"
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
        real matrix residuals, XX, Xy, M, E, y, X
        real scalar obs
        real vector beta, pre_event_dates, event_beta 

        // Setup Matrices for pre and post period
        st_view(E, . , (time_column,panel_column,regcols), event_sample_col)
        // If nobalance option selected and data is not already sorted by time,
        //  import as a data matrix. Otherwise import as a view
        // 09/20/23 need to think about this
        if (nobalance == "nobalance" & st_local("sortvar1") != timevar) {
            M = st_data(. , (time_column,panel_column,regcols), pre_event_sample_col)
            M = sort(M,(1,2))
        }
        else {
            //st_view(M, . , (time_column,panel_column,regcols), pre_event_sample_col)
            M = st_data(. , (time_column,panel_column,regcols), pre_event_sample_col)
        }

        pre_event_dates = uniqrows(M[.,1])
        firms = rows(uniqrows(M[.,2]))
        ndates = rows(pre_event_dates)
        pre_event_beta = J(ndates, 1 + cols(regcols), .)
        residuals = J(ndates, firms, .)
        
        //calculate pre-event betas
        startobs = 1
        i = 1
        last_obs_flag = 0

        //This works for balanced and for non-balanced
        for (obs=1; obs <= rows(M); obs++) {

            if (obs == rows(M)) last_obs_flag = 1
            else if (M[obs,1] != M[obs+1,1]) last_obs_flag = 1
            
            if (last_obs_flag) {
                    st_subview(y, M, (startobs,obs), 3)
                    st_subview(X, M, (startobs,obs), 4\.)
                    X = X,J(rows(X),1,1)
                    XX = quadcross(X,X)
                    Xy = quadcross(X,y)
                    // ------------ COMPUTE COEFFICIENTS --------------------
                    beta    = cholsolve(XX,Xy)
                    // ------------ STORE OUTPUT ----------------------------
                    st_subview(curdate,M,obs,1)
                    pre_event_beta[i,.] = curdate, beta'

                    if (gls == "gls") {    
                        residuals[i,.] = (y-X*beta)'
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
            Omega =  makesymmetric(pca_coeff * diag(variance(pca_score)) * pca_coeff' + diag(sig2_e))
            
            L = cholesky(Omega)


            //redo pre-event betas using GLS
            startobs = 1
            i = 1
            last_obs_flag = 0
            
            for (obs=1; obs <= rows(M); obs++) {
               
                if (obs == rows(M)) last_obs_flag = 1
                else if (M[obs,1] != M[obs+1,1]) last_obs_flag = 1
                
                if (last_obs_flag) {
                
                    st_subview(y, M, (startobs,obs), 3)
                    st_subview(X, M, (startobs,obs), 4\.)
                    X = X,J(rows(X),1,1)             
                    X = solvelower_wrapper(L,X)
                    y = solvelower_wrapper(L,y)
                    XX = quadcross(X,X)
                    Xy = quadcross(X,y)
                    // ------------ COMPUTE COEFFICIENTS --------------------
                    beta    = cholsolve(XX,Xy)
                    // ------------ STORE OUTPUT ----------------------------
                    //st_subview(curdate,M,obs,1)
                    pre_event_beta[i,.] = M[obs,1], beta'
                    startobs = obs + 1
                    last_obs_flag = 0
                    i = i + 1
                    
                }
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
    void is_balanced( string scalar varlist, //
        string scalar panelid, string scalar timeid, //
        string scalar touse) {
        
        real matrix X, unique_ids
        real scalar valid_time, i

        st_view(X,.,(panelid, timeid, tokens(varlist)), touse)
        unique_ids = uniqrows(X[,2],1)
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
