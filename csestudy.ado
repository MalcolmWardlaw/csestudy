*! version 0.1  31oct2022
program define csestudy,  eclass
    syntax varlist [if] [in], EVentdate(string) STARTpreeventdate(string) ///
    ENDpreeventdate(string) ///
    [tsols npc(integer 100) PRESAMPLEmarker(name) nobalcheck timeit]  

    local cmdline "csestudy `0'"

    time_section, clear `timeit'
    time_section, label(Setup) `timeit'

    capture confirm variable `presamplemarker'
    if _rc ==0 {
        di as error "Pre-sample marker variable name (`presamplemarker') already exists."
        exit
    }


    qui tsset, noquery
    local panelvar = r(panelvar) 
    local timevar = r(timevar)

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
    qui replace `touse' = 0 if !`event' == 1
    
    tokenize `varlist'
    local lhsvar `1'
    macro shift
    local rhsvars `*'
    scalar Dim = wordcount("`rhsvars' constant")


    * Get min and max pre-dates
    tempvar tcount
    qui by `panelvar' (`timevar'): gen `tcount' = sum(`preevent')
    qui by `panelvar' (`timevar'): replace `tcount' = . if _n != _N
    sum `tcount', meanonly
    local n_pre_dates = r(max)
    
    time_section, label(Balancing Routine) `timeit'
    
    * Generate max number of trading dates assuming no gaps
    tempvar touse_pre_event
    mark `touse_pre_event' if !mi(`lhsvar') & `preevent'

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

    qui count if `touse'
    local nobs = r(N)

    qui count if `touse_pre_event'
    local ndates = r(N)/`nobs'

    **********************************************************************
    *                              GLS                                   *
    **********************************************************************
    
    * Turn gls method on or off
    if mi("`tsols'") {
        time_section, label(GLS Setup) `timeit'
        tempname b V

        * Check that the number of pca components is less than the number of observations
        capture assert `nobs' > `npc'
        if _rc {
            di as error "After balancing pre-period, you have more pca components than test observations."
            exit
        }

        di as text "Running PCA and GLS estimation..." _n

        time_section, label(Mata Code) `timeit'
        mata: _tsgls("`varlist'", "`touse_pre_event'", "`touse'", `nobs', `ndates',"`b'", "`V'", `npc')
        local cnames `rhsvars' _cons
        matrix colnames `b' = `cnames'
        matrix colnames `V' = `cnames'
        matrix rownames `V' = `cnames'
        time_section, label(Post results to e) `timeit'

        ereturn clear
        ereturn post `b' `V', depname("`lhsvar'") esample(`touse')
        ereturn local  cmd     "csestudy"
        ereturn local  cmdline  "`cmdline'"
        ereturn scalar N = `nobs'
        ereturn scalar n_pe_dates = `ndates'
        ereturn scalar npc = `npc'

        di as text "Time Series GLS Estimates"
        di _col(53) as text "Number of obs  = " as result %9.0fc `nobs'
        di _col(41) as text "Number of pre-period dates = " as result %9.0fc `ndates'
        di _col(37) as text "Number of Principal Components = " as result %9.0fc `npc'
        
        ereturn display
    }


    **********************************************************************
    *                        Time Series OLS                             *
    **********************************************************************
    else {
        time_section, label(Preserve) `timeit'
        preserve
        qui keep if `touse' | `touse_pre_event' 

        time_section, label(regressby setup) `timeit'

        * Sort by time variable
        sort `timevar'
        * Generate a single counting variable for time
        tempvar grp
        gen `grp' = sum(`timevar'!=`timevar'[_n-1])

        di as text "Estimating TS-OLS error structure..." _n
        
        time_section, label(Mata Code Regressby) `timeit'
        * Perform regressions on each by-group, store in dataset
        mata: _regressby("`varlist'", "`grp'", "`timevar'")

        time_section, label(Mata Cleanup) `timeit'

        rename cons constant

        gen byte `event' = `timevar'==`eventdate'
        time_section, label(mvreg) `timeit'
        
        qui mvreg `rhsvars' constant = `event'
        time_section, label(Get information for post) `timeit'
        * Get Results
        tempname b V ts_z pcdf ts_se tsBeta
        mkmat `rhsvars' constant if `event'==1, matrix(`b')
        * Get the odd columns which are coefficients
        mata st_matrix("`tsBeta'",select(st_matrix("e(b)"),mod(1..cols(st_matrix("e(b)")),2)))
        mata st_matrix("`ts_se'", sqrt(diagonal(select(select(st_matrix("e(V)"),mod(1..cols(st_matrix("e(V)")),2)),mod(1::cols(st_matrix("e(V)")),2))))')

        * Assign new names
        local eqnames: coleq e(b)
        local bcolnames: list uniq eqnames

        * local rhsvars w_ltdd w_me w_btm
        matrix rownames `b' = y1
        matrix colnames `b' = `rhsvars' :_cons
        matrix rownames `ts_se' = y1
        matrix colnames `ts_se' = `rhsvars' :_cons



        * Calculate alternate Z-Score
        scalar Dim = wordcount("`rhsvars' constant")
        matrix `ts_z' = J(1,Dim,0)
        matrix rownames `ts_z' = y1
        matrix colnames `ts_z' = `rhsvars' :_cons
        matrix `pcdf' = J(1,Dim,0)
        matrix rownames `pcdf' = y1
        matrix colnames `pcdf' = `rhsvars' :_cons

        foreach coef in `rhsvars' constant {
            qui sum `coef' if !`event'
            local mu_os_beta = r(mean)
            local mu_os_sd = r(sd)
            qui sum `coef' if `event', meanonly
            local mu_event_beta = r(mean)
            local zscore = (`mu_event_beta' - `mu_os_beta')/`mu_os_sd'
            * Handle constant
            if "`coef'" == "constant" {
                local colnm "_cons" 
            }
            else {
                local colnm `coef'      
            }
            matrix `ts_z'[1, colnumb(`ts_z',"`colnm'") ] = `zscore'
        }

        * Calculate finite sample adjustment to avoid deterministic p-val 
        qui count if !`event'
        local fsadj = (1/r(N))/2
        tempvar pctile
        foreach coef in `rhsvars' constant {
            qui sum `coef' if `event', meanonly
            local eventbeta = r(mean)
            qui sum `coef' if !`event', meanonly
            local samplebeta = r(mean)

            qui egen `pctile' = mean(abs(`coef' - `samplebeta') > abs(`eventbeta' - `samplebeta')) if !`event' & !mi(`coef')
            qui sum `pctile', meanonly
            local pval = cond(r(mean) > .5, r(mean) - `fsadj', r(mean) + `fsadj')
            * Handle constant
            if "`coef'" == "constant" {
                local colnm "_cons" 
            }
            else {
                local colnm `coef'      
            }
            matrix `pcdf'[1, colnumb(`pcdf',"`colnm'") ] = `pval'
            drop `pctile'
        }


        * Dark Magic
        matrix `V' = J(Dim,Dim,0)
        mata st_matrix("`V'", diag(-abs(st_matrix("`b'")):*invnormal(st_matrix("`pcdf'")/2):^-1):^2)
        matrix rownames `V' = `rhsvars' :_cons
        matrix colnames `V' = `rhsvars' :_cons

        time_section, label(restore) `timeit'
        restore

        time_section, label(Post results to e) `timeit'
        ereturn clear

        di as text "OLS Estimates with Time Series Corrected Errors"
        di _col(36) as text "Number of obs  = " as result %9.0fc `nobs'
        di _col(24) as text "Number of pre-period dates = " as result %9.0fc `ndates'

        di as text "{hline 13}{c TT}{hline 47}"
        di as text %12s abbrev("`lhsvar'",12)  " {c |}  Coefficient" _col(29) %~12s  "CDF p-val" _col(41)  %~12s  "TS Z-Score" _col(53) %~12s  "TS-SE" 
        di as text "{hline 13}{c +}{hline 47}"
        foreach colnm in `rhsvars' _cons {
            di as text %12s abbrev("`colnm'",12) " {c |}" ///
            _col(17) as result %9.0g `b'[1, colnumb(`b',"`colnm'") ] ///
            _col(29) %9.3f `pcdf'[1, colnumb(`pcdf',"`colnm'") ] ///
            _col(41) as result %9.3f `ts_z'[1, colnumb(`ts_z',"`colnm'") ] ///
            _col(53) %9.0g `ts_se'[1, colnumb(`ts_se',"`colnm'") ] ///

        }
        di as text "{hline 13}{c BT}{hline 47}" _n
        di as result "Note: The e(V) matrix is reverse-engineered from the CDF p-values"
        di as result "to be compatible with table output and should not be used literally."


        ereturn post `b' `V', depname("`lhsvar'") esample(`touse')
        ereturn local  cmd  "csestudy"
        ereturn local  cmdline  "`cmdline'"
        ereturn scalar N = `nobs'
        ereturn scalar n_pe_dates = `ndates'
        ereturn matrix ts_z = `ts_z'
        ereturn matrix pcdf = `pcdf'
        ereturn matrix ts_se = `ts_se'

    }

    if !mi("`presamplemarker'") {
        gen byte `presamplemarker' = `touse_pre_event'
        la var `presamplemarker' "Pre-Event Sample Used"
    }

    time_section, off `timeit'
    time_section, list `timeit'
end


*-------------------------------------------------------------------------------
* Mata program: _regressby
* Inputs:
*   - A y-var and x-var for an OLS regression
*   - A group var, for which each value represents a distinct by-group. 
*       This var MUST be in ascending order.
*   - A list of numeric by-variables, whose groups correspond to the group var.
* Outputs:
*   - dataset of coefficients from OLS regression for each by-group
*-------------------------------------------------------------------------------

version 13.1
mata:
    void _regressby(string scalar regvars, string scalar grpvar, ///
        string scalar byvars) {

        // Convert variable names to column indices
        real rowvector regcols, bycols
        real scalar grpcol
        regcols     = st_varindex(tokens(regvars))
        bycols      = st_varindex(tokens(byvars))
        grpcol      = st_varindex(grpvar)

        // Fetch number of groups
        real scalar numgrp, startobs, curgrp
        numgrp      = _st_data(st_nobs(),grpcol)
        startobs    = 1  
        curgrp      = _st_data(1,grpcol)

        // Preallocate matrices for output
        real matrix groups, coefs, nobs
        groups      = J(numgrp, cols(bycols), .)
        coefs       = J(numgrp, cols(regcols), .)
        nobs        = J(numgrp, 1, .)

        // Preallocate regression objects
        real matrix XX, Xy, XX_inv, M, y, X
        real scalar N, k, obs
        real vector beta
        string vector covariates
        string scalar covName
        // -----------------------------------------------------------------------------
        // Iterate over groups
        // -----------------------------------------------------------------------------

        // Iterate over groups
        for (obs=1; obs<=st_nobs(); obs++) {
            if (_st_data(obs+1,grpcol) != curgrp) {
                st_view(M, (startobs,obs), regcols, 0)
                st_subview(y, M, ., 1)
                st_subview(X, M, ., (2\.))
                N    = rows(X)
                // Augment x with either column of 1's or weights
                // Define matrix products
                X = X,J(N,1,1)
                XX      = quadcross(X,X)
                Xy      = quadcross(X,y)
                XX_inv  = invsym(XX)
                // ------------ COMPUTE COEFFICIENTS --------------------
                beta    = (XX_inv*Xy)'
                // ------------ STORE OUTPUT ----------------------------
                coefs[curgrp,.]     = beta
                nobs[curgrp,1]      = N
                groups[curgrp,.]    = st_data(startobs,bycols)
                // ------------ WRAP UP BY ITERATING COUNTERS -----------
                curgrp   = _st_data(obs + 1,grpcol)
                startobs = obs + 1
            }
        }  

        // -----------------------------------------------------------------------------
        // Gather output and pass back into Stata
        // -----------------------------------------------------------------------------

        // Store group identifiers in dataset
        stata("qui keep in 1/"+strofreal(numgrp, "%18.0g"))
        stata("keep " + byvars)
        st_store(.,tokens(byvars),groups)

        // Store coefficients in dataset:

        // ... Number of observations,
        (void) st_addvar("long", "N")
        st_store(., ("N"), nobs)

        // ... And then looping over covariates,
        covariates = (cols(regcols)>1) ? tokens(regvars)[|2 \ .|], "cons" : ("cons")
        for (k=1; k<=length(covariates); k++) {
            covName = covariates[k]
            // ... Coefficients and standard errors,
            (void) st_addvar("float", covName)
            st_store(., covName,  coefs[., k])
        }
    }
end


mata:
    void _tsgls(string scalar regvars, string scalar pre_event_sample, ///
        string scalar event_sample, real scalar N, ///
        real scalar ndates, string scalar bmat, string scalar Vmat, ///
        real scalar npc) {
        // Convert variable names to column indices
        real rowvector regcols
        real scalar escol, pescol
        regcols = st_varindex(tokens(regvars))
        escol = st_varindex(event_sample)
        pescol = st_varindex(pre_event_sample)
        
        // Preallocate regression objects
        real matrix return_matrix, A, U, Vt, M, X, XBX_inv, svd_ev, pca_coeff, pca_score
        real matrix sig2_e, Omega, L, VCV, tilde_X, tilde_y, XX, Xy
        real vector s, y, coef
        
        return_matrix = colshape(st_data(.,(regcols[1]),pescol),ndates)'
        return_matrix = return_matrix:-mean(return_matrix')' // De-mean returns by equal weighted daily avg
        A = return_matrix:-mean(return_matrix)
        fullsvd(A,U,s,Vt)
        svd_ev = (s[1..npc]:^2)/(ndates-1)
        pca_coeff = Vt'[,1..npc]
        pca_score = A * pca_coeff
        sig2_e = variance(return_matrix - pca_score * pca_coeff')
        // Omega should be symmetric, but Mata doesn't recognize that it is
        Omega =  makesymmetric(pca_coeff * diag(variance(pca_score)) * pca_coeff' + diag(sig2_e))
        st_view(M, . , regcols, escol)
        st_subview(y, M, ., 1)
        st_subview(X, M, ., (2\.))
        X = X,J(rows(X),1,1)
        L = cholesky(Omega)
        tilde_X = solvelower(L,X)
        tilde_y = solvelower(L,y)
        XX = quadcross(tilde_X,tilde_X)
        Xy = quadcross(tilde_X,tilde_y)
        coef = cholsolve(XX,Xy)
        VCV = cholinv(XX)
        st_matrix(bmat, coef')
        st_matrix(Vmat, VCV)

    }
end

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