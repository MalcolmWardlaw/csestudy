*! version 1.6  30July2025


capture program drop csestudy
program define csestudy, eclass
    syntax varlist [if], EVENTstartdate(string) ///
        FIRSTPREeventdate(string) LASTPREeventdate(string) ///
        [gls npc(real 100) coefsonly]


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
    marksample touse
    qui replace `touse' = 0 if `timevar' != `eventstartdate'

    // If GLS is specified, mark gls_window containing event and pre-period
    if "`gls'" == "gls" {
        tempvar gls_window
        mark `gls_window' if !mi(`lhsvar') & ///
            (inrange(`timevar', `firstpreeventdate', `lastpreeventdate') | ///
            `timevar' == `eventstartdate')
    }


    ****************************************************************************
    *                                 Run Test                                 *
    ****************************************************************************


    tempname b nobs

    //pre-allocate matrices
    local colnames `rhsvars' :_cons
    local ncols: word count `colnames'
    matrix `b' = J(1,`ncols',.)

    // Get event period coefficients
    mata _get_coefficients( ///
        "`varlist'", ///
        "`panelvar'", ///
        "`timevar'", ///
        "`touse'", ///
        "`gls_window'", ///
        `n_pre_event_days', ///
        `npc', ///
        "`gls'", ///
        "`b'", ///
        "`nobs'" )

    // Label beta matrix
    local colnames `rhsvars' :_cons
    matrix rownames `b' = y1
    matrix colnames `b' = `colnames'
  

    if mi("`coefsonly'") {
        tempname all_betas all_nobs
        // Allocate all_betas matrix and store 
        // event period coefficient at start of matrix
        matrix `all_betas' = J(`n_pre_event_days'+1,`ncols',.)    
        local all_betas_colnames `eventstartdate'
        forval i = `lastpreeventdate' (-1) `firstpreeventdate' {
            local all_betas_colnames `all_betas_colnames' "`i'"
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
        di "Percent complete = 0% " _continue

        // Create new touse_pre_event variable
        // Create subsample that can be quickly marked if a large
        // amount of data is retained in the Stata dataset
        marksample marked_all
        tempvar touse_pre_event estimation_window
        gen byte `touse_pre_event' = 0
        if "`gls'" == "gls" {
            // GLS needs to test data back to n_pre_events prior
            local all_data_start_date = `firstpreeventdate' - (`eventstartdate'-`firstpreeventdate')
            // GLS also requires an ID vector to test whether !mi(lhsvar)
            tempvar marked_y
            gen byte `marked_y' = !mi(`lhsvar')

        }
        else {
            local all_data_start_date = `firstpreeventdate'            
        }

        // Mark subset of data which should be tested for valid inclusion
        gen byte `estimation_window' = inrange(`timevar',`all_data_start_date',`eventstartdate')

        preserve
        qui keep if `estimation_window'
        // Loop through each pre-event date and run regression
        tempname pre_event_b pre_event_nobs
        forval noevent_date = `lastpreeventdate'(-1)`firstpreeventdate' {

            mata _set_touse("`touse_pre_event'", "`marked_all'", ///
                "`timevar'", "`estimation_window'", `noevent_date')

            if "`gls'" == "gls" {
                local noevent_lastpreeventdate = `noevent_date' - (`eventstartdate' - `lastpreeventdate')
                local noevent_firstpreeventdate = `noevent_date' - (`eventstartdate' - `firstpreeventdate')

                mata _set_gls_window(`noevent_date', ///
                    `noevent_lastpreeventdate', ///
                    `noevent_firstpreeventdate', ///
                    "`timevar'", ///
                    "`gls_window'", ///
                    "`marked_y'", ///
                    "`estimation_window'"
                )
            }
            
            mata _get_coefficients( ///
                "`varlist'", ///
                "`panelvar'", ///
                "`timevar'", ///
                "`touse_pre_event'", ///
                "`gls_window'", ///
                `n_pre_event_days', ///
                `npc', ///
                "`gls'", ///
                "`pre_event_b'", ///
                "`pre_event_nobs'" )

            local j =  `lastpreeventdate' - `noevent_date' + 2
            matrix `all_betas'[`j',1] = `pre_event_b'
            matrix `all_nobs'[`j',1] = `pre_event_nobs'

            // Report percentage of gls regressions finished
            local pe_regno =  `lastpreeventdate' - `noevent_date' + 1
            local pe_total = `lastpreeventdate' - `firstpreeventdate' + 1
            local percent_complete = `pe_regno'/`pe_total' * 100
            if `percent_complete' - `percent_complete_last' >= 10 {
                di "... " %-2.0f `percent_complete' "% " _continue
                local percent_complete_last = `percent_complete'
            }            
        }
        restore

        tempname pcdf ts_z
        mata _get_significance_stats("`all_betas'", "`pcdf'", "`ts_z'")    
        matrix rownames `pcdf' = y1 
        matrix colnames `pcdf' = `rhsvars' :_cons

        matrix rownames `ts_z' = y1
        matrix colnames `ts_z' = `rhsvars' :_cons

        di _n
        if !mi("`gls'") {
            di as text "GLS Estimates with Time Series Corrected Errors"
        }
        else {
            di as text "OLS Estimates with Time Series Corrected Errors"
        }

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
        ereturn matrix all_N = `all_nobs'
        ereturn matrix pcdf = `pcdf'
        ereturn matrix ts_z = `ts_z'
    }    

end



capture mata mata drop _get_coefficients()
mata:
    void _get_coefficients(string scalar regvar_names, ///
        string scalar panelvar_name, ///
        string scalar timevar_name, ///
        string scalar touse_name, ///
        string scalar gls_window_name, ///
        real scalar pca_window_length, ///
        real scalar num_principal_components, ///
        string scalar gls_flag, ///
        string scalar b_macro, ///
        string scalar nobs_macro) {

        real matrix AllData, EventData, gls_outputs, X
        real colvector pre_event_y_rect, y_all, y, panelvar, touse, b
        real scalar nobs

        pragma unused timevar_name

        if (gls_flag == "gls") {
        // If GLS is enabled, balance the panel

            pragma unset AllData
            st_view(AllData, . , (panelvar_name, touse_name, tokens(regvar_names)[1]), gls_window_name)

            pragma unset panelvar
            st_subview(panelvar, AllData,.,1)

            pragma unset touse
            st_subview(touse, AllData,.,2)

            pragma unset y_all
            st_subview(y_all, AllData,.,3)

            touse_pre_event = balance_y(panelvar,touse,y_all,pca_window_length)
            pre_event_y_rect = (colshape(select(y_all,touse_pre_event), pca_window_length))'
        }

        pragma unset EventData
        st_view(EventData,., (tokens(regvar_names)), touse_name)

        pragma unset y
        st_subview(y, EventData,.,1)

        pragma unset X
        st_subview(X, EventData,.,(2\.))
        X = X, J(rows(y), 1, 1)

        if (gls_flag == "gls") {
            gls_outputs = gls_mat(y, X, pre_event_y_rect, num_principal_components)
            y = gls_outputs[.,1]
            X = gls_outputs[., (2..cols(gls_outputs))]
        }

        b = beta_coefficients(y, X)
        nobs = rows(y)

        // Post the coefficients
        st_matrix(b_macro, b')
        st_numscalar(nobs_macro, nobs)
    }     
end

capture mata mata drop _get_significance_stats()
mata:
    void _get_significance_stats(string scalar betas_mat, ///
        string scalar pcdf_mat, ///
        string scalar ts_zmat) {


        betas = st_matrix(betas_mat)
        event_coefs = betas[1,.]        
        pre_event_coefs = betas[2..rows(betas),.]
        sd = sqrt(diagonal(quadvariance(pre_event_coefs)))'
        mean_coefs = mean(pre_event_coefs)
    
        event_pctile = (colsum(abs(betas:-mean_coefs):>= abs(event_coefs:-mean_coefs))):/(rows(betas))


        ts_z =  abs(event_coefs - mean_coefs):/ (sd :* sqrt((rows(betas))/(rows(betas)-1)))
        ts_z = 2:*ttail(rows(betas)-2, ts_z)

        st_matrix(ts_zmat, ts_z)
        st_matrix(pcdf_mat, event_pctile)
    }
end



capture mata mata drop beta_coefficients()
mata:
    real colvector beta_coefficients(real colvector y, real matrix X) {
        // allocate matrices
        real matrix XX, Xy
        real colvector beta

        XX = quadcross(X,X)
        Xy = quadcross(X,y)
        beta = cholsolve(XX,Xy)
        return(beta)
    }
end


capture mata mata drop CholOmega()
mata:
    real matrix CholOmega(real matrix pre_event_y_rect, real scalar npc) {
        real matrix A, U, Vt, pca_coeff, pca_score
        real matrix sig2_e, Omega
        real vector s
        A = pre_event_y_rect :- mean(pre_event_y_rect)
        
        fullsvd(A,U,s,Vt)

        pca_coeff = Vt'[,1..npc]
        pca_score = A*pca_coeff
        sig2_e = variance(A - pca_score*pca_coeff')
        // Omega should be symmetric, but Mata doesn't recognize that it is
        // so we coerce it to be symmetric 
        Omega =  makesymmetric(pca_coeff * diag(variance(pca_score)) * // 
            pca_coeff' + diag(sig2_e))

        return (cholesky(Omega))
    }
end



capture mata mata drop _set_touse()
mata:
    void _set_touse(string scalar touse_name, ///
        string scalar marked_all_name, ///
        string scalar timevar_name, ///
        string scalar estimation_window_name, ///
        real scalar eventstartdate) {
        
        real matrix AllData

        st_view(AllData=., . , (touse_name, marked_all_name, timevar_name), ///
        estimation_window_name)
        AllData[.,1] = AllData[.,2] :* (AllData[.,3]:==eventstartdate)
    }
end


capture mata mata drop _set_gls_window()
mata:
    void _set_gls_window(real scalar noevent_date, ///
        real scalar noevent_lastpreeventdate, ///
        real scalar noevent_firstpreeventdate, ///
        string scalar timevar_name, ///
        string scalar gls_window_name, ///
        string scalar marked_y_name, ///
        string scalar estimation_window_name) {
        
        real matrix AllData


        st_view(AllData=., . , ///
            (gls_window_name, marked_y_name, timevar_name), ///
            estimation_window_name)

        AllData[.,1] = AllData[.,2] :* (AllData[.,3]:==noevent_date :| (AllData[.,3]:<=noevent_lastpreeventdate :& AllData[.,3]:>=noevent_firstpreeventdate ))
    }
end




capture mata mata drop balance_y()
mata: 

    real colvector balance_y ( ///
        real colvector panelvar, ///
        real colvector touse, ///
        real colvector y_all, ///
        real scalar pca_window_length) {

        // This function balances the pre-event PCA y data
        // edits the touse vector and returns a selector for the pre-event data

        // Allocate matrices and scalars
        real colvector start_position_index, end_position_index, touse_pre_event
        real scalar i, nobs, invalid_panel

        // Initialize touse_pre_event
        touse_pre_event = J(rows(touse), 1, 0)

        // Since the observations are pre-sorted, we can use the
        // take the last element for each paneld
        end_position_index = selectindex((panelvar[1::rows(panelvar)-1] :!= panelvar[2::rows(panelvar)]) \ 1)

        // start_position_index --> first element of each panelid
        start_position_index =  (1 \ end_position_index[|1 \ rows(end_position_index)-1|] :+ 1)

        // Reset all touse and touse_pre_event values to 0 where there
        // are not a full number of observations
        for (i = 1; i <= rows(end_position_index); i++) {
            nobs = end_position_index[i] - start_position_index[i] + 1
            invalid_panel = 0
            // If the number of observations is less than the pca_window_length + 1
            if (nobs < pca_window_length + 1) invalid_panel = 1
            // Or if the last observation in the panel is marked out
            if (touse[end_position_index[i]] == 0) invalid_panel = 1
            // Or if all y observations are 0
            if (sum(y_all[start_position_index[i]::end_position_index[i]]) == 0) invalid_panel = 1
            if (invalid_panel == 1) {
                // Reset the touse values. May already be 0 if event_date obs is missing
                touse[end_position_index[i]] = 0
            }
            else {
                // Set the touse_pre_event values to 1
                touse_pre_event[start_position_index[i]::end_position_index[i]-1] = J(nobs-1,1,1)
            }
        }
        // return touse_pre_event, which selects the valid pre-event y
        return(touse_pre_event)
    }
end


capture mata mata drop gls_matrices()
mata:
    struct gls_matrices {
        real matrix X
        real colvector y
    }
end

capture mata mata drop gls_mat()
mata:
    real matrix gls_mat (
        real colvector y, ///
        real matrix X, ///
        real matrix pre_event_y_rect, ///
        real scalar npc
        ) {
        
        L = CholOmega(pre_event_y_rect, npc)
        y = solvelower_wrapper(L,y)
        X = solvelower_wrapper(L,X)
        
        return(y,X)
    }
end

if c(stata_version) >=17 {
    capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return(solvelowerlapacke(A,B))
}
else {
    capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return (solvelower(A,B))
}

