*! version 1.7  6July2025


capture program drop csestudy_coefs
program define csestudy_coefs
    syntax varlist [if], firstdate(string) lastdate(string) ///
        [pewindowlength(real 200) gap(real 1) gls npc(real 100)]
    _xt, trequired
    local panelvar = r(ivar) 
    local timevar = r(tvar)

    // Evaluate eventstartdate, firstpreeventdate, and lastpreeventdate
    local firstdate = `firstdate'
    local lastdate = `lastdate'

    local n_pre_event_days = `pewindowlength'

    local cmdline "_csestudy_coefs `0'"
    
    tokenize `varlist'
    local lhsvar `1'
    macro shift
    local rhsvars `*'
    scalar Dim = wordcount("`rhsvars' constant")

    // Mark all valid event and pre-event observations and create touse
    marksample touse_all

    // Set Data View
    tempvar data_window
    gen byte `data_window' = 1

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
    
    mata all_betas_matrix = J(0,`ncols',.)
    mata all_N_matrix = J(0,1,.)
    mata all_date_matrix = J(0,1,.)

    if "`gls'" == "gls" {
        mata pre_event_date_matrix = J(0,2,.)
    }

    di "Finished date " _continue
    tempname pre_event_b pre_event_nobs
    forval pseudo_event_date = `lastdate'(-1)`firstdate' {

        local pseudo_event_lastpreeventdate = `pseudo_event_date' - `gap'
        local pseudo_event_firstpreeventdate = `pseudo_event_date' - `gap' - `n_pre_event_days' + 1

        mata _get_coefficients(long_data, full_index, `pseudo_event_date', `pseudo_event_lastpreeventdate', `pseudo_event_firstpreeventdate', `npc', "`b'", "`nobs'")
        mata all_betas_matrix = all_betas_matrix \ st_matrix("`b'")
        mata all_N_matrix = all_N_matrix \ st_numscalar("`nobs'")
        mata all_date_matrix = all_date_matrix \ (`pseudo_event_date')
        if "`gls'" == "gls" {
            mata pre_event_date_matrix = pre_event_date_matrix \ (`pseudo_event_firstpreeventdate', `pseudo_event_lastpreeventdate')
        }
        di %-2.0f `pseudo_event_date' " " _continue
    }
    
    clear
    foreach var of local rhsvars {
        local beta_varnames `beta_varnames' _b_`var'
    }
    local beta_varnames `beta_varnames' _b_cons

    getmata event_date = all_date_matrix /// 
        (`beta_varnames') = all_betas_matrix _nobs = all_N_matrix
    if "`gls'" == "gls" {
        getmata (pe_firstdate pe_lastdate) = pre_event_date_matrix
    }


end


findfile "csestudy.mata"
include "`r(fn)'"
