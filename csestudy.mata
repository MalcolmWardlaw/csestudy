
// capture mata mata drop _get_coefficients()
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

// capture mata mata drop _get_significance_stats()
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



// capture mata mata drop beta_coefficients()
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


// capture mata mata drop CholOmega()
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



// capture mata mata drop _set_touse()
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


// capture mata mata drop _set_gls_window()
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




// capture mata mata drop balance_y()
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


// capture mata mata drop gls_mat()
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
    // capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return(solvelowerlapacke(A,B))
}
else {
    // capture mata mata drop solvelower_wrapper()
    mata numeric matrix solvelower_wrapper(numeric matrix A, numeric matrix B) return (solvelower(A,B))
}

