


capture mata mata drop data_views()
capture mata mata drop get_data_views()
mata
struct data_views {
    real matrix panelids, timeids, touse, y_data, X_data
    real scalar nrows
}

struct data_views scalar get_data_views(real matrix A) {
    struct data_views scalar long_data

    // Matrix A should contain the following in each column:
    // 1. panelid
    // 2. timeid
    // 3. touse (1 if the observation is valid for on any date, 0 otherwise)
    // 4. y (the independent variable, can be missing)
    // 5-end.  X variables  (optional)

    st_subview(long_data.panelids, A, .,1)
    st_subview(long_data.timeids, A, .,2)
    st_subview(long_data.touse, A, .,3)
    st_subview(long_data.y_data, A, .,4)
    st_subview(long_data.X_data, A, .,5\.)
    long_data.nrows = rows(A)
    return(long_data)
}
end


capture mata mata drop data_indexes()
capture mata mata drop get_data_indexes()
mata:
    struct data_indexes {
    real scalar index_date, gls_flag
    real matrix data_row_index, valid_touse, valid_y, rect_y
    }

    struct data_indexes scalar get_data_indexes(struct data_views scalar long_data, 
        string scalar gls) {
        struct data_indexes scalar full
        real scalar min_time, max_time, n_time, n_panels
        real colvector a, sequential_panelids
        real matrix data_row_index, valid_touse, valid_y, rect_y

        if (gls == "gls") {
            full.gls_flag = 1
        }
        else {
            full.gls_flag = 0
        }        

        // Step 1: Basic setup
        min_time = min(long_data.timeids)
        max_time = max(long_data.timeids)
        n_time = max_time - min_time + 1

        a = 1\ (long_data.panelids[2::long_data.nrows] :!= long_data.panelids[1::long_data.nrows-1])
        sequential_panelids = runningsum(a)

        n_panels = sequential_panelids[rows(sequential_panelids)]

        data_row_index = J(n_panels, n_time, .)
        valid_touse = J(n_panels, n_time, 0)
        if (full.gls_flag == 1) {
            // If GLS is used, we need to track valid y values
            // for the pre-event period
            valid_y = J(n_panels, n_time, 0)
            rect_y = J(n_panels, n_time, .)
        }
        else {
            valid_y = J(0,0,.)
            rect_y =  J(0,0,.)
        }


        // Step 2: Fill data_row_index, valid_touse, valid_y, rect_y in a loop
        for (i = 1; i <= long_data.nrows; i++) {
            row = sequential_panelids[i]    // since panelid is 1-based sequential
            col = long_data.timeids[i] - min_time + 1
            data_row_index[row, col] = i
            if (long_data.touse[i] == 1) valid_touse[row, col] = 1
            // If GLS is used, we need to track valid y values
            if (full.gls_flag == 1) {
                rect_y[row, col] = long_data.y_data[i]
                if (long_data.y_data[i] !=.) valid_y[row, col] = 1
            }
        }
        // Step 3: Fill in the full rectangular lookups

        full.index_date = min_time
        full.data_row_index = data_row_index
        full.valid_touse = valid_touse
        full.valid_y = valid_y
        full.rect_y = rect_y

        return(full)
    }
end


capture mata mata drop current_data_indexes()
capture mata mata drop get_current_indexes()
mata:
    struct current_data_indexes {
    real matrix touse_index, pre_event_touse_index
    }

    struct current_data_indexes scalar get_current_indexes( 
        struct data_indexes scalar full, real scalar current_date,
           | real scalar pe_start_date, real scalar pe_end_date) {
        
        struct current_data_indexes scalar current
        real colvector current_valid_y, current_touse, nonzero_ys
        real rowvector col_selection
        real matrix current_rect_y
        real scalar current_col_number, pe_start_col, pe_end_col, total_cols
        
        current_col_number = current_date - full.index_date + 1

        if (full.gls_flag == 1) {
            pe_start_col = pe_start_date - full.index_date + 1
            pe_end_col = pe_end_date - full.index_date + 1
            col_selection = (pe_start_col..pe_end_col, current_col_number)
            total_cols = cols(col_selection)

            // Check if y is valid for all pre-event days
            current_valid_y = rowsum(full.valid_y[., col_selection]) :== 
                total_cols
            // check if y is 0 or close to 0 for all pre-event days
            current_rect_y = abs(full.rect_y[.,pe_start_col..pe_end_col])
            nonzero_ys = rowsum(current_rect_y) :>= .01
            current_valid_y = nonzero_ys :& current_valid_y
            
            
            current_touse = full.valid_touse[., current_col_number] :& current_valid_y
            current.touse_index = select(full.data_row_index[. , current_col_number],
                current_touse)
            current.pre_event_touse_index = vec((select(
                full.data_row_index[. , (pe_start_col..pe_end_col)],current_touse))')
        }
        else {
            current_touse = full.valid_touse[., current_col_number]
            current.touse_index = select(full.data_row_index[. , current_col_number],
                current_touse)
            
            current.pre_event_touse_index = J(0,0,.)
        }
        return(current)
    }
end



capture mata mata drop _get_coefficients()
mata:
    void _get_coefficients( real matrix A, ///
        struct data_indexes scalar full, ///
        real scalar current_date, ///
        real scalar pe_end_date, ///
        real scalar pe_start_date, ///
        real scalar num_principal_components, ///
        string scalar b_macro, ///
        string scalar nobs_macro) {
        
        struct current_data_indexes scalar current
        real matrix X, pre_event_y_rect, gls_outputs
        real colvector y, pre_event_y
        real scalar pre_event_window_length, nobs

        
        current = get_current_indexes(full, current_date, pe_start_date, pe_end_date) 
        st_subview(y=., A, current.touse_index, 4)
        st_subview(X=., A, current.touse_index, 5\.)
        X = X, J(rows(X), 1, 1)
        
        if (full.gls_flag == 1) {
            pre_event_y = A[current.pre_event_touse_index,4]
            pre_event_window_length = pe_end_date - pe_start_date + 1
            pre_event_y_rect = (colshape(pre_event_y,pre_event_window_length))'
            gls_outputs = gls_mat(y, X, pre_event_y_rect, num_principal_components)
            y = gls_outputs[.,1]
            X = gls_outputs[., (2..cols(gls_outputs))]
        }

        b = beta_coefficients(y, X)
        nobs = rows(X)

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

