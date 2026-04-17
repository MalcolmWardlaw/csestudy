


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
    real scalar pe_start_date, pe_end_date, pre_event_window_length
    }

    struct current_data_indexes scalar get_current_indexes( 
        struct data_indexes scalar full, real scalar current_date,
           | real scalar pe_end_date, real scalar pe_start_date ) {
        
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
        current.pe_start_date = pe_start_date
        current.pe_end_date = pe_end_date
        current.pre_event_window_length = pe_end_date - pe_start_date + 1
        return(current)
    }
end



capture mata mata drop _get_coefficients()
mata:
    void _get_coefficients(struct data_views scalar long_data, ///
        struct current_data_indexes scalar current, ///
        string scalar b_macro, ///
        string scalar nobs_macro, | ///
        string scalar gls_flag, ///
        real scalar num_principal_components, ///
        string scalar woodbury_flag
        ) {

        real matrix X, pre_event_y_rect, gls_outputs
        real colvector y, pre_event_y
        real scalar pre_event_window_length, nobs

        st_subview(y, long_data.y_data, current.touse_index, .)
        st_subview(X, long_data.X_data, current.touse_index, .)
        X = X, J(rows(X), 1, 1)

        nobs = rows(X)

        if (gls_flag == "gls") {
            pre_event_y = long_data.y_data[current.pre_event_touse_index]
            pre_event_y_rect = (colshape(pre_event_y, current.pre_event_window_length))'
            if (woodbury_flag == "woodbury") {
                b = gls_beta_woodbury(y, X, pre_event_y_rect, num_principal_components)
            }
            else {
                gls_outputs = gls_mat(y, X, pre_event_y_rect, num_principal_components)
                y = gls_outputs[.,1]
                X = gls_outputs[., (2..cols(gls_outputs))]
                b = beta_coefficients(y, X)
            }
        }
        else {
            b = beta_coefficients(y, X)
        }

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


capture mata mata drop gls_beta_woodbury()
mata:
    real colvector gls_beta_woodbury (
        real colvector y, ///
        real matrix X, ///
        real matrix pre_event_y_rect, ///
        real scalar npc
        ) {
        // Woodbury identity approach to GLS.
        // Computes beta = (X' Omega^{-1} X)^{-1} X' Omega^{-1} y directly
        // without forming or factoring the full N x N Omega matrix.
        //
        // Omega = V Lambda V' + D  where D = diag(sig2_e), V is N x k PCA loadings
        // Omega^{-1} = D^{-1} - D^{-1} V M V' D^{-1}
        //   where M = (Lambda^{-1} + V' D^{-1} V)^{-1}   (k x k)
        //
        // The only matrix inversion is k x k (k = npc), vs N x N for Cholesky.
        // Speedup is ~50-66% but with slightly less numerical precision due to
        // the wide range of idiosyncratic variances entering D^{-1}.

        real matrix A, U, Vt, V, pca_score
        real matrix Lambda, VtDinv, M, Oinv_X
        real colvector d, D_inv, Oinv_y, beta
        real vector s

        A = pre_event_y_rect :- mean(pre_event_y_rect)

        fullsvd(A, U, s, Vt)

        V = Vt'[, 1..npc]
        pca_score = A * V
        Lambda = diag(variance(pca_score))

        // Idiosyncratic variance (diagonal only)
        d = diagonal(variance(A - pca_score * V'))
        D_inv = 1 :/ d

        // Core Woodbury:  M = (Lambda^{-1} + V' D^{-1} V)^{-1}
        VtDinv = V' :* D_inv'
        M = invsym(invsym(Lambda) + VtDinv * V)

        // Apply Omega^{-1} to y and X:
        //   Omega^{-1} z = D^{-1} z - (D^{-1} V) M (V' D^{-1}) z
        Oinv_y = D_inv :* y - (D_inv :* V) * M * (VtDinv * y)
        Oinv_X = D_inv :* X - (D_inv :* V) * M * (VtDinv * X)

        // GLS beta = (X' Omega^{-1} X)^{-1} X' Omega^{-1} y
        beta = cholsolve(quadcross(X, Oinv_X), quadcross(X, Oinv_y))

        return(beta)
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

