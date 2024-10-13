%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings);

% evaluating output for 1000 evenly-spaced points
xvals = linspace(full_sol.x(1), full_sol.x(end), 1000);
yvals = deval(full_sol, xvals);

s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];

function dydt = ddefull(y, params, hist, lags)
    for i=1:1000
        nlags = length(lags);
        Z = zeros(nlags, 1);
        for k = 1:nlags
            if lags(k) > i  
                Z(k) = hist(k);
            else
            Z(k) = y(k, i);
            end
        end

        % default computation for Y 
        y_infected = params.beta*y(2, i)*y(1, i) - params.d_y*y(3, i) - params.k_iy*y(3, i)*y(5, i)-params.k_tey*y(3, i)*y(7, i);
    
        %% threshold for T and B cells
        % default values below threshold
        t_h = -params.d_th*y(6, i);
        t_e = -params.d_te*y(7, i);
        t_m = 0;

        b_ll = 0;
        b_e = -params.d_be*y(10, i);

        % threshold conditions -- T cells
        if y(1) > params.v_star
            t_h = t_h + params.k_th*y(5, i)*y(3, i);
            t_e = t_e + params.k_tetm*y(8, i)*y(3, i);
            t_m = t_m - params.k_tetm*y(8, i)*y(3, i);
        end

        if ylag11(1, i) > params.v_star
            t_e = t_e + params.k_te*ylag11(3)*(1-y(7)/(params.t_total-y(8, i)))*y(5, i)*y(6, i); 
        end

        % generation of memory T cells post-infection
        if y(3) < settings.tolerance
            t_e = t_e - params.d_tepi*y(7, i) - params.d_tmpi*y(7, i);
            t_m = params.d_tmpi*y(7);
  
            if settings.manual == true && y(6, i) < settings.tolerance
                t_h = -y(6, i);
            end

            if settings.manual == true && y(7, i) < settings.tolerance
                t_e = -y(7, i);
            end
        end
        
        % threshold conditions -- B cells
        if ylag8(1, i) > params.v_star
            b_ll = b_ll + params.k_bll*ylag8(1, i)*ylag8(6,i);
        end

        if ylag9(1, i) > params.v_star
            b_ll = b_ll - params.k_bllbe*ylag9(1, i)*ylag9(9, i);
            b_e = b_e + params.k_bllbe*ylag9(1, i)*ylag9(9, i);
        end

        if ylag10(1, i) > params.v_star
            b_e = b_e + params.k_be*ylag10(1, i)*ylag10(6, i);
        end

        dydt = [params.k*ylag1(3, i)/(params.k1_tilde^params.n_1 + ylag3(5, i)^params.n_1) - params.d_v*y(1, i) - params.rhoV*y(1, i)*y(11, i)
                params.mu - params.d_x*y(2) - params.beta*y(2)*y(1) - params.k_ix*y(2)*y(5)
                y_infected
                params.k_ix*y(2, i)*y(5, i) + params.k_iy*y(3, i)*y(5, i) - params.d_r*y(4, i)
                params.k_i * ylag2(3, i) + (params.b_2*ylag4(5, i)^params.n_2)/(params.k_2^params.n_2 + ylag4(5, i)^params.n_2) + params.k_ite*ylag6(7, i) + params.k_ith*ylag7(6, i) - params.d_i*y(5, i)
                t_h
                t_e
                t_m
                b_ll
                b_e
                params.k_blla*y(9, i) + params.k_bea*y(10, i) - params.rhoA*y(1, i)*y(11, i)- params.d_a*y(11, i)];
    end
end

    function variation_p = sensitivity_analysis(x1, x2, delta_p)
    delta_x = abs(x1-x2)/(x1);
    variation_p = delta_x/delta_p;
end
