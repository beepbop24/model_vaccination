
%% OPTION 1
IIS_sol_ex1 = IIS_ex1("default_params_IIS_big.mat", @ddeIIShist);

function IIS_sol_ex1 = IIS_ex1(params_file, hist)
    params_struct = load(params_file);
    params = params_struct.params;

    lag_vals = params_struct.lags;
    lags = [lag_vals.tau_1, lag_vals.tau_2, lag_vals.tau_35, lag_vals.tau_4];

    IIS_sol_ex1 = dde23(@ddeIIS, lags, hist, params_struct.tspan);

    function dydt = ddeIIS(t, y, Z)
    % Differential equations function for VXYRI.
        ylag1 = Z(:,1); %tau_1
        ylag2 = Z(:,2); %tau_2
        ylag3 = Z(:,3); %tau_35
        ylag4 = Z(:,4); %tau_4

        %condition = params.beta*params.k*params.mu-params.d_x*params.k1_tilde*params.d_v*params.d_y;

        %if condition >= 0
        %  disp('Parameter Condition Met')
        % disp(condition)
        %else
        %disp('Parameter Condition Not Met')
        %disp(condition)
        %end

        dydt = [params.k*ylag1(3)/(params.k1_tilde + ylag3(5)^params.n_1) - params.d_v*y(1)
                params.mu - params.d_x*y(2) - params.beta*y(2)*y(1) - params.k_ix*y(2)*y(5)
                params.beta*y(2)*y(1) - params.d_y*y(3) - params.k_iy*y(3)*y(5)
                params.k_ix*y(2)*y(5) + params.k_iy*y(3)*y(5) - params.d_r*y(4)
                params.k_i * ylag2(3) + (params.b_2*ylag4(5)^params.n_2)/(params.k_2^params.n_2 + ylag4(5)^params.n_2) - params.d_i*y(5)];
    end
end

function s = ddeIIShist(t)
    % constant history function for VXYRI.
    s = [0; 7*10^7; 7; 0; 0];
end
