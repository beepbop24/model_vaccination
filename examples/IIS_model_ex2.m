%% OPTION 2
params = struct();

% the user can specify certain values or none at all
params.k = 5;

IIS_sol_ex2 = IIS_ex2(params, @ddeIIShist);

function IIS_sol_ex2 = IIS_ex2(params, hist)

    % parameter values for V
    defaults.k = 2000; % 
    defaults.k1_tilde = 500; %*
    defaults.n_1 = 1.875; %*
    defaults.d_v = 2;

    % parameter values for X
    defaults.mu = 2.625*10^9*log(2)/4320;
    defaults.d_x = 37.5*log(2)/4320;
    defaults.beta = 0.000000009375;
    defaults.k_ix = 2.5*10^(-11); 

    % parameter values for Y
    defaults.d_y = 1/12;
    defaults.k_iy = 8*10^(-11);

    % parameter values for R
    defaults.d_r = 37.5*log(2)/4320;

    % parameter values for I
    defaults.k_i = 1;
    defaults.b_2 = 0;
    defaults.k_2 = 0.5;
    defaults.n_2 = 1;
    defaults.d_i = 0.5;

    % lag values
    defaults.tau_1 = 10;
    defaults.tau_2 = 6;
    defaults.tau_35 = 12; % estimated
    defaults.tau_4 = 9; % estimated

    % duration of simulation
    defaults.tspan = [0, 20];


    % setting param values with default if not specified in user input
    default_names = fieldnames(defaults);
    param_names = fieldnames(params);
    missing = find(~ismember(default_names, param_names));

    for k = 1:length(missing)
        params.(default_names{missing(k)}) = defaults.(default_names{missing(k)});
    end

    lags = [params.tau_1, params.tau_2, params.tau_35, params.tau_4];

    IIS_sol_ex2 = dde23(@ddeIIS, lags, hist, params.tspan);

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