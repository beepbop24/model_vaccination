%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 20);
k_vals = linspace(10^(-3), 1, 20);

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhist, options, beta_vals, k_vals, 2, 'default_betak_2');

%% DEFAULT RUN -- NO B CELL COMPARTMENT
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% parameter values for B_LL
params.k_bll = 0; 
params.k_bllbe = 0;
% parameter values for B_E
params.k_be = 0;
params.d_be = 0;
% parameter values for A
params.d_a = 0;
params.k_blla = 0;
params.k_bea = 0;

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 20);
k_vals = linspace(10^(-3), 1, 20);

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhistnoB, options, beta_vals, k_vals, 2, 'noB_betak');

%% DEFAULT RUN -- NO T CELL COMPARTMENT
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% parameter values for T_E
params.k_te = 0;
params.k_tetm = 0;
params.d_tepi = 0; 
params.d_te = 0;
params.d_tmpi = 0;

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 20);
k_vals = linspace(10^(-3), 1, 20);

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhistnoT, options, beta_vals, k_vals, 2, 'noT_betak');

%% DEFAULT RUN -- NO IFN MODEL
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% parameter values for I
params.k_i = 0;
params.k_ite = 0;
params.k_ith = 0;

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 20);
k_vals = linspace(10^(-3), 1, 20);

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhist, options, beta_vals, k_vals, 2, 'noIFN_betak');

%% HISTORY FUNCTIONS FOR VARIOUS SCENARIOS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddefullhistnoB(t)
    % constant history function for full model -- no B
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; 0; 0; 0];
end

function s = ddefullhistnoT(t)
    % constant history function for full model -- no T
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 0; log10(90.5); 0; log10(30);];
end

