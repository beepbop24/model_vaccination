%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 100);
k_vals = linspace(0.1, 1, 100);

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhist, options, beta_vals, k_vals, 1, 'default_betak');

%% DEFAULT RUN FOR INFECTION -- IIS MODEL
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 100);
k_vals = linspace(0.1, 1, 100);

param_contour_plot("iis", death_criteria, "default_params_full.mat", params, @ddeIIShist, options, beta_vals, k_vals, 1, 'iis_betak');


%% DEFAULT RUN FOR INFECTION WITH NO IMMUNE SYSTEM AT ALL -- IIS MODEL
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

% parameter values for I
params.k_i = 0;

% contour plot for beta, k_V
params.tspan = [0, 300];

beta_vals = linspace(10^(-3), 5*10^(-2), 100);
k_vals = linspace(0.1, 1, 100);

param_contour_plot("iis", death_criteria, "default_params_full.mat", params, @ddeIIShist, options, beta_vals, k_vals, 1, 'viral_betak');

%% HISTORY FUNCTIONS FOR VARIOUS SCENARIOS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddeIIShist(t)
    % constant history function for VXYRI.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0];
end

