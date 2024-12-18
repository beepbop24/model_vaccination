%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
%settings.plot = false;

[full_sol, I_full, M_full] = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'default_model');

%% DEFAULT RUN FOR FIRST INFECTION WITH NO PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
full_sol = full_model("default_params_full.mat", params, @ddefullhist0, options, settings,'default0_model');

% contour plot for beta, k_V
params.tspan = [0, 300];
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

beta_vals = linspace(10^(-2), 5*10^(-2), 10);
k_vals = linspace(0.1, 1, 10);
death_criteria = 0.9;

param_contour_plot("full", death_criteria, "default_params_full.mat", params, @ddefullhist0, options, beta_vals, k_vals, 1, 'default0_betak');

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO B CELL COMPARTMENT
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

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

full_noB_sol = full_model("default_params_full.mat", params, @ddefullhistnoB, options, settings, 'noB_model');

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO T CELL COMPARTMENT
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

% parameter values for T_E
params.k_te = 0;
params.k_tetm = 0;As
params.d_tepi = 0; 
params.d_te = 0;
params.d_tmpi = 0;

full_noT_sol = full_model("default_params_full.mat", params, @ddefullhistnoT, options, settings, 'noT_model');

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO IFN MODEL
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

% parameter values for I
params.k_i = 0;
params.k_ite = 0;
params.k_ith = 0;

full_noIFN_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'noIFN_model');

%% DEFAULT RUN FOR INFECTION -- IIS MODEL
clearvars;
params = struct();
params.tspan = [0, 100];

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

% UU
IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, settings, 'iis_model');

% US
params.b_2 = 9;
%IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, settings);
inset_plots("default_params_full.mat", params, @ddeIIShist, options, settings, "I star", "iis_DF", 'iis_sim_us')

% SU
params.beta = 2.05*10^(-4);
params.b_2 = 0.2;
%IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, settings);
inset_plots("default_params_full.mat", params, @ddeIIShist, options, settings, "I zero", "iis_DF", 'iis_sim_su')

% SS -- converges to I^* not 0
params.beta = 1.5*10^(-3);
params.b_2 = 0.9;
%IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, settings);
inset_plots("default_params_full.mat", params, @ddeIIShist, options, settings, "I star", "iis_DF", 'iis_sim_ss')


%% DEFAULT RUN FOR INFECTION WITH NO IMMUNE SYSTEM AT ALL -- IIS MODEL
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

% parameter values for I
params.k_i = 0;
params.tspan = [0, 100];
IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, settings, 'viral_model');

% simulation for R0 = 1 with inset
params.tspan = [0, 5000];
default_params_full
params.d_x = 10^(-1); 
params.mu = log10(5.25*10^9)*params.d_x; 
params.beta = defaults.d_v * defaults.d_y*params.d_x/(defaults.k * params.mu);

inset_plots("default_params_full.mat", params, @ddeIIShist, options, settings, "", "viral_DF", 'viral_sim_r0')


%% HISTORY FUNCTIONS FOR VARIOUS SCENARIOS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddefullhist0(t)
    % constant history function for full model -- Generation 0
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 0; 0; 0; 0];
end

function s = ddefullhistnoB(t)
    % constant history function for full model -- no B
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; 0; 0; 0];
end

function s = ddefullhistnoT(t)
    % constant history function for full model -- no T
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 0; log10(90.5); 0; log10(30);];
end

function s = ddeIIShist(t)
    % constant history function for VXYRI.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0];
end