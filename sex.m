%% male
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
%settings.plot = false;

% IFN
params.k_i = 0.09*1.005;
params.b_2 = 0.15*0.75;
params.tau_35 = 5+3;

params.k_th = 0.003765*0.905;
params.d_th = 0.006*0.8;
params.k_te = 0.001485*0.985;
params.k_bll = 0.001757485336*0.95;

full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'default_model');

%% male reinfections
params.tspan = [0,2000];
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.ifi = zeros(1, 10);
re_specs.type = ones(1, 10);
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

%reinfections("default_params_full.mat", params, options, settings, init, iter, re_specs)

% "random"
re_specs.cons_epi = [0 0.64 0.90 0 0 0.80 0 0 0.87 0.66];
re_specs.type = [0 1 1 0 0 1 0 0 1 1];
re_specs.ifi = [39 47 43 39 50 12 35 10 25 28];
re_specs.season = [1 0 0 0 0 0 0 0 0 0 0];
reinfections("default_params_full.mat", params, options, settings, init, iter, re_specs)


%% female
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-12, 'RelTol', 1e-12);
settings = struct();
%settings.plot = false;

% IFN
params.k_i = 0.09*0.65;
params.b_2 = 0.15*4.45;
params.tau_35 = 5-3;

params.k_th = 0.003765*1.19;
params.d_th = 0.006*1.35;
params.k_te = 0.001485*1.07;
params.k_bll = 0.001757485336*1.1;

full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'default_model');

%% female reinfections
params.tspan = [0,2000];
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.ifi = zeros(1, 10);
re_specs.type = ones(1, 10);
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

%preinfections("default_params_full.mat", params, options, settings, init, iter, re_specs)

% "random"
re_specs.cons_epi = [0 0.64 0.90 0 0 0.80 0 0 0.87 0.66];
re_specs.type = [0 1 1 0 0 1 0 0 1 1];
re_specs.ifi = [39 47 43 39 50 12 35 10 25 28];
re_specs.season = [1 0 0 0 0 0 0 0 0 0 0];
reinfections("default_params_full.mat", params, options, settings, init, iter, re_specs)

%% HISTORY FUNCTIONS FOR VARIOUS SCENARIOS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end
