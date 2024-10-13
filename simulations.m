%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings);

%% DEFAULT RUN FOR FIRST INFECTION WITH NO PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
full_sol = full_model("default_params_full.mat", params, @ddefullhist0, options, settings);

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO B CELL COMPARTMENT
clearvars;
params = struct();
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

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
full_noB_sol = full_model("default_params_full.mat", params, @ddefullhistnoB, options, settings);

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO B CELL COMPARTMENT
clearvars;
params = struct();
settings = struct();

% parameter values for T_E
params.k_te = 0;
params.k_tetm = 0;
params.d_tepi = 0; 
params.d_te = 0;
params.d_tmpi = 0;

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
full_noT_sol = full_model("default_params_full.mat", params, @ddefullhistnoT, options, settings);

%% DEFAULT RUN FOR INFECTION WITH PREVIOUS IMMUNITY -- NO IFN MODEL
clearvars;
params = struct();
settings = struct();

% parameter values for I
params.k_i = 0;
params.k_ite = 0;
params.k_ith = 0;

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
full_noT_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings);

%% DEFAULT RUN FOR INFECTION WITH NO PREVIOUS IMMUNITY -- IIS MODEL
clearvars;
params = struct();

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, true);

%% DEFAULT RUN FOR INFECTION WITH NO IMMUNE SYSTEM AT ALL -- IIS MODEL
clearvars;
params = struct();

% parameter values for I
params.k_i = 0;

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
IIS_sol = IIS_model("default_params_full.mat", params, @ddeIIShist, options, true);

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