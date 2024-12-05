% initial conds
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
%settings.plot = false;

%full_sol1 = full_model("default_params_full.mat", params, @ddefullhist1, options, settings, 'init1_model');
%full_sol2 = full_model("default_params_full.mat", params, @ddefullhist2, options, settings, 'init2_model');

%settings.manual = false;
%full_sol3 = full_model("default_params_full.mat", params, @ddefullhist3, options, settings, 'initv_model');


full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'init1_model');
xvals = linspace(full_sol.x(1), full_sol.x(end), 1000);
yvals = deval(full_sol, xvals);
beta = 2.175*10^(-2);
plot(yvals(1,:).*yvals(2,:))
hold on
plot(yvals(1,:))
hold on
plot(yvals(2,:))

%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();
%settings.plot = false;

full_sol = full_model("default_params_full.mat", params, @ddefullhist, options, settings, 'test');

%% DEFAULT RUN -- NO B CELL COMPARTMENT
%clearvars;
params = struct();
settings = struct();
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
full_noB_sol = full_model("default_params_full.mat", params, @ddefullhistnoB, options, settings, 'noB_model');


%%
default_params_full

xvals = linspace(full_sol.x(1), full_sol.x(end), 10000);
yvals = deval(full_sol, xvals);
new_cells = 10.^(defaults.beta*(yvals(1,:).*yvals(2,:)));

xvalsB = linspace(full_noB_sol.x(1), full_noB_sol.x(end), 10000);
yvalsB = deval(full_noB_sol, xvals);
new_cellsB = 10.^(defaults.beta*(yvalsB(1,:).*yvalsB(2,:)));

figure();
%semi(yvals(1,:))
hold on
%semilogy(10.^yvals(2,:))
%hold on
semilogy(new_cells)
hold on
semilogy(new_cellsB)

%% HISTORY FUNCTIONS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddefullhistnoT(t)
    % constant history function for full model -- no T
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 0; log10(90.5); 0; log10(30);];
end

function s = ddefullhistnoB(t)
    % constant history function for full model -- no B
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; 0; 0; 0];
end

function s = ddefullhist1(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5.25); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddefullhist2(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5.25*10^6); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end

function s = ddefullhist3(t)
    % constant history function for full model.
    s = [log10(7.5*10^2); log10(5.25*10^9); 0; 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
end