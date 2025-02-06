% simulations for various initial conditions
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
settings = struct();

full_sol1 = full_model("default_params_full.mat", params, @ddefullhist1, options, settings, 'init1_model');
full_sol2 = full_model("default_params_full.mat", params, @ddefullhist2, options, settings, 'init2_model');

settings.manual = false;
full_sol3 = full_model("default_params_full.mat", params, @ddefullhist3, options, settings, 'initv_model');

settings.manual = false;
full_sol4 = full_model("default_params_full.mat", params, @ddefullhistjump, options, settings, 'initjump_model');

%% HISTORY FUNCTIONS
function s = ddefullhist(t)
    % constant history function for full model.
    s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
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

function s = ddefullhistjump(t)
    if t==0
        s = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
    else
        s = [0; log10(5.25*10^9); 0; 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
    end
end