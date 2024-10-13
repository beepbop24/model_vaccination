
init = [0; log10(5.25*10^9); log10(52.5); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];

options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
params = struct();
full_sol0 = reinfections("default_params_full.mat", params, options, model_settings, init, 10, '1over', 2);

% params_file, params, options, model_settings, init, iter, variable, time