%% DEFAULT RUN FOR FIRST INFECTION WITH PREVIOUS IMMUNITY
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];

full_sol0 = reinfections("default_params_full.mat", params, options, model_settings, init, 10, '1over', 2);

%%