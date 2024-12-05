%% reinfection immediately after primary infection, default params
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.ifi = zeros(1, 10);
re_specs.type = ones(1, 10);
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

reinfections("default_params_full.mat", params, options, model_settings, init, iter, re_specs)

%% randomizing each of the specifications independently
% random cons_epi
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
re_specs = struct();
re_specs.ifi = zeros(1, 10);
re_specs.type = ones(1, 10);
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

reinfections("default_params_full.mat", params, options, model_settings, init, iter, re_specs)

% random ifi
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.type = ones(1, 10);
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

reinfections("default_params_full.mat", params, options, model_settings, init, iter, re_specs)

% random type -- with previous = 0
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.ifi = zeros(1, 10);
re_specs.previous = 0;
re_specs.season = zeros(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

reinfections("default_params_full.mat", params, options, model_settings, init, iter, re_specs)


% random season 
clearvars;
params = struct();
options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);
model_settings = struct();
re_specs = struct();
re_specs.cons_epi = ones(1, 10);
re_specs.ifi = zeros(1, 10);
re_specs.type = ones(1, 10);

init = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; log10(90.5); 0; log10(30);];
iter = 10;

reinfections("default_params_full.mat", params, options, model_settings, init, iter, re_specs)
