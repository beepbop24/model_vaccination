% parameter values for V
params.k = 100; % 
params.k1_tilde = 50; %*
params.n_1 = 1; %*
params.d_v = 0.75;

% parameter values for X
params.mu = 7*10^6*log(2)/4320;
params.d_x = log(2)/4320;
params.beta = 0.0001;
params.k_ix = 0; % 0.000000000000005;

% parameter values for Y
params.d_y = 1/12;
params.k_iy = 0; %0.00005;

% parameter values for R
params.d_r = log(2)/4320;

% parameter values for I
params.k_i = 0;
params.b_2 = 0;
params.k_2 = 0;
params.n_2 = 0;
params.d_i = 0.7;

% lag values
lags.tau_1 = 8;
lags.tau_2 = 6;
lags.tau_35 = 11; % estimated
lags.tau_4 = 9; % estimated

% span
tspan = [0, 200];
save('param_values.mat','params','lags', 'tspan');