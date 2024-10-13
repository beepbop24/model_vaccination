% parameter values for V
defaults.k = 12.15; 
defaults.k1_tilde = 5; 
defaults.n_1 = 1.775; %OK
defaults.d_v = 0.9; %OK

% parameter values for X
defaults.d_x = 2.725; %OK
defaults.mu = log10(5.25*10^9)*defaults.d_x; %OK
defaults.k_ix = 2.5*10^(-1); 

% parameter values for Y
defaults.d_y = 1/18; %OK
defaults.k_iy = 4*10^(-2);

% parameter values for R
defaults.d_r = defaults.d_x; %OK

% parameter values for I
defaults.k_i = 0;
defaults.b_2 = 0.15;
defaults.k_2 = 0.8;
defaults.n_2 = 1;
defaults.d_i = 0.675; %OK

% lag values
defaults.tau_1 = 8; %OK
defaults.tau_2 = 6; %OK
defaults.tau_35 = 5; 
defaults.tau_4 = 6; 

% computed values
defaults.beta = defaults.d_x*defaults.d_v*defaults.d_y/(defaults.k*defaults.mu)-0.00008;
r0 = defaults.beta*defaults.k*defaults.mu/(defaults.d_x*defaults.d_v*defaults.d_y);
disp(defaults.beta)
disp(r0)

save('default_params_VXY.mat','defaults');