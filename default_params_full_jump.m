% parameter values for V
defaults.k = 0.84; 
defaults.k1_tilde = 3; 
defaults.n_1 = 1.775; 
defaults.d_v = 1; 
defaults.rhoV = 0.03;

% parameter values for X
defaults.d_x = 2.725; 
defaults.mu = log10(5.25*10^9)*defaults.d_x; 
defaults.beta = 2.11*10^(-2);
defaults.k_ix = 2.5*10^(-1);

% parameter values for Y
defaults.d_y = 1/18; %OK
defaults.k_iy = 3.75*10^(-2);
defaults.k_tey = 0.055;

% parameter values for R
defaults.d_r = defaults.d_x; %OK

% parameter values for I
defaults.k_i = 0.091;
defaults.b_2 = 0.15;
defaults.k_2 = 0.8;
defaults.n_2 = 1;
defaults.d_i = 0.675; %OK
defaults.k_ite = 0.05;
defaults.k_ith = 0.003;

% threshold
defaults.v_star = 4; %O 

% parameter values for T_H
defaults.k_th = 0.003935; %OK
defaults.d_th = 0.006; %OK

% parameter values for T_E
defaults.k_te = 0.00155;
defaults.t_total = log10(5.2*10^4); %OK
defaults.k_tetm = 6.7699494710*10^(-5);
defaults.d_tepi = 0.90*(log10(52000)/336); %OK
defaults.d_te = 0.002/24; %OK
defaults.d_tmpi = 0.1*(log10(52000)/336); %OK

% parameter values for B_LL
defaults.k_bll = 0.001776696245; %OK 
defaults.k_bllbe = 0.00075;

% parameter values for B_E
defaults.k_be = 0.0023;
defaults.d_be = 0.01; %OK

% parameter values for A
defaults.d_a = 0.04/24; %OK
defaults.k_blla = log10(30)*defaults.d_a/log10(90.5); %OK
defaults.k_bea = defaults.k_blla*2.08; %OK
defaults.rhoA = 0.000005;

% lag values
defaults.tau_1 = 8; %OK
defaults.tau_2 = 8; %OK
defaults.tau_35 = 5; 
defaults.tau_4 = 4; 
defaults.tau_6 = 5; %OK
defaults.tau_7 = 7;
defaults.tau_8 = 10;
defaults.tau_9 = 8;
defaults.tau_10 = 20;
defaults.tau_11 = 72; %OK

% duration of simulation
defaults.tspan = [0, 2000];

save('default_params_full_jump.mat','defaults','defaults');
