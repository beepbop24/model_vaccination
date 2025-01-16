%% CONTOUR PLOTS FOR IIS DISEASE-FREE STEADY STATES
% contour plots for disease-free steady states (IIS)
clearvars;
params = struct();
[X, Y] = meshgrid(linspace(-1.2, 1.5, 300), linspace(-3, 3, 300));


% contour plot with default parameters for IIS model (UU)
eq_contour_plot('default_params_full.mat', params, 'DF', 'I zero', X, Y, 'uu_zero');
eq_contour_plot('default_params_full.mat', params, 'DF', 'I star', X, Y, 'uu_star');

% contour plot with default parameters for IIS model (US)
params.b_2 = 4;
eq_contour_plot('default_params_full.mat', params, 'DF', 'I zero', X, Y, 'us_zero');
eq_contour_plot('default_params_full.mat', params, 'DF', 'I star', X, Y, 'us_star');

% contour plot with default parameters for IIS model (SU)
params.beta = 2.05*10^(-4);
params.b_2 = 0.15;
eq_contour_plot('default_params_full.mat', params, 'DF', 'I zero', X, Y, 'su_zero');
eq_contour_plot('default_params_full.mat', params, 'DF', 'I star', X, Y, 'su_star');


%% contour plot with default parameters for IIS model (SS) -- does not exist
params.beta = 2.05*10^(-12);
params.b_2 = 0.6;
%eq_contour_plot('default_params_full.mat', params, 'I zero', X, Y, 'ss_zero', 'DF');
%eq_contour_plot('default_params_full.mat', params, 'I star', X, Y, 'ss_star', 'DF');

X = linspace(-0.7, 0.5, 5000);
eq_realline('default_params_full.mat', params, 'I zero', X, 'test_izero')
eq_realline('default_params_full.mat', params, 'I star', X, 'test_istar')


%% CONTOUR PLOTS FOR VXY DISEASE-FREE STEADY STATES
clearvars;
default_params_full;
params = struct();
params.k_i = 0;
params.k_ix = 0;
params.k_iy = 0;
params.b_2 = 0;

[X, Y] = meshgrid(linspace(-1.5, 1.5, 300), linspace(-3, 3, 300));

params.beta = 0.05*defaults.d_v*defaults.d_y*defaults.d_x/(defaults.k * defaults.mu);
eq_contour_plot('default_params_full.mat', params, 'DF', 'I zero', X, Y, 'viral_contour_df');

params.beta = 10*defaults.d_v*defaults.d_y*defaults.d_x/(defaults.k * defaults.mu);
eq_contour_plot('default_params_full.mat', params, 'E', '', X, Y, 'viral_contour_e');

%% DISEASE-FREE EQUILIBRIUM CONTOUR PLOT
% function that produces contour plot of characteristic polynomial about
% linearized steady state for specific parameters (params), steady state value I^*,
% mesh (X,Y) and name of file that will be saved (title), and state (DF or E).
% If no parameters are provided, function will run with default values found in params_file

function eq_contour_plot(params_file, params, state, steady_state_I, X, Y, title)
    % loading defaults
    params_struct = load(params_file);
    defaults = params_struct.defaults;

    % setting param values with defaults if not specified in user input
    default_names = fieldnames(defaults);
    param_names = fieldnames(params);
    missing = find(~ismember(default_names, param_names));

    for k = 1:length(missing)
        params.(default_names{missing(k)}) = defaults.(default_names{missing(k)});
    end

    if strcmp(state, 'DF')
        % for which DF steady state are we plotting the contour plot
        if strcmp(steady_state_I, 'I zero')
            I_star = 0;
        elseif strcmp(steady_state_I, 'I star')
            I_star = params.b_2/params.d_i - params.k_2;
        else
            error('Not a valid steady state definition');
        end
    end
    
    % contour plot
    h = figure();
    if strcmp(state, 'DF')
        charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, I_star, params), X+1i*Y);
    elseif strcmp(state, 'E')
        charpoly_Df2 = arrayfun(@(lambda) characEqE(lambda, params), X+1i*Y);
    else
        error('Not a valid steady state definition');
    end
    contour(X, Y, charpoly_Df2, [0.00001 0.0001 0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10 25 50 100], 'LineWidth', 1.5)
    %hold on
    %plot(-params.d_r, 0, '.')
    %hold on
    %plot(-params.d_x, 0, '.')
    %hold on
    %plot(-params.d_i, 0, '.')
    hold on
    xline(0, 'LineWidth', 1.5)
    hold off
    ax = gca;
    ax.FontSize = 16;
    xlabel('$\Re(\lambda)$','Interpreter','latex', 'FontSize', 24)
    ylabel('$\Im(\lambda)$','Interpreter','latex', 'FontSize', 24)
    saveas(h, fullfile('./contour_plots', title), 'png')
end


%% FUNCTION THAT PLOTS CHARACTERISTIC POLYNOMIAL ALONG THE REAL LINE
% function that produces plot of Re(lambda) of characteristic polynomial
% about linearized steady state for specific parameters (params), steady
% state, range (X) and name of file that will be saved (title). if no
% parameters are provided, function will run with default values found in
% params_file
function eq_realline(params_file, params, steady_state_I, X, title)
    % loading defaults
    params_struct = load(params_file);
    defaults = params_struct.defaults;

    % setting param values with defaults if not specified in user input
    default_names = fieldnames(defaults);
    param_names = fieldnames(params);
    missing = find(~ismember(default_names, param_names));

    for k = 1:length(missing)
        params.(default_names{missing(k)}) = defaults.(default_names{missing(k)});
    end

    % for which DF steady state are we plotting the contour plot
    if strcmp(steady_state_I, 'I zero')
        I_star = 0;
    elseif strcmp(steady_state_I, 'I star')
        I_star = params.b_2/params.d_i - params.k_2;
    else
        error('Not a valid steady state definition');
    end

    % plot along real axis
    h = figure();
    charpoly_real = arrayfun(@(lambda) characEq(lambda, I_star, params), X);
    plot(X, charpoly_real, 'LineWidth', 1.5)
    ax = gca;
    ax.FontSize = 16;
    xlabel('$\lambda$','Interpreter','latex', 'FontSize', 24)
    ylabel('$\Delta(\lambda)$','Interpreter','latex', 'FontSize', 24)
    saveas(h, fullfile('./contour_plots', title), 'jpeg')
end


%% FUNCTION THAT COMPUTES ABSOLUTE VALUE OF CHARACTERISTIC EQUATION
% functions that computes the characteristic polynomial about linearized
% steady states for specific parameters (params), steady state, and lambda

% disease-free steady states
function charpoly_Df = characEq(lambda, I_star, params)
    if I_star == 0
        Df_x0 = [-params.d_v, 0, params.k*exp(-lambda*params.tau_1), 0, 0;
                 -params.beta*params.mu/params.d_x, -params.d_x, 0, 0, -params.k_ix*params.mu/params.d_x;
                 params.beta*params.mu/params.d_x, 0, -params.d_y, 0, 0;
                 0, 0, 0, -params.d_r, params.k_ix*params.mu/params.d_x;
                 0, 0, params.k_i*exp(-lambda*params.tau_2), 0, params.b_2/params.k_2*exp(-lambda*params.tau_4)-params.d_i];
    else
        Df_x0 = [-params.d_v, 0, params.k*params.k1_tilde^params.n_1/(params.k1_tilde^params.n_1+I_star^params.n_1)*exp(-lambda*params.tau_1), 0, 0;
                 -params.beta*params.mu/(params.d_x+params.k_ix*I_star), - params.d_x-params.k_ix*I_star, 0, 0, -params.k_ix*params.mu/(params.d_x+params.k_ix*I_star);
                 params.beta*params.mu/(params.d_x+params.k_ix*I_star), 0, -params.d_y - params.k_iy*I_star, 0, 0;
                 0, params.k_ix*I_star, params.k_iy*I_star, -params.d_r, params.k_ix*params.mu/(params.d_x+params.k_ix*I_star);
                 0, 0, params.k_i*exp(-lambda*params.tau_2), 0, params.k_2*params.b_2/(params.k_2+I_star)^2*exp(-lambda*params.tau_4)-params.d_i];
    end
    charpoly_Df = abs(det(Df_x0-lambda*eye(5)));

     %charpoly_Df = abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
        %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
        %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));
end

% endemic disease steady state for viral model (no IFN)
function charpoly_E = characEqE(lambda, params)
    Df_x0 = [-params.d_v, 0, params.k*exp(-lambda*params.tau_1);
             -params.d_v*params.d_y/params.k, -params.beta*params.k*params.mu/(params.d_v*params.d_y), 0;
             params.d_v*params.d_y/params.k, params.beta*params.k*params.mu/(params.d_v*params.d_y)-params.d_x, -params.d_y;];
   
    %charpoly_E = abs((params.d_v+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
            %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
            %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));

    charpoly_E = abs(det(Df_x0-lambda*eye(3)));
end

