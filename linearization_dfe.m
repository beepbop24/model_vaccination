%% CONTOUR PLOTS FOR IIS DISEASE-FREE STEADY STATES
% contour plot -- overview
[X, Y] = meshgrid(linspace(-1.5, 1.5, 300), linspace(-3, 3, 300));
params = struct();

% contour plot with default parameters for IIS model (UU)
eq_contour_plot('default_params_full.mat', params, 'I zero', X, Y, 'uu_zero');
eq_contour_plot('default_params_full.mat', params, 'I star', X, Y, 'uu_star');

% contour plot with default parameters for IIS model (US)
params.b_2 = 9;
eq_contour_plot('default_params_full.mat', params, 'I zero', X, Y, 'us_zero');
eq_contour_plot('default_params_full.mat', params, 'I star', X, Y, 'us_star');

% contour plot with default parameters for IIS model (SU)
params.beta = 2.05*10^(-4);
params.b_2 = 0.2;
eq_contour_plot('default_params_full.mat', params, 'I zero', X, Y, 'su_zero');
eq_contour_plot('default_params_full.mat', params, 'I star', X, Y, 'su_star');

% contour plot with default parameters for IIS model (SS)
params.beta = 1.5*10^(-3);
params.b_2 = 0.9;
eq_contour_plot('default_params_full.mat', params, 'I zero', X, Y, 'ss_zero');
eq_contour_plot('default_params_full.mat', params, 'I star', X, Y, 'ss_star');


%% CONTOUR PLOTS FOR VXY DISEASE-FREE STEADY STATES
params = struct();
[X, Y] = meshgrid(linspace(-1.5, 1.5, 300), linspace(-3, 3, 300));
eq_contour_plot('default_params_VXY.mat', params, 'I zero', X, Y, 's_contour_VXY');

X = linspace(-0.75, 0.1, 100);
eq_realline('default_params_VXY.mat', params, 'I zero', X, 's_realline_VXY');


%% DISEASE-FREE EQUILIBRIUM CONTOUR PLOT
% function that produces contour plot of characteristic polynomial about
% linearized steady state for specific parameters (params), steady state,
% mesh (X,Y) and name of file that will be saved (title). if no parameters
% are provided, function will run with default values found in params_file

function eq_contour_plot(params_file, params, steady_state, X, Y, title)
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
    if strcmp(steady_state, 'I zero')
        I_star = 0;
    elseif strcmp(steady_state, 'I star')
        I_star = params.b_2/params.d_i - params.k_2;
    else
        error('Not a valid steady state definition');
    end
    
    % contour plot
    h = figure();
    charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, I_star, params), X+1i*Y);
    contour(X, Y, charpoly_Df2, [0.0001 0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10 25 50 100], 'LineWidth', 1.5)
    %hold on
    %plot(-params.d_r, 0, '.')
    %hold on
    %plot(-params.d_x, 0, '.')
    %hold on
    %plot(-params.d_i, 0, '.')
    hold on
    xline(0, 'LineWidth', 1.5)
    hold off
    xlabel('$\Re(\lambda)$','Interpreter','latex', 'FontSize', 14)
    ylabel('$\Im(\lambda)$','Interpreter','latex', 'FontSize', 14)
    saveas(h, fullfile('./contour_plots', title), 'jpeg')

    %char_eq = @(lambda) abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
        %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
        %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));

    %test_value = fsolve(char_eq, 0.01+0.01*1i);
end

%% FUNCTION THAT PLOTS CHARACTERISTIC POLYNOMIAL ALONG THE REAL LINE
% function that produces plot of Re(lambda) of characteristic polynomial
% about linearized steady state for specific parameters (params), steady
% state, range (X) and name of file that will be saved (title). if no
% parameters are provided, function will run with default values found in
% params_file
function eq_realline(params_file, params, steady_state, X, title)
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
    if strcmp(steady_state, 'I zero')
        I_star = 0;
    elseif strcmp(steady_state, 'I star')
        I_star = params.b_2/params.d_i - params.k_2;
    else
        error('Not a valid steady state definition');
    end

    % plot along real axis
    h = figure();
    charpoly_real = arrayfun(@(lambda) characEq(lambda, I_star, params), X);
    plot(X, charpoly_real, 'LineWidth', 1.5)
    xlabel('$\lambda$','Interpreter','latex', 'FontSize', 14)
    ylabel('$\Delta(\lambda)$','Interpreter','latex', 'FontSize', 14)
    saveas(h, fullfile('./contour_plots', title), 'jpeg')
end


%% FUNCTION THAT COMPUTES ABSOLUTE VALUE OF CHARACTERISTIC EQUATION
% function that computes the characteristic polynomial about linearized
% steady state for specific parameters (params), steady state, and lambda
function charpoly_Df = characEq(lambda, I_star, params)
    if I_star == 0
        Df_x0 = [-params.d_v, 0, params.k/params.k1_tilde*exp(-lambda*params.tau_1), 0, 0;
                 -params.beta*params.mu/params.d_x, - params.d_x, 0, 0, -params.k_ix*params.mu/params.d_x;
                 params.beta*params.mu/params.d_x, 0, -params.d_y, 0, 0;
                 0, 0, 0, -params.d_r, params.k_ix*params.mu/params.d_x;
                 0, 0, params.k_i*exp(-lambda*params.tau_2), 0, -params.b_2/params.k_2*exp(-lambda-params.tau_4)-params.d_i];
            %charpoly_Df = abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
            %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
            %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));
    else
        Df_x0 = [-params.d_v, 0, params.k/(params.k1_tilde+I_star^params.n_1)*exp(-lambda*params.tau_1), 0, 0;
                 -params.beta*params.mu/(params.d_x+params.k_ix*I_star), - params.d_x-params.k_ix*I_star, 0, 0, -params.k_ix*params.mu/(params.d_x+params.k_ix*I_star);
                 params.beta*params.mu/(params.d_x+params.k_ix*I_star), 0, -params.d_y - params.k_iy*I_star, 0, 0;
                 0, params.k_ix*I_star, params.k_iy*I_star, -params.d_r, params.k_ix*params.mu/(params.d_x+params.k_ix*I_star);
                 0, 0, params.k_i*exp(-lambda*params.tau_2), 0, params.k_2*params.d_i/params.b_2*exp(-lambda*params.tau_4)-params.d_i];
    end

    charpoly_Df = abs(det(Df_x0-lambda*eye(5)));
end