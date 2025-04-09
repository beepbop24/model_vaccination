function inset_plots(params_file, params, hist, dde_options, settings, steady_state_I, state, title)
    % function that returns numerical simulations of innate immune system
    % with inset representing convergence to DF or E (for viral dynamics)
    % steady state 
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % hist:  history function of the DDE (function of t) for the first
    % infection
    % dde_options = ddeset options
    % model_settings.stiff: boolean -- if true use stiff solver
    % model_settings.manual: boolean -- if manual true, when number of infected
    % cells < tolerance, it is set to 0 (to ensure infection doesn't
    % artifically blow up)
    % model_settings.tolerance: value at which infection is considered to 
    % be eliminated 
    % model_settings can be empty -- if so will run with defaults
    % steady_state_I: 'I zero' or 'I star'
    % state: 'iis_DF', 'viral_DF' or 'viral_E'
    % title: name of png file that will be saved.

    settings.plot = false;
    IIS_sol = IIS_model(params_file, params, hist, dde_options, settings);

    % evaluating output for 1000 evenly-spaced points
    xvals = linspace(IIS_sol.x(1), IIS_sol.x(end), 1000);
    yvals = deval(IIS_sol, xvals);

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
   

    if strcmp(state, 'iis_DF')
        if strcmp(steady_state_I, 'I zero')
            i_star = 0;
        elseif strcmp(steady_state_I, 'I star')
            i_star = params.b_2/params.d_i - params.k_2;
        else
            error('Not a valid steady state definition');
        end
        
        x_star = params.mu/(params.d_x+params.k_ix*i_star);
        r_star = params.mu*params.k_ix*i_star/(params.d_r*(params.d_x+params.k_ix*i_star));
        star = [0; x_star; 0; r_star; i_star;];


    % ony for no immune system dynamics    
    elseif strcmp(state,'viral_DF')
         star = [0; params.mu/params.d_x; 0; 0; 0;];

    elseif strcmp(state, 'E')
        v_star = params.k*params.mu/(params.d_v*params.d_y)-params.d_x/params.beta;
        x_star = params.d_v*params.d_y/(params.k*params.beta);
        y_star = params.mu/params.d_y-params.d_v*params.d_x/(params.k*params.beta);

        star = [v_star; x_star; y_star; 0; 0;];

    else
        error('Not a valid steady state definition');
    end

    FontSize = 18; Width = 1.5;

    h = figure();
    if strcmp(state, 'iis_DF')
    semilogy(xvals, 10.^(yvals), 'Linewidth', Width)
    else
        semilogy(xvals, 10.^(yvals(1:3,:)), 'Linewidth', Width)
    end
    set(gca, 'FontSize', FontSize);
    colororder(["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"])
    xlabel('Time (h)');
    ylabel('Magnitude of Infection and IR');

    if strcmp(state, 'iis_DF')
        legend('$V$', '$X$', '$Y$', '$R$', '$I$', 'Interpreter', 'latex', 'Position', [0.775 0.725 0.1 0.1])
    else
        legend('$V$', '$X$', '$Y$', 'Interpreter', 'latex', 'Position', [0.775 0.725 0.1 0.1])
    end

    % inset
    axes('Position', [0.6 0.275 0.25 0.3]) % numbers are adjustable: https://www.mathworks.com/help/matlab/ref/axes.html
    box on; 
    semilogy(xvals, vecnorm(yvals-star), 'Linewidth', Width, 'Color', '#A2142F'); % difference to steady state y*
    legend('$| x-x^* |$', 'Interpreter', 'latex')
    set(gca,'FontSize', FontSize); 
    saveas(h, fullfile('./insets', title), 'png')
end