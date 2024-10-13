function reinfection_sol = reinfections(params_file, params, dde_options, model_settings, init, iter, variable, time)
    % function that returns numerical simulations of the full immune system
    % for reinfections with INPUTS -- 
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % dde_options = ddeset options
    % model_settings.stiff: boolean -- if true use stiff solver
    % model_settings.manual: boolean -- if manual true, when number of infected cells <
    % tolerance, it is set to 0 (to ensure infection doesn't artifically
    % blow up
    % model_settings.tolerance: value at which infection is considered to be eliminated 
    % model_settings.plot: boolean -- if true, function produces plot of output
    % model_settings can be empty -- if so will run with defaults
    % init:  history function of the DDE (function of t) for the first
    % infection
    % iter: number of reinfections (total number of infections includes
    % original infection so it's iter + 1)
    % variable: "all" (each infection = 1 plot with all variables) or 
    % "1over" (each variable = 1 plot with all infections + reinfections)
    % time: time between reinfections (immunity decreases over long periods
    % of time)

    
    % setting setting values with default if not specified in user input
    default_settings.stiff = true;
    default_settings.manual = true;
    default_settings.tolerance = 5*10^(-5);
    default_settings.plot = true;
    
    default_setting_names = fieldnames(default_settings);
    setting_names = fieldnames(model_settings);
    missing = find(~ismember(default_setting_names, setting_names));

    for k = 1:length(missing)
        model_settings.(default_setting_names{missing(k)}) = default_settings.(default_setting_names{missing(k)});
    end
    
   % for which variables are we plotting the output
   % will result in separate plots for all the variables (11 plots)
    if strcmp(variable, '1over')
        model_settings.plot = false;
        V_plot = zeros(iter+1, 1000);
        X_plot = zeros(iter+1, 1000);
        Y_plot = zeros(iter+1, 1000);
        R_plot = zeros(iter+1, 1000);
        I_plot = zeros(iter+1, 1000);
        TH_plot = zeros(iter+1, 1000);
        TE_plot = zeros(iter+1, 1000);
        TM_plot = zeros(iter+1, 1000);
        BLL_plot = zeros(iter+1, 1000);
        BE_plot = zeros(iter+1, 1000);
        A_plot = zeros(iter+1, 1000);
    
    % will result in a single plot for all the variables
    elseif strcmp(variable, 'all')
        %
       
    elseif strcmp(variable, '1sep')
        %
    
    else
        error('Not a valid choice of variables');
    end

    % same number of infected cells for all subsequent infections
    y_infected = init(3);

    % G0 -- Initial Infection
    full_sol0 = full_model(params_file, params, @ddefullhist, dde_options, model_settings);
    init = full_sol0.y(:, end);
    init(3) = y_infected;

    if strcmp(variable, '1over')
        % evaluating output for 1000 evenly-spaced points
        xvals = linspace(full_sol0.x(1), full_sol0.x(end), 1000);
        yvals = deval(full_sol0, xvals);

        V_plot(1,:) = yvals(1,:);
        X_plot(1,:) = yvals(2,:);
        Y_plot(1,:) = yvals(3,:);
        R_plot(1,:) = yvals(4,:);
        I_plot(1,:) = yvals(5,:);
        TH_plot(1,:) = yvals(6,:);
        TE_plot(1,:) = yvals(7,:);
        TM_plot(1,:) = yvals(8,:);
        BLL_plot(1,:) = yvals(9,:);
        BE_plot(1,:) = yvals(10,:);
        A_plot(1,:) = yvals(11,:);
    end
    
    for i=1:iter
        disp([num2str(i),' infection'])
        full_soli = full_model(params_file, params, @ddefullhist, dde_options, model_settings);
        init = full_soli.y(:, end);
        init(3) = y_infected;

        if strcmp(variable, '1over')
            % evaluating output for 1000 evenly-spaced points
            xvals = linspace(full_soli.x(1), full_soli.x(end), 1000);
            yvals = deval(full_soli, xvals);

            V_plot(i+1,:) = yvals(1,:);
            X_plot(i+1,:) = yvals(2,:);
            Y_plot(i+1,:) = yvals(3,:);
            R_plot(i+1,:) = yvals(4,:);
            I_plot(i+1,:) = yvals(5,:);
            TH_plot(i+1,:) = yvals(6,:);
            TE_plot(i+1,:) = yvals(7,:);
            TM_plot(i+1,:) = yvals(8,:);
            BLL_plot(i+1,:) = yvals(9,:);
            BE_plot(i+1,:) = yvals(10,:);
            A_plot(i+1,:) = yvals(11,:);
        end

    end

    reinfection_sol = full_soli;

    if strcmp(variable, '1over')
        plot_var(xvals, V_plot)
        plot_var(xvals, X_plot)
        plot_var(xvals, Y_plot)
        plot_var(xvals, R_plot)
        plot_var(xvals, I_plot)
        plot_var(xvals, TH_plot)
        plot_var(xvals, TE_plot)
        plot_var(xvals, TM_plot)
        plot_var(xvals, BLL_plot)
        plot_var(xvals, BE_plot)
        plot_var(xvals, A_plot)  
    end

    function s = ddefullhist(t)
        % constant history function for full model.
        if model_settings.manual
            init(init<model_settings.tolerance) = 0;
        end
        
        s = init;
    end

    function plot_var(xvals, yvals)
        i = 0:iter;
        figure();
        semilogy(xvals, 10.^(yvals))
        xlabel('Time (h)');
        ylabel('Number of Cells');
        legend (sprintfc('G%i', i));
    end
end