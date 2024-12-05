function reinfection_sol = reinfections(params_file, params, dde_options, model_settings, init, iter, re_specs)
    % function that returns numerical simulations of the full immune system
    % for reinfections with INPUTS -- 
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % dde_options = ddeset options
    % model_settings.stiff: boolean -- if true use stiff solver
    % model_settings.manual: boolean -- if manual true, when number of infected cells <
    % tolerance, it is set to 0 (to ensure infection doesn't artifically blow up
    % model_settings.tolerance: value at which infection is considered to be eliminated 
    % model_settings.plot: boolean -- if true, function produces plot of output
    % model_settings can be empty -- if so will run with defaults
    % init:  history function of the DDE (function of t) for the first
    % infection
    % iter: number of reinfections (total number of infections includes
    % original infection so it's iter + 1)
    % plot_type: "all" (each infection = 1 plot with all variables) or 
    % "over" (each variable = 1 plot with all infections + reinfections)
    % time: time between reinfections (immunity decreases over long periods
    % of time)
    % re_specs.cons_epi: float percentage (e.g. 0.7) of conserved epitopes
    % between infections or "rand"
    % re_specs.ifi: inter flu interval: integer number of months between
    % infections or "rand"
    % re_specs.previous: integer number of previous infections (for vax --
    % set to 0 for natural infections)
    % re_specs.season: 0 (winter) or 1 (summer) or "rand"

    % setting setting values with default if not specified in user input
    default_settings.stiff = true;
    default_settings.manual = true;
    default_settings.tolerance = 5*10^(-5);
    default_settings.plot = false;
    
    default_setting_names = fieldnames(default_settings);
    setting_names = fieldnames(model_settings);
    missing = find(~ismember(default_setting_names, setting_names));

    for k = 1:length(missing)
        model_settings.(default_setting_names{missing(k)}) = default_settings.(default_setting_names{missing(k)});
    end

    % loading default param values
    params_struct = load(params_file);
    defaults = params_struct.defaults;

    % setting param values with defaults if not specified in user input
    default_names = fieldnames(defaults);
    param_names = fieldnames(params);
    missing = find(~ismember(default_names, param_names));

    for k = 1:length(missing)
        params.(default_names{missing(k)}) = defaults.(default_names{missing(k)});
    end
    
    % setting reinfection specifications with default if not specified in user input
    default_specs.cons_epi = 100;
    default_specs.type = 100;
    default_specs.ifi = 100;
    default_specs.previous = 0;
    default_specs.season = 100;
    
    default_specs_names = fieldnames(default_specs);
    re_spec_names = fieldnames(re_specs);
    missing = find(~ismember(default_specs_names, re_spec_names));

    for k = 1:length(missing)
        re_specs.(default_specs_names{missing(k)}) = default_specs.(default_specs_names{missing(k)});
    end

    % setting up number of previous infections
    if isscalar(re_specs.cons_epi) && re_specs.previous == 100
        previous = randi([1, 7]);
    elseif mod(re_specs.previous, 1) == 0 && re_specs.previous < 7 && re_specs.previous >= 0
        previous = re_specs.previous;
    else
        error("Not a valid number of previous infections")
    end

    % set up -- same number of infected cells for all subsequent infections
    y_infected = init(3);

    % arrays for plot values for all iters
    model_settings.plot = false;
    V_plot = zeros(iter,1000);
    X_plot = zeros(iter, 1000);
    Y_plot = zeros(iter,  1000);
    R_plot = zeros(iter,  1000);
    I_plot = zeros(iter,  1000);
    TH_plot = zeros(iter,  1000);
    TE_plot = zeros(iter,  1000);
    TM_plot = zeros(iter,  1000);
    BLL_plot = zeros(iter,  1000);
    BE_plot = zeros(iter,  1000);
    A_plot = zeros(iter,  1000);


    %% iter # of REINFECTIONS
    for i=1:iter

        % time since previous infection -- if not specified, follows a normal
        % distribution with mu = 36, sigma = 12; otherwise must be in
        % {1, .., 60}
        if isscalar(re_specs.ifi) && re_specs.ifi == 100
            ifi = normrnd(36, 12);

        elseif mod(re_specs.ifi(i), 1) == 0 && re_specs.ifi(i) < 60 && re_specs.ifi(i) >= 0
            ifi = re_specs.ifi(i);
        else
            error("Not a valid 'months since last infection' value")
        end

        % type of "infection" (50-50 natural infection or vax)
        if isscalar(re_specs.type) && re_specs.type == 100
            x = rand;
            type = (x > 0.5);

        elseif re_specs.type(i) == 0 || re_specs.type(i) == 1
            type = re_specs.type(i);

        else
            error("Not a valid type of infection")
        end

        % natural infection
        if type
            % conserved epitopes -- if not specified, follows a uniform
            % distribution [0.3, 1] (30% to 100%); otherwise must be in
            % [0.3, 1]
            if isscalar(re_specs.cons_epi) && re_specs.cons_epi == 100
                cons_epi = 0.7*rand+0.3;

            elseif re_specs.cons_epi(i) <= 1 && re_specs.cons_epi(i) >= 0.3
                cons_epi = re_specs.cons_epi(i);
            else
                error("Not a valid percentage of conserved epitopes")
            end

            % "similar" infection
            if cons_epi > 0.7
                ax = -30*0.5/1095;
            else
                ax = -30*0.00125;
            end
        
        % vax    
        else
            cons_epi = 0;

            switch(previous)
                case 0
                    ax = -0.025;
                case 1
                    ax = -0.025;
                case 2
                    ax = -0.03;
                case 3
                    ax = -0.0371791;
                case 4
                    ax = -0.0413043;
                case 5
                    ax = -0.0457608;
                case 6
                    ax = -0.0522826;
                case 7
                    ax = -0.0571739;
            end
        end

        % season -- if not specified, follows a distribution W (0.8) and S
        % (0.2)
        % ; otherwise must be 0 (S) or 1 (W) or 
        if isscalar(re_specs.season) && re_specs.season == 100
            x = rand;
            season = (x > 0.2);

        elseif re_specs.season(i) == 0 || re_specs.season(i) == 1
            season = re_specs.season(i);
        else
            error("Not a valid season")
        end
    
        % if winter, IFN inhibition of virus less effective
        if ~season 
            params.n_1 = 0.95*1.775;
        else
            params.n_1 = 1.775;
        end
    
        % initial number of cells
        init(2) = log10(5.25*10^9); %target cells
        init(3) = y_infected; %infected cells
        init(8) = cons_epi*init(8); %memory T cells
        init(11) = init(11)+ax*ifi; %antibody
    
        % calculated from antibody
        init(9) = init(11)*params.d_a/params.k_blla; %long-lived b cells

        full_soli = full_model(params_file, params, @ddefullhist, dde_options, model_settings);
        xvals = linspace(full_soli.x(1), full_soli.x(end), 1000);
        yvals = deval(full_soli, xvals);
     
        % init for next infection
        init = yvals(:, end);
        
        % only plotting natural infection
        if type
         V_plot(i,:) = yvals(1,:);
         X_plot(i,:) = yvals(2,:);
         Y_plot(i,:) = yvals(3,:);
         R_plot(i,:) = yvals(4,:);
         I_plot(i,:) = yvals(5,:);
         TH_plot(i,:) = yvals(6,:);
         TE_plot(i,:) = yvals(7,:);
         BE_plot(i,:) = yvals(10,:);
        
        else
            previous = previous + 1;
        end

        TM_plot(i,:) = yvals(8,:);
        BLL_plot(i,:) = yvals(9,:);
        A_plot(i,:) = yvals(11,:);

        % summary of specs display
        variables = {'Conserved Epitopes', 'Type of Infection', 'Number of Months Since Previous', 'Number of Previous', 'Season', 'Iteration'}';
        columns = {'Spec', 'Value'};
        values = [cons_epi; type; ifi; previous; season; i];
        summary = cell2table([variables, num2cell(values)], 'VariableNames', columns);
        summary.(1) = categorical(summary.(1));
        disp(summary)
    end

    % plot for each variable
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
    
    function plot_var(xvals, yvals)
        j = 0:iter;
        figure();
        a = gca;
        a.FontSize = 16;
        semilogy(xvals, 10.^(yvals))
        xlabel('Time (h)', 'FontSize',18);
        ylabel('Number of Cells', 'FontSize',18);
        legend (sprintfc('G%i', j));
    end
    
    function s = ddefullhist(t)
        % making sure initial conditions which should be 0 are actually 0
        if model_settings.manual
            init(init < model_settings.tolerance) = 0;
        end
        
        % constant history function for full model.
        s = init;
    end
end

