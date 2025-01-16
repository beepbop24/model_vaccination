function reinfections(params_file, params, dde_options, model_settings, init, iter, re_specs)
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
        previous = randi([0, 6]);
    elseif mod(re_specs.previous, 1) == 0 && re_specs.previous < 6 && re_specs.previous >= 0
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

    % arrays for max (xvals)
    V_xvals_plot = zeros(iter, 1);
    X_xvals_plot = zeros(iter, 1);
    Y_xvals_plot = zeros(iter, 1);
    R_xvals_plot = zeros(iter, 1);
    I_xvals_plot = zeros(iter, 1);
    TH_xvals_plot = zeros(iter, 1);
    TE_xvals_plot = zeros(iter, 1);
    TM_xvals_plot = zeros(iter, 1);
    BLL_xvals_plot = zeros(iter, 1);
    BE_xvals_plot = zeros(iter, 1);
    A_xvals_plot = zeros(iter, 1);

    % arrays for max (yvals)
    V_yvals_plot = zeros(iter, 1);
    X_yvals_plot = zeros(iter, 1);
    Y_yvals_plot = zeros(iter, 1);
    R_yvals_plot = zeros(iter, 1);
    I_yvals_plot = zeros(iter, 1);
    TH_yvals_plot = zeros(iter, 1);
    TE_yvals_plot = zeros(iter, 1);
    TM_yvals_plot = zeros(iter, 1);
    BLL_yvals_plot = zeros(iter, 1);
    BE_yvals_plot = zeros(iter, 1);
    A_yvals_plot = zeros(iter, 1);

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
            previous = previous + 1;
            cons_epi = 0;
            tm_init = init(8);

            switch(previous)
                case 1
                    ax = -0.0225;
                case 2
                    ax = -0.03;
                case 3
                    ax = -0.034;
                case 4
                    ax = -0.045;
                case 5
                    ax = -0.0525;
                case 6
                    ax = -0.06;
                case 7
                    ax = -0.0675;
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
        init(8) = log10(cons_epi)+init(8); %memory T cells
        init(11) = init(11)+ax*ifi; %antibody
    
        % calculated from antibody
        init(9) = init(11)*params.d_a/params.k_blla; %long-lived b cells

        [full_soli, I_fulli, M_fulli] = full_model(params_file, params, @ddefullhist, dde_options, model_settings);
        xvals = linspace(full_soli.x(1), full_soli.x(end), 1000);
        yvals = deval(full_soli, xvals);

        I_fulli =  xvals(I_fulli);
     
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
         TM_plot(i,:) = yvals(8,:);
         BE_plot(i,:) = yvals(10,:);
        
        else
            init(8) = tm_init;
        end
        
        BLL_plot(i,:) = yvals(9,:);
        A_plot(i,:) = yvals(11,:);

          if type
            % max -- xvals
            V_xvals_plot(i) = I_fulli(1)-72.0720720720721;
            X_xvals_plot(i) = I_fulli(2);
            Y_xvals_plot(i) = I_fulli(3)-66.0660660660661;
            R_xvals_plot(i) = I_fulli(4)-76.0760760760761;
            I_xvals_plot(i) = I_fulli(5)-76.0760760760761;
            TH_xvals_plot(i) = I_fulli(6)-152.152152152152;
            TE_xvals_plot(i) = I_fulli(7)-228.228228228228;
            TM_xvals_plot(i) = I_fulli(8)-1171.17117117117;
            BE_xvals_plot(i) = I_fulli(10)-174.174174174174;
            
            % max -- yvals
            V_yvals_plot(i) = M_fulli(1)-5.97800884303669;
            X_yvals_plot(i) = M_fulli(2)-9.72015930340596;
            Y_yvals_plot(i) = M_fulli(3)-9.43598703288781;
            R_yvals_plot(i) = M_fulli(4)-1.25376219956616;
            I_yvals_plot(i) = M_fulli(5)-1.41419787061428;
            TH_yvals_plot(i) = M_fulli(6)-3.71458365203661;
            TE_yvals_plot(i) = M_fulli(7)-4.71652536931099;
            TM_yvals_plot(i) = M_fulli(8)-4.18255795817746;
            BE_yvals_plot(i) = M_fulli(10)-2.96339000918377;
            

        else
            TM_yvals_plot(i) = tm_init-4.18255795817746;
        end

        
        BLL_xvals_plot(i) = I_fulli(9)-166.166166166166;
        A_xvals_plot(i) = I_fulli(11)-540.540540540541;

        BLL_yvals_plot(i) = M_fulli(9)-3.55169736695018;
        A_yvals_plot(i) = M_fulli(11)-2.8033450921344;

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

    % plot for xvals (timing shift)
    x = linspace(1, iter, iter);
    figure();
    a = gca;
    a.FontSize = 16;
    colororder(["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"])
    isNZ = (~V_xvals_plot == 0);      
    scatter(x(isNZ), V_xvals_plot(isNZ), "filled");
    hold on
    isNZ = (~X_xvals_plot == 0);      
    scatter(x(isNZ), X_xvals_plot(isNZ), "filled");
    hold on
    isNZ = (~Y_xvals_plot == 0);      
    scatter(x(isNZ), Y_xvals_plot(isNZ), "filled");
    hold on
    isNZ = (~R_xvals_plot == 0);      
    scatter(x(isNZ), R_xvals_plot(isNZ), "filled");
    hold on
    isNZ = (~I_xvals_plot == 0);      
    scatter(x(isNZ), I_xvals_plot(isNZ), "filled");
    hold off
    legend('$V$', '$X$', '$Y$', '$R$', '$I$', 'Interpreter', 'latex')
    xlabel('Iteration', 'FontSize', 18);
    ylabel('Time (h)', 'FontSize', 18);

    figure();
    a = gca;
    a.FontSize = 16;
    colororder(["#CE7E00" "#C90076" "#6A329F" "#4DBEEE" "#A2142F" "#000000"])
    isNZ = (~TH_xvals_plot == 0);      
    scatter(x(isNZ), TH_xvals_plot(isNZ), "filled");
    hold on
    isNZ = (~TE_xvals_plot == 0);      
    scatter(x(isNZ), TE_xvals_plot(isNZ), "filled");
    hold on
    scatter(x, TM_xvals_plot, "filled");
    hold on
    scatter(x, BLL_xvals_plot, "filled");
    hold on
    isNZ = (~BE_xvals_plot == 0);      
    scatter(x(isNZ), BE_xvals_plot(isNZ), "filled");
    hold on
    scatter(x, A_xvals_plot, "filled");
    hold off
    legend('$T_H$', '$T_E$', '$T_M$', '$B_{LL}$', '$B_E$', '$A$', 'Interpreter', 'latex')
    xlabel('Iteration', 'FontSize',18);
    ylabel('Time (h)', 'FontSize',18);

    % plot for yvals (peak shift)
    x = linspace(1, iter, iter);
    figure();
    a = gca;
    a.FontSize = 16;
    colororder(["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"])
    isNZ = (~V_yvals_plot == 0);      
    scatter(x(isNZ), V_yvals_plot(isNZ), "filled");
    hold on
    isNZ = (~X_yvals_plot == 0);      
    scatter(x(isNZ), X_yvals_plot(isNZ), "filled");
    hold on
    isNZ = (~Y_yvals_plot == 0);      
    scatter(x(isNZ), Y_yvals_plot(isNZ), "filled");
    hold on
    isNZ = (~R_yvals_plot == 0);      
    scatter(x(isNZ), R_yvals_plot(isNZ), "filled");
    hold on
    isNZ = (~I_yvals_plot == 0);      
    scatter(x(isNZ), I_yvals_plot(isNZ), "filled");
    hold off
    legend('$V$', '$X$', '$Y$', '$R$', '$I$', 'Interpreter', 'latex')
    xlabel('Iteration', 'FontSize', 18);
    ylabel('log10 Fold Change', 'FontSize', 18);
    
    figure();
    a = gca;
    a.FontSize = 16;
    colororder(["#CE7E00" "#C90076" "#6A329F" "#4DBEEE" "#A2142F" "#000000"])
    isNZ = (~TH_yvals_plot == 0);      
    scatter(x(isNZ), TH_yvals_plot(isNZ), "filled");
    hold on
    isNZ = (~TE_yvals_plot == 0);      
    scatter(x(isNZ), TE_yvals_plot(isNZ), "filled");
    hold on
    scatter(x, TM_yvals_plot, "filled");
    hold on
    scatter(x, BLL_yvals_plot, "filled");
    hold on
    isNZ = (~BE_yvals_plot == 0);      
    scatter(x(isNZ), BE_yvals_plot(isNZ), "filled");
    hold on
    scatter(x, A_yvals_plot, "filled");
    hold off
    legend('$T_H$', '$T_E$', '$T_M$', '$B_{LL}$', '$B_E$', '$A$', 'Interpreter', 'latex')
    xlabel('Iteration', 'FontSize', 18);
    ylabel('log10 Fold Change', 'FontSize', 18);

    
    function plot_var(xvals, yvals)
        j = 1:iter;
        figure();
        a = gca;
        a.FontSize = 16;
        semilogy(xvals, 10.^(yvals), 'LineWidth', 1.5)
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

