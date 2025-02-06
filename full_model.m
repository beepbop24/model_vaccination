function [full_sol, I_full, M_full] = full_model(params_file, params, hist, dde_options, settings, title)
    % function that returns numerical simulations of full immune system
    % with INPUTS -- 
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % hist = history function of the DDE (function of t)
    % dde_options = ddeset options
    % settings.stiff: boolean -- if true use stiff solver
    % settings.manual: boolean -- if manual true, when number of infected cells <
    % tolerance, it is set to 0 (to ensure infection doesn't artifically
    % blow up
    % settings.tolerance: value at which infection is considered to be eliminated 
    % settings.plot: boolean -- if true, function produces plot of output
    % settings can be empty -- if so will run with defaults

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

    % setting settings values with default if not specified in user input
    default_settings.stiff = true;
    default_settings.manual = true;
    default_settings.tolerance = 5*10^(-5);
    default_settings.plot = true;
    default_settings.dashed = false;
    
    default_setting_names = fieldnames(default_settings);
    setting_names = fieldnames(settings);
    missing = find(~ismember(default_setting_names, setting_names));

    for k = 1:length(missing)
        settings.(default_setting_names{missing(k)}) = default_settings.(default_setting_names{missing(k)});
    end

    % setting lags for dde solver
    lags = [params.tau_1, params.tau_2, params.tau_35, params.tau_4, ...
            params.tau_6, params.tau_7, params.tau_8, params.tau_9, ...
            params.tau_10, params.tau_11];

    % using stiff DDE solver
    if settings.stiff
        stiff_start = tic;
        full_sol = dde15s_updated(@ddefull, lags, hist, params.tspan, dde_options);  
        stiff_end = toc(stiff_start);
        disp(['Stiff Time: ', num2str(stiff_end)])

    % using regular DDE solver
    else
        reg_start = tic;
        full_sol = dde23(@ddefull, lags, hist, params.tspan, dde_options);
        reg_end = toc(reg_start);
        disp(['Regular Time: ', num2str(reg_end)])
    end
    
    % evaluating output for 1000 evenly-spaced points
    xvals = linspace(full_sol.x(1), full_sol.x(end), 1000);
    yvals = deval(full_sol, xvals);

    % Peak Values
    [M_full, I_full] = max(yvals.');
  
    % Peak Values Display
    variables = {'V', 'X', 'Y', 'R', 'I', 'T_H', 'T_E', 'T_M', 'B_{LL}', 'B_E', 'A'}';
    columns = {'Variable', 'Time of Peak', 'Peak Value'};
    values = cat(2, xvals(I_full).', M_full.');

    peak_summary = cell2table([variables, num2cell(values)], 'VariableNames', columns);
    peak_summary.(1) = categorical(peak_summary.(1));
    disp(peak_summary)

    % plot of model simulation
    if settings.plot
        h = figure();
        semilogy(xvals, 10.^(yvals), 'LineWidth', 1.5)

        if settings.dashed
            hold on
            yline(5.25*10^9, '--', 'LineWidth', 1.5, 'Color', 'black')
            xval = find(yvals(3, :) >= log10(5.25*10^9), 1);
            xline(xvals(xval), '--', 'LineWidth', 1.5, 'Color', 'black')
        end
        
        hold off
        ax = gca;
        ax.FontSize = 16;
        colororder(["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30" "#CE7E00" "#C90076" "#6A329F" "#4DBEEE" "#A2142F" "#000000"])
        xlabel('Time (h)', 'FontSize',18);
        %ylabel('Number of Cells', 'FontSize',18);
        legend('$V$', '$X$', '$Y$', '$R$', '$I$', '$B_{LL}$', '$B_E$', '$A$', '$T_E$', '$T_M$', '$T_H$', 'Interpreter', 'latex')
        saveas(h, fullfile('./simulations', title), 'png')
    end
    
% -------------------------------------------------------------------------

    function dydt = ddefull(t, y, Z)
        ylag1 = Z(:,1); %tau_1
        ylag2 = Z(:,2); %tau_2
        ylag3 = Z(:,3); %tau_35
        ylag4 = Z(:,4); %tau_4
        ylag6 = Z(:,5); %tau_6
        ylag7 = Z(:,6); %tau_7
        ylag8 = Z(:,7); %tau_8
        ylag9 = Z(:,8); %tau_9
        ylag10 = Z(:,9); %tau_10
        ylag11 = Z(:,10); %tau_11

        % default computation for Y 
        y_infected = params.beta*y(2)*y(1) - params.d_y*y(3) - params.k_iy*y(3)*y(5)-params.k_tey*y(3)*y(7);
        %disp(y_infected)

        % if variables smaller than tolerance force it to 0 to prevent
        % "reinfections" from occuring
        if settings.manual == true && y(3) < settings.tolerance
            y_infected = -y(3);
        end
        %% threshold for T and B cells
        % default values below threshold
        t_h = -params.d_th*y(6);
        t_e = -params.d_te*y(7);
        t_m = 0;

        b_ll = 0;
        b_e = -params.d_be*y(10);

        % threshold conditions -- T cells
        if y(1) > params.v_star
            t_h = t_h + params.k_th*y(5)*y(3);
            t_e = t_e + params.k_tetm*y(8)*y(3);
            t_m = t_m - params.k_tetm*y(8)*y(3);
        end

        if ylag11(1) > params.v_star
            %P = log10(params.t_total-10^(y(10)));
            P = log10(params.t_total^10-ylag11(8)^10);
            t_e = t_e + params.k_te*ylag11(3)*(1-ylag11(7)/P)*ylag11(5)*ylag11(6); 
        end

        % generation of memory T cells post-infection
        if y(3) < settings.tolerance
            t_e = t_e - params.d_tepi*y(7) - params.d_tmpi*y(7);
            t_m = params.d_tmpi*y(7);
  
            if settings.manual == true && y(6) < settings.tolerance
                t_h = -y(6);
            end

            if settings.manual == true && y(7) < settings.tolerance
                t_e = -y(7);
            end
        end
        
        % threshold conditions -- B cells
        if ylag8(1) > params.v_star
            b_ll = b_ll + params.k_bll*ylag8(1)*ylag8(6);
        end

        if ylag9(1) > params.v_star
            b_ll = b_ll - params.k_bllbe*ylag9(1)*ylag9(9);
            b_e = b_e + params.k_bllbe*ylag9(1)*ylag9(9);
        end

        if ylag10(1) > params.v_star
            b_e = b_e + params.k_be*ylag10(1)*ylag10(6);
        end

        dydt = [params.k*ylag1(3)*params.k1_tilde^params.n_1/(params.k1_tilde^params.n_1 + ylag3(5)^params.n_1) - params.d_v*y(1) - params.rhoV*y(1)*y(11)
                params.mu - params.d_x*y(2) - params.beta*y(2)*y(1) - params.k_ix*y(2)*y(5)
                y_infected
                params.k_ix*y(2)*y(5) + params.k_iy*y(3)*y(5) - params.d_r*y(4)
                params.k_i * ylag2(3) + (params.b_2*ylag4(5)^params.n_2)/(params.k_2^params.n_2 + ylag4(5)^params.n_2) + params.k_ite*ylag6(7) + params.k_ith*ylag7(6) - params.d_i*y(5)
                t_h
                t_e
                t_m
                b_ll
                b_e
                params.k_blla*y(9) + params.k_bea*y(10) - params.rhoA*y(1)*y(11)- params.d_a*y(11)];
    end
end