
function [IIS_sol, M_full] = IIS_model(params_file, params, hist, options, settings, title)
    % function that returns numerical simulations of innate immune system
    % with INPUTS -- 
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % hist = history function of the DDE (function of t)
    % dde_options = ddeset options
    % stiff: boolean -- if true use stiff solver

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

    default_settings.stiff = true;
    default_settings.plot = true;
    
    default_setting_names = fieldnames(default_settings);
    setting_names = fieldnames(settings);
    missing = find(~ismember(default_setting_names, setting_names));

    for k = 1:length(missing)
        settings.(default_setting_names{missing(k)}) = default_settings.(default_setting_names{missing(k)});
    end

    % setting lags for dde solver
    lags = [params.tau_1, params.tau_2, params.tau_35, params.tau_4];

    % using stiff DDE solver
    if settings.stiff
        stiff_start = tic;
        IIS_sol = dde15s_updated(@ddeIIS, lags, hist, params.tspan, options);  
        stiff_end = toc(stiff_start);
        disp(['Stiff Time: ', num2str(stiff_end)])
    % using regular DDE solver
    else
        reg_start = tic;
        IIS_sol = dde23(@ddeIIS, lags, hist, params.tspan, options);
        reg_end = toc(reg_start);
        disp(['Regular Time: ', num2str(reg_end)])
    end
   
    % evaluating output for 1000 evenly-spaced points
    xvals = linspace(IIS_sol.x(1), IIS_sol.x(end), 1000);
    yvals = deval(IIS_sol, xvals);

    % Peak Values
    [M_full, I_full] = max(yvals.');
  
    % Peak Values Display
    variables = {'V', 'X', 'Y', 'R', 'I'}';
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
        colororder(["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30"])
        xlabel('Time (h)', 'FontSize', 18);
        ylabel('Number of Cells', 'FontSize', 18);
        legend('$V$', '$X$', '$Y$', '$R$', '$I$', 'Interpreter', 'latex')
        saveas(h, fullfile('./simulations', title), 'png')
    end

% -------------------------------------------------------------------------

    function dydt = ddeIIS(t, y, Z)
    % Differential equations function for VXYRI.
        ylag1 = Z(:,1); %tau_1
        ylag2 = Z(:,2); %tau_2
        ylag3 = Z(:,3); %tau_35
        ylag4 = Z(:,4); %tau_4

        dydt = [params.k*ylag1(3)*params.k1_tilde^params.n_1/(params.k1_tilde^params.n_1 + ylag3(5)^params.n_1) - params.d_v*y(1)
                params.mu - params.d_x*y(2) - params.beta*y(2)*y(1) - params.k_ix*y(2)*y(5)
                params.beta*y(2)*y(1) - params.d_y*y(3) - params.k_iy*y(3)*y(5)
                params.k_ix*y(2)*y(5) + params.k_iy*y(3)*y(5) - params.d_r*y(4)
                params.k_i * ylag2(3) + (params.b_2*ylag4(5)^params.n_2)/(params.k_2^params.n_2 + ylag4(5)^params.n_2) - params.d_i*y(5)];
    end
end