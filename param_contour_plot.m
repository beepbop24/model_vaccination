function heat_matrix = param_contour_plot(model, params_file, params, hist, options, beta_vals, k_vals, option, title)
    % function that returns numerical simulations of the full immune system
    % for reinfections with INPUTS -- 
    % model: 'full' or 'iis'
    % params_file: default parameters (.mat file)
    % params = struct() of specified parameters (can be empty)
    % options = ddeset options
    % hist:  history function of the DDE (function of t)
    % beta_vals, k_vals: linspace of parameter values for beta, k_v
    % option: 1 (peak number of infected cells) or 2 (total number of infected cells)
    % title: name of png file generated

    settings = struct();
    settings.plot = false;

    m = length(beta_vals);
    n = length(k_vals);

    heat_matrix = zeros(m,n);

    for i=1:m
        params.beta = beta_vals(i);
        for j=1:n
            params.k = k_vals(j);
            if strcmp(model, "full")
                [sol, I_full, M_full] = full_model(params_file, params, hist, options, settings);
            elseif strcmp(model, "iis")
                [sol, M_full] = IIS_model(params_file, params, hist, options, settings);
            else
                error('Not a valid model');
            end

            if option==1
                heat_matrix(i, j) = (10^M_full(3))/(5.25*10^9);
                
            elseif option==2
                new_cells = params.beta.*sol.y(1,:).*sol.y(2,:);
                heat_matrix(i, j) = trapz(new_cells);
            else
                error("Not a valid option")
            end
        end
    end

    % contour plot
    h = figure();
    if option==1
        contourf(beta_vals, k_vals, heat_matrix, [0 0.01 0.05 0.1 0.25 0.4 0.5 0.75 0.9 1], 'LineWidth', 1.5);
    elseif option==2
        % [0 3 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]
        contourf(beta_vals, k_vals, heat_matrix,[0 500 1000 2000 3000 4000 5000 6000 7000 8000 9000], 'LineWidth', 1.5);
    end
    colormap(turbo);
    ax = gca;
    ax.FontSize = 16;
    xlabel('$\beta$','Interpreter','latex', 'FontSize', 24)
    ylabel('$k_V$','Interpreter','latex', 'FontSize', 24)
    xticks([0 0.01 0.02 0.03 0.04 0.05])
    yticks([0 0.2 0.4 0.6 0.8 1])
    %contourcmap("turbo",[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], "Colorbar","on","Location", "horizontal", "TitleString","% of Lungs Infected")
    saveas(h, fullfile('./betak', title), 'png')
end