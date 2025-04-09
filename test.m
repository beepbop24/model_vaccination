beta_k_test

%%
load('dataM.mat')

cMap=jet(256); %set the colomap using the "jet" scale
F2=figure(1);
[c,h]=contourf(xrow,ycol,BDmatrix,50);
set(h, 'edgecolor','none');

xlim([0.0352 0.3872]);
ylim([0.0352 0.3872]);

colormap(cMap);
cb=colorbar;
caxis([0.7 0.96]);
% box on;
hold on;

%%
beta_vals = linspace(0.015, 0.025, 10);

a0_vals = linspace(0, log10(300), 10);
bm0_vals = linspace(0, log(905), 10);
t0_vals = linspace(0, 7, 10);

h = figure();
[c, k] = contourf(a0_vals, beta_vals, heat_matrix, 15, 'LineWidth', 1.5);
set(k, 'edgecolor','none');
colormap(turbo);
colorbar;
ax = gca;
ax.FontSize = 16;
xlabel('$A^0$','Interpreter','latex', 'FontSize', 24)
ylabel('$\beta$','Interpreter','latex', 'FontSize', 24)
%%
function heat_matrix = beta_k_test
    params = struct();
    params.tspan = [0, 800];
    options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-9);

    beta_vals = linspace(0.015, 0.025, 100);

    a0_vals = linspace(0, log10(300), 100);
    bm0_vals = linspace(0, log(905), 100);
    t0_vals = linspace(0, 7, 100);

    settings = struct();
    settings.plot = false;

    m = length(beta_vals);
    n = length(a0_vals);

    heat_matrixa0 = zeros(m,n);
    for i=1:m
        params.beta = beta_vals(i);
        for j=1:n
            init_hist = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; 4; bm0_vals(j); 0; a0_vals(j);];
            [sol, I_full, M_full] = full_model("default_params_full.mat", params, @ddefullhist, options, settings);
            
            heat_matrixa0(i, j) = min((10^M_full(3))/(5.25*10^9), 1);
        end
    end
    writematrix(heat_matrixa0, 'a0beta.txt')

    heat_matrixtm0 = zeros(m,n);
    for i=1:m
        params.beta = beta_vals(i);
        for j=1:n
            init_hist = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; t0_vals(j); log(90.5); 0; log(30);];
            [sol, I_full, M_full] = full_model("default_params_full.mat", params, @ddefullhist, options, settings);
            
            heat_matrixtm0(i, j) = min((10^M_full(3))/(5.25*10^9), 1);
        end
    end
    writematrix(heat_matrixtm0, 'tm0beta.txt')

    heat_matrixatm = zeros(m,n);
    params.beta = 2.175*10^(-2);
    for i=1:m
        for j=1:n
            init_hist = [0; log10(5.25*10^9); log10(5250); 0; 0; 0; 0; t0_vals(i); bm0_vals(j); 0; a0_vals(j);];
            [sol, I_full, M_full] = full_model("default_params_full.mat", params, @ddefullhist, options, settings);
            
            heat_matrixatm(i, j) = min((10^M_full(3))/(5.25*10^9), 1);
        end
    end
    writematrix(heat_matrixatm, 'a0tm0.txt')

    % contour plot
    %h = figure();
    %contourf(a0_vals, beta_vals, heat_matrix, [0 0.01 0.05 0.1 0.25 0.4 0.5 0.75 0.9 1], 'LineWidth', 1.5);
    %set(h, 'edgecolor','none');
    %colormap(turbo);
    %ax = gca;
    %ax.FontSize = 16;
    %xlabel('$T_M^0$','Interpreter','latex', 'FontSize', 24)
    %ylabel('$\beta$','Interpreter','latex', 'FontSize', 24)
    %xticks([0 0.01 0.02 0.03 0.04 0.05])
    %yticks([0 0.2 0.4 0.6 0.8 1])
    %contourcmap("turbo",[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1], "Colorbar","on","Location", "horizontal", "TitleString","% of Lungs Infected")
    %saveas(h, fullfile('./betak', title), 'png')

    %%
    function s = ddefullhist(t)
        s = init_hist;
    end
end