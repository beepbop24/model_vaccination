% parameter values for V
params.k = 40;
params.k1_tilde = 39.6;
params.n_1 = 2;
params.d_v = 1/2;

% parameter values for X
params.mu = 1/312;
params.d_x = 1/312;
params.beta = 0.5;
params.k_ix = 0.5;

% parameter values for Y
params.d_y = 1/12;
params.k_iy = 0.5;

% parameter values for R
params.d_r = 1/50;

% parameter values for I
params.k_i = 0.3;
params.b_2 = 80;
params.k_2 = 0.1;
params.n_2 = 2;
params.d_i = 0.7;

% delays
params.tau_1 = 10;
params.tau_2 = 8;
params.tau_35 = 12;
params.tau_4 = 9;

%% DFE with I^* = 0
params.n1 = 1;
params.n2 = 1;

% contour plot -- zoomed out
h1 = figure(1);
[X, Y] = meshgrid(linspace(-1, 3, 300), linspace(-3, 3, 300));
I_star = 0;
charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, params, I_star), X+1i*Y);
contour(X, Y, charpoly_Df2, [0.0001 0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10 25 50 100])
hold on
plot(-params.d_r, 0, '.')
hold on
plot(-params.d_x, 0, '.')
%hold on
%plot(-params.d_i, 0, '.')
hold off
title('Value of the Characteristic Polynomial of the Linearized System Around DFE for $n_2=1$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$Re(\lambda)$','Interpreter','latex')
ylabel('$Im(\lambda)$','Interpreter','latex')
saveas(h1, 'contour_DFE', 'jpeg')

% contour plot -- zoomed in
h2 = figure(2);
[X, Y] = meshgrid(linspace(-0.5, 0.5, 300), linspace(-0.5, 0.5, 300));
charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, params, I_star), X+1i*Y);
contour(X, Y, charpoly_Df2, [0.00001 0.0005 0.0001 0.0005 0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10])
hold on
plot(-params.d_r, 0, '.')
hold on
plot(-params.d_x, 0, '.')
%hold on
%plot(-params.d_i, 0, '.')
hold off
legend('', '$\lambda_2 = - d_R$', '$\lambda_1 = -d_X$', 'Interpreter', 'latex', 'Fontsize', 16)
title('Value of the Characteristic Polynomial of the Linearized System Around DFE for $n_2=1$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$Re(\lambda)$','Interpreter','latex')
ylabel('$Im(\lambda)$','Interpreter','latex')
saveas(h2, 'contour_DFE_zoomedin', 'jpeg')

% plot along real axis
h3 = figure(3);
charpoly_real = arrayfun(@(lambda) characEq(lambda, params, I_star), linspace(-1, 3, 50));
plot(linspace(-1, 3, 50), charpoly_real)
axis([-1 3 0 500])


%% DFE WITH I^* = b2/dI - K2
params.n_1 = 1;
params.n_2 = 1;
% contour plot
h4 = figure(4);
[X, Y] = meshgrid(linspace(-0.8, 0.2, 300), linspace(-0.8, 0.8, 300));
I_star = params.b_2/params.d_i - params.k_2;
charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, params, I_star), X+1i*Y);
contour(X, Y, charpoly_Df2, [0.0001 0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10 25 50 100])
title('Value of the Characteristic Polynomial of the Linearized System Around DFE$^*$ for $n_2=1$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$Re(\lambda)$','Interpreter','latex')
ylabel('$Im(\lambda)$','Interpreter','latex')
saveas(h4, 'contour_DFE_star', 'jpeg')

% plot along real axis
h5 = figure(5);
charpoly_real = arrayfun(@(lambda) characEq(lambda, params, I_star), linspace(-0.8, 0.2, 50));
plot(linspace(-0.8, 0.2, 50), charpoly_real)
%axis([-1 3 0 500])


%% FUNCTION THAT COMPUTES ABSOLUTE VALUE OF CHARACTERISTIC EQUATION
function charpoly_Df = characEq(lambda, params, I_star)
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


%char_eq = @(lambda) abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
        %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
        %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));

%test_value = fsolve(char_eq, 0.01+0.01*1i);
