% parameter values for V
params.k = 40;
params.k1_tilde = 39.6;
params.n_1 = 1;
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
params.n_2 = 1;
params.d_i = 0.7;

% delays
params.tau_1 = 10;
params.tau_2 = 8;
params.tau_35 = 12;
params.tau_4 = 9;

% contour plot
h1 = figure(1);
[X, Y] = meshgrid(linspace(-1, 3, 300), linspace(-3, 3, 300));
charpoly_Df2 = arrayfun(@(lambda) characEq(lambda, params), X+1i*Y);
contour(X, Y, charpoly_Df2, [0.001 0.01 0.025 0.05 0.1 0.25 0.5 1 1.25 1.5 2 2.5 5 10 25 50 100])

% plot along real axis
h2 = figure(2);
[X, Y] = meshgrid(linspace(-1, 3, 300), linspace(-3, 3, 300));
charpoly_real = arrayfun(@(lambda) characEq(lambda, params), X+1i*Y);
plot(linspace(-1, 3, 300), charpoly_real)
axis([-1 3 0 500])

function charpoly_Df = characEq(lambda, params)
    Df_x0 = [-params.d_v, 0, params.k/params.k1_tilde*exp(-lambda*params.tau_1), 0, 0;
            -params.beta*params.mu/params.d_x, - params.d_x, 0, 0, -params.k_ix*params.mu/params.d_x;
            params.beta*params.mu/params.d_x, 0, -params.d_y, 0, 0;
            0, 0, 0, -params.d_r, params.k_ix*params.mu/params.d_x;
            0, 0, params.k_i*exp(-lambda*params.tau_2), 0, -params.d_i];

    charpoly_Df = abs(det(Df_x0-lambda*eye(5)));
    %charpoly_Df = abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
        %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
        %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));
end


%char_eq = @(lambda) abs((params.d_x+lambda)*(params.d_r+lambda)*(params.d_i+lambda)*...
        %(-(params.d_v+lambda)*(params.d_y+lambda)+params.k/params.k1_tilde*...
        %exp(-lambda*params.tau_1)*params.beta*params.mu/params.d_x));

%test_value = fsolve(char_eq, 0.01+0.01*1i);
