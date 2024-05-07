
function model_test
tau_1 = 12;
tau_2 = 8;
tau_35 = 12;
tau_4 = 9;

sol = dde23(@ddex1de, [tau_1, tau_2, tau_35, tau_4], @ddex1hist, [0, 1000]);
disp(sol.y)
figure;
plot(sol.x, sol.y)
ylim([0 1])
title('Innate Immune System Dynamics');
xlabel('time (h)');
ylabel('y (# of cells)');
legend('V', 'X', 'Y', 'R', 'I')

% --------------------------------------------------------------------------

function s = ddex1hist(t)
% Constant history function for VXYRI.
s = [0; 100; 1; 0; 0];

% --------------------------------------------------------------------------

function dydt = ddex1de(t, y, Z)
% Differential equations function for VXYRI.
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
params.b_2 = 80; %0.07;
params.k_2 = 0.1;
params.n_2 = 1;
params.d_i = 0.7;

ylag1 = Z(:,1); %tau_1
ylag2 = Z(:,2); %tau_2
ylag3 = Z(:,3); %tau_35
ylag4 = Z(:,4); %tau_4

dydt = [params.k*ylag1(3)/(params.k1_tilde + ylag3(5)^params.n_1) - params.d_v*y(1)
        params.mu - params.d_x*y(2) - params.beta*y(2)*y(1) - params.k_ix*y(2)*y(5)
        params.beta*y(2)*y(1) - params.d_y*y(3) - params.k_iy*y(3)*y(5) 
        params.k_ix*y(2)*y(5) + params.k_iy*y(3)*y(5) - params.d_r*y(4)
        params.k_i * ylag2(3) + (params.b_2*ylag4(5)^params.n_2)/(params.k_2 + ylag4(5)^params.n_2) - params.d_i*y(5)];
