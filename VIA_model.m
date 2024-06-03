function VIA_model
    tau_1 = 12;
    tau_2 = 8;
    tau_3 = 7;
    tau_4 = 9;
    tau_5 = 5;

    sol = dde23(@ddex1de, [tau_1, tau_2, tau_3, tau_4, tau_5], @ddex1hist, [0, 100]);
    figure;
    plot(sol.x, sol.y)
    title('VIA Dynamics');
    xlabel('time (h)');
    ylabel('y (# of cells)');
    legend('V', 'I', 'A')

% --------------------------------------------------------------------------

    function s = ddex1hist(t)
    % Constant history function for VIA.
        s = [1000; 10; 2];

% --------------------------------------------------------------------------

    function dydt = ddex1de(t, y, Z)
    params.k_1 = 4;
    params.k_2 = 0.3;
    params.k_3 = 0.1;
    params.d_1 = 0.1;
    params.d_2 = 0.7;
    params.d_3 = 0.12;
    params.b_1 = 10;
    params.b_2 = 80;
    params.K_1 = 33;
    params.K_2 = 0.1;
    params.n_1 = 1;
    params.n_2 = 1;

    ylag1 = Z(:,1); %tau_1
    ylag2 = Z(:,2); %tau_2
    ylag3 = Z(:,3); %tau_3
    ylag4 = Z(:,4); %tau_4
    ylag5 = Z(:,5); %tau_5

    dydt = [params.k_1*ylag1(1)*params.b_1*params.K_1^params.n_1/...
            (params.K_1^params.n_1 + ylag5(3)^params.n_1) - params.d_1*y(1)
            params.k_2*ylag2(1)+ params.b_2*ylag4(2)^params.n_2/...
            (params.K_2^params.n_2 + ylag4(2)^params.n_2) - params.d_2*y(2)
            params.k_3*ylag3(2)-params.d_3*y(3)];
