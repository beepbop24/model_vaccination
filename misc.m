% endemic disease steady state for IIS expressed as a function of I^*
i_star = linspace(0, 5, 1000);
x_star = defaults.d_v*(defaults.d_y+defaults.k_iy*i_star).*(defaults.k1_tilde^defaults.n_1+i_star.^defaults.n_1)/(defaults.k1_tilde^defaults.n_1*defaults.beta*defaults.k);
full_expression = defaults.mu-defaults.d_x*x_star-1/defaults.k_i*(defaults.d_y+defaults.k_iy*i_star).*(defaults.d_i*i_star-defaults.b_2*i_star./(defaults.k_2+i_star))-defaults.k_ix*x_star.*i_star;

plot(i_star, full_expression)
hold on
yline(0, 'LineWidth', 1.5, 'Color', 'red')
hold off
