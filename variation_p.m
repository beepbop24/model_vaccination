function variation_p = sensitivity_analysis(x1, x2, delta_p)
    delta_x = abs(x1-x2)/(x1);
    variation_p = delta_x/delta_p;
end