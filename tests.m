
IIS_sol = IIS_model("param_values.mat", @ddeIIShist);

function s = ddeIIShist(t)
    % Constant history function for VXYRI.
    s = [0; 7*10^6; 1; 0; 0];
end