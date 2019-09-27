function G = analytical_solution2(aGmax, a_T, b_G, b_T, c)
    y_m = 0.2; %degradation mRNA (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min) 
    y_G = log(2)/30; %degradation GOI
    
    K_d = k_off./k_on;
    
    alpha = b_G*aGmax/(y_G*y_m); 
    beta = b_T*a_T/(y_T*y_m);
    G = alpha*(c*K_d/(K_d+beta*c));
end
