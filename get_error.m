function error = get_error(a_T, b_T, c_min)
    y_m = 0.2; %degradation mRNA (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min) 
    
    beta = b_T*a_T/(y_T*y_m);
    
    K_d = k_off./k_on;
    error = K_d/(beta*c_min);
end