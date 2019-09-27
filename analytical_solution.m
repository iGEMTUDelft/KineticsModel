function G = analytical_solution(aGmax, a_T, b_G, b_T)
    y_m = 0.2; %degradation mRNA (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min) 
    y_G = log(2)/30; %degradation GOI
    
    K_d = k_off./k_on;
     
    aG = aGmax.*K_d;
    
    G = (aG.*b_G.*y_T*y_m)/(a_T.*b_T.*y_G.*y_m);
    
end
