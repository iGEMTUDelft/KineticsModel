function dxdt = full_solution(x, aGmax, a_T, b_G, b_T)
    y_m = 0.2; %degradation mRNA (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min)
    aGmin = 0; %transcription when TALE bound
    y_G = log(2)/30; %degradation GFP
    
    mT = x(1);
    T = x(2);
    P_G = x(3);
    P_G_T = x(4);
    mG = x(5);
    G = x(6);
    
    n = 1; %cooperativity of binding
    c = P_G + P_G_T;
    
    dmT = c*a_T - y_m*mT; %change in mRNA TALE
    dT = b_T*mT - y_T*T - n*k_on*T^n*P_G + n*k_off*P_G_T; %translation TALE
    dP_G = n*k_off*P_G_T - n*k_on*T^n*P_G + n*y_T*P_G_T; %change in free promoter GFP
    dP_G_T = k_on*T^n*P_G - k_off*P_G_T - n*y_T*P_G_T; %change in bound promoter GFP
    dmG  = aGmax*P_G+aGmin*P_G_T-y_m*mG; %transcription GFP
    dG = b_G*mG-y_G*G;  
        
    dxdt = [dmT; dT; dP_G; dP_G_T; dmG; dG];
end
