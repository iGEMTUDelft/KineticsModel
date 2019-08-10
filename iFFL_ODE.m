clear all, clc

tspan = [0 60*210];
x0 = zeros(6,1);
steady_states = zeros(100,1);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model(x), tspan, x0);
    steady_states(i,1) = x(end,6);    
end

plot(1:100,steady_states)
hold on

x02 = zeros(4,1);
steady_states = zeros(100,1);
for i = 1:100
    c = i;
    [t, x] = ode15s(@(t,x) TALE_model_reduced(x, c), tspan, x02);
    steady_states(i,1) = x(end,4);    
end

plot(1:100,steady_states)
%ylim([0 2])
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
title('GFP expression vs copy number')

function dxdt = TALE_model(x)
    
    a_T = 1.03; %transcription TALE (nM/min)
    y_m = 0.2; %degradation mRNA (1/min)
    b_T = 0.019; %translation TALE (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min)
    aGmax = 3.78; %transcription when no TALE bound (nM/min)
    aGmin = 0; %transcription when TALE bound
    b_G = 3.65; %translation GFP
    y_G = log(2)/30; %degradation GFP
    
    mT = x(1);
    T = x(2);
    P_G = x(3);
    P_G_T = x(4);
    mG = x(5);
    G = x(6);
    
    n = 1;
    c = P_G + P_G_T;
    
    dmT = c*a_T - y_m*mT; %change in mRNA TALE
    dT = b_T*mT - y_T*T - n*k_on*T^n*P_G + n*k_off*P_G_T + (n-1)*n*y_T*P_G_T; %translation TALE
    dP_G = k_off*P_G_T - k_on*T^n*P_G + n*y_T*P_G_T; %change in free promoter GFP
    dP_G_T = k_on*T^n*P_G - k_off*P_G_T - n*y_T*P_G_T; %change in bound promoter GFP
    dmG  = aGmax*P_G+aGmin*P_G_T-y_m*mG; %transcription GFP
    dG = b_G*mG-y_G*G;  
        
    dxdt = [dmT; dT; dP_G; dP_G_T; dmG; dG];
end

function dxdt = TALE_model_reduced(x, c)
    
    a_T = 1.03; %transcription TALE (nM/min)
    y_m = 0.2; %degradation mRNA (1/min)
    b_T = 0.019; %translation TALE (1/min)
    k_on = 100; %1/(nM*min)
    k_off = 0.1; %1/min
    y_T = log(2)/30; %degradation TALE (1/min) is only dependent on dilution
    aGmax = 3.78; %transcription when no TALE bound (nM/min)
    aGmin = 0; %transcription when TALE bound
    b_G = 3.65; %translation GFP
    y_G = log(2)/30; %degradation GFP is only dependent on dilution
    
    mT = x(1);
    T = x(2);
    mG = x(3);
    G = x(4);
    
    K_d = k_off/k_on;
    
    n = 1;
    
    dmT = c*a_T - y_m*mT; %change in mRNA TALE
    dT = b_T*mT - y_T*T;
    dmG = c*(aGmin + (aGmax - aGmin)*(K_d^n/(K_d^n + T^n))) - y_m*mG;
    dG = b_G*mG-y_G*G;  
        
    dxdt = [dmT; dT; dmG; dG];
end
