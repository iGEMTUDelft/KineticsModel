clear all, clc
tspan = [0 60*250];
x0 = ones(6,1).*0.1;

%normal system
steady_states_normal = zeros(1,100);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model_normal(x), tspan, x0);
    steady_states_normal(1,i) = x(end,6);    
end

%leaky system
steady_states_leaky = zeros(1,100);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model_leaky(x), tspan, x0);
    steady_states_leaky(1,i) = x(end,6);    
end
    

N = 1;
M = 100000;
steady_states = zeros(M,100);
steady_states_avg = zeros(N,100);
for j = 1:N
    for i = 1:M
        r = 612.*rand(1,1);
        if r<=1
            steady_states(i,:) = steady_states_leaky;
        else
            steady_states(i,:) = steady_states_normal;
        end
    end
    steady_states_avg(j,:) = sum(steady_states,1)./M;
end

plot(1:100,steady_states_avg)
hold on
plot(1:100,steady_states_normal)
hold off
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
legend("leaky terminator", 'perfect terminator')
title('Insulator efficiency')


function dxdt = TALE_model_normal(x)
    
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

function dxdt = TALE_model_leaky(x)
    
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
    
    a = 0;
    b = 612;
    r = (b-a).*rand(1,1) + a;
    
    if r<=1
        leaky = 1;
    else
        leaky = 0;
    end
    
    dmT = c*a_T - y_m*mT; %change in mRNA TALE
    dT = b_T*mT - y_T*T - n*k_on*T^n*P_G + n*k_off*P_G_T + (n-1)*n*y_T*P_G_T; %translation TALE
    dP_G = k_off*P_G_T - k_on*T^n*P_G + n*y_T*P_G_T; %change in free promoter GFP
    dP_G_T = k_on*T^n*P_G - k_off*P_G_T - n*y_T*P_G_T; %change in bound promoter GFP
    dmG  = aGmax*P_G+aGmin*P_G_T-y_m*mG+c*a_T*1/612;
    dG = b_G*mG-y_G*G;  
        
    dxdt = [dmT; dT; dP_G; dP_G_T; dmG; dG];
end