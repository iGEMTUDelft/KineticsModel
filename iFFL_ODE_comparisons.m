clear all, clc
%%

tspan = [0 60];
x0 = ones(6,1).*0.1;
steady_states = zeros(100,1);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model(x, 0), tspan, x0);
    steady_states(i,1) = x(end,6);    
end

tspan = [0 60*250];
x0 = ones(6,1).*0.1;
steady_states_leaky = zeros(100,1);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model(x, 1), tspan, x0);
    steady_states_leaky(i,1) = x(end,6);    
end

x02 = zeros(4,1);
steady_states_reduced = zeros(100,1);
for i = 1:100
    c = i;
    [t, x] = ode15s(@(t,x) TALE_model_reduced(x, c), tspan, x02);
    steady_states_reduced(i,1) = x(end,4);    
end

analytical_sol = zeros(100,1);
for i = 1:100
    analytical_sol(i) = analytical_solution(i,1);
end

figure(1),
subplot(2,2,1);
plot(1:100,steady_states)
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
ylim([0 1.5])
title('Full model')
subplot(2,2,2);
plot(1:100,steady_states_reduced)
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
ylim([0 1.5])
title('Reduced model')
subplot(2,2,4)
plot(1:100,analytical_sol)
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
ylim([0 1.5])
title('Steady state solution')
subplot(2,2,3)
plot(1:100,steady_states_leaky)
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
title('Full Model with leaky terminator')

figure(3),
plot(1:100,steady_states_leaky)
hold on
plot(1:100,steady_states)
hold off
legend({"With leaky terminator", "Without leaky terminator"}, 'Location','northwest')
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
grid on
%%
tspan = [0 60*250];
x0 = ones(6,1).*0.1;
steady_states_leaky = zeros(100,1);
for i = 1:100
    x0(3) = i;
    [t, x] = ode15s(@(t,x) TALE_model(x, 0), tspan, x0);
    steady_states_leaky(i,1) = x(end,6);    
end
plot(1:100,steady_states_leaky)
xlabel('Copynumber of plasmid')
ylabel('GFP (nM)')
ylim([0 1.5])

%%
clc, clear all
n = 0.5:0.5:2;
c = 1:100;
G_values = zeros(length(c),length(n));
for i = 1:length(n)
    for ii = 1:length(c)
        G_values(ii,i) = analytical_solution(c(ii), n(i));
    end
end

figure(2),
plot(c, G_values)
title('Variation of n')
xlabel('copynumber')
ylabel('GFP (nM)')
ylim([0 5])
legend('n = 0.5','n = 1', 'n = 1.5', 'n = 2')
function dxdt = TALE_model(x, leaky)
    
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
    if leaky == 0
        dmG  = aGmax*P_G+aGmin*P_G_T-y_m*mG; %transcription GFP
    else
        dmG  = aGmax*P_G+aGmin*P_G_T-y_m*mG+c*a_T*1/60012;
    end
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

function G = analytical_solution(c,n)
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
    
    
    K_d = k_off/k_on;
     
    aG = aGmax*K_d;
    
    G = (c/c.^n).*((aG.*b_G.*y_T.^n*y_m.^n)/(a_T.^n.*b_T.^n.*y_G.*y_m));
end