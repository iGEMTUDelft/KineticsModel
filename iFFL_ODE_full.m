clear all, clc

tspan = [0 16*60];
x0 = zeros(6,1);
copynumbers = 1:100;

steady_states = zeros(length(copynumbers),1);
for i = 1:length(copynumbers)
    x0(3) = copynumbers(i);
    [t, x] = ode15s(@(t,x) TALE_model(x), tspan, x0);
    steady_states(i) = x(end,6);
end

analytical = analytical_solution().*ones(length(copynumbers),1);

figure(1),
hold on

model = line(copynumbers, steady_states); 
anal = line(copynumbers, analytical);

set(model, 'Color', [0 0 0])

set(model, 'LineWidth', 2)
set(anal, 'Color', [.75 .75 1])
set(anal, 'LineWidth', 2)

hLegend = legend('Full model solution', 'Analytical Solution with simplifications');

hTitle = title('Steady state solution GOI vs copy number');
hXLabel = xlabel('Copy number');
hYLabel = ylabel('GFP steady-state (nM)');

set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 10)
set([hXLabel, hYLabel], 'FontSize', 16)
set(hTitle, 'FontSize', 18, 'FontWeight' , 'bold')
set(gca, 'YGrid', 'on', 'XGrid', 'on')
ylim([10 18])
hold off

function dxdt = TALE_model(x)
    
    a_T = 1.03; %transcription TALE (nM/min)
    y_m = 0.2; %degradation mRNA (1/min)
    b_T = 0.44; %translation TALE (1/min)
    k_on = 9.85; %1/(nM*min)
    k_off = 2.19; %1/min
    y_T = 0.047; %degradation TALE (1/min)
    aGmax = 3.78; %transcription when no TALE bound (nM/min)
    aGmin = 0; %transcription when TALE bound
    b_G = 3.65; %translation GFP
    y_G = 0.019; %degradation GFP
    
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
function G = analytical_solution()
   a_T = 1.03; %transcription TALE (nM/min)
    y_m = 0.2; %degradation mRNA (1/min)
    b_T = 0.44; %translation TALE (1/min)
    k_on = 9.85; %1/(nM*min)
    k_off = 2.19; %1/min
    y_T = 0.047; %degradation TALE (1/min)
    aGmax = 3.78; %transcription when no TALE bound (nM/min)
    aGmin = 0; %transcription when TALE bound
    b_G = 3.65; %translation GFP
    y_G = 0.019; %degradation GFP
    K_d = k_off/k_on;
     
    aG = aGmax*K_d;
    
    G = (aG.*b_G.*y_T*y_m)/(a_T.*b_T.*y_G.*y_m);
end

