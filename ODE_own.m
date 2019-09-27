
tspan = [0 60*16];
x0 = zeros(9,1);
copynumbers = 1:200;

steady_states = zeros(length(copynumbers),1);
for i = 1:length(copynumbers)
        x0(1) = copynumbers(i);
        x0(4) = copynumbers(i);
        [t, x] = ode15s(@(t,x) ODE(x), tspan, x0);
        steady_states(i) = x(end,8);
end
plot(copynumbers, steady_states)
xlabel('Copy number');
ylabel('GFP (nM)');
% plot(t, x)
% legend("T7p","mT", "T", "P_G", "P_G_T7", "P_G_T", "mG", "G", "T7p_T")
function dxdt = ODE(x)
t7_pool = 10^4;

y_m = 0.2; %degradation mRNA (1/min)
b_T = 0.019; %translation TALE (1/min)
y_T = log(2)/30; %degradation TALE (1/min)
y_G = log(2)/30; %degradation TALE (1/min)
b_G = 3.65; %translation GFP

n = 1;

k_on = 100; %1/(nM*min)
k_off = 0.1; %1/min

k_on_T7 = 0.000131; %T7RNAP binding rate to T7 promoter (nM^-1 min^-1)
k_off_T7 = 0.0003311; %T7RNAP dissociation rate from T7 promoter (min^-1)

TALE_len = 2640;
GFP_len = 705; 

k_elong = 16293; %nt/min

T7p = x(1);
mT = x(2);
T = x(3);
P_G = x(4);
P_G_T7 = x(5);
P_G_T = x(6);
mG = x(7);
G = x(8);
T7p_T7 = x(9);

T7RNAP = t7_pool - P_G_T7 - T7p_T7;

dT7p_T7 = k_on_T7*T7p*T7RNAP - k_off_T7*T7p_T7 - k_elong/TALE_len * T7p_T7; %binding of T7RNAP to TALE promoter
dT7p = -dT7p_T7;
dmT = k_elong/TALE_len * T7p_T7 - mT*y_m; %transcription of TALE
dT = b_T*mT - y_T*T - n*k_on*T^n*P_G + n*k_off*P_G_T; %translation of TALE
dP_G = n*k_off*P_G_T - n*k_on*T^n*P_G + n*y_T*P_G_T - k_on_T7*P_G*T7RNAP + k_off_T7*P_G_T7; %change in free promoter GFP
dP_G_T7 = -k_off_T7*P_G_T7 + k_on_T7*P_G*T7RNAP - k_elong/GFP_len * P_G_T7; 
dP_G_T = k_on*T^n*P_G - k_off*P_G_T - n*y_T*P_G_T; %change in promoter GFP bound by TALE
dmG  = k_elong/GFP_len * P_G_T7 - y_m*mG; %transcription GFP
dG = b_G*mG-y_G*G; %translation of GFP

dxdt = [dT7p; dmT; dT; dP_G; dP_G_T7; dP_G_T; dmG; dG; dT7p_T7];

end

