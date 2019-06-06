
tspan = [0 1500];
x0 = zeros(8,1);
x0(3,1) = 490;
x0(4,1) = 7.5556e-08;
x0(5,1) = 0.5; % definitely wrong
x0(6,1) = 0.26;
x0(7,1) = 210;

[t, x] = ode15s(@(t,x) ODE(x), tspan, x0);

figure(1);
plot(t, x(:,1))
title("T7RNAP")
figure(2);
plot(t, x(:,2))
title("mRNA T7RNAP")

function dxdt = ODE(x)
    mRNA_T7RNAP = x(1);
    T7RNAP = x(2);
    eRNAP = x(3); 
    Ribosome = x(4);
    LacRepressor = x(5); 
    Pr_T7RNAP_free = x(6);
    C1 = x(7);
    IPTG = x(8);
  
    %parameters
    k0 = 2.76e12; %1/(uM min)
    k1 = 12; %1/min
    k2 = 1042; %initiations per min
    k4 = 3120; %nt/min
    k5 = 1200; %aa/min
    k6 = 0.000131; %1/(uM min)
    k7 = 0.000331; %1/min
    k8 = 2580; %nt/min
    k12 = 0.0059; %uM/min (ln(2)/30)
    k13 = 0.036; %uM/min (ln(2)/5)
    k14 = 1200; %aa/min
    L_T7RNAP = 2697; %nts
    
    r0 = k0*LacRepressor*IPTG;
    r1 = k1*Pr_T7RNAP_free;
    r2 = k2*Pr_T7RNAP_free*eRNAP;
    r4 = k4*C1/L_T7RNAP;
    r5 = k5*mRNA_T7RNAP*Ribosome;
    r13 = k12*T7RNAP;
    r14 = k13*mRNA_T7RNAP;
    r18 = k12*IPTG;
    
    %ODEs
    d_mRNA_T7RNAP = r4-r14;
    d_T7RNAP = r5-r13;
    d_eRNAP = -r2 + r4; 
    d_Ribosome = 0;
    d_LacRepressor = -r0+r1;
    d_Pr_T7RNAP_free = r0 - r1 + r4;
    d_C1 = r2 - r4;
    d_IPTG = -r0 + r1 - r18;
    
    dxdt = [d_mRNA_T7RNAP; d_T7RNAP; d_eRNAP; d_Ribosome; d_LacRepressor; d_Pr_T7RNAP_free; d_C1; d_IPTG];
end