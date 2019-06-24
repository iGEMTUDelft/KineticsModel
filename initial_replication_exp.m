clear all
close all

%%Parameters
%Basics
tspan = [0 200]; %Time span in minutes
numElements = 14; %Number of elements in the model
names = {'T7RNAP initiation complex', 'mRNA T7RNAP', 'Ribosomes T7RNAP', ...
    'T7RNAP', 'DNAP initiation complex', 'mRNA DNAP', 'Ribosome DNAP', ...
    'DNAP', 'GFP initiation complex', 'mRNA GFP', 'Ribosomes GFP', ...
    'GFP', 'Replication complex', 'Plasmid'}; %Names of elements in the model

%Initial concentrations
IPTG = true; %Boolean indicating whether or not IPTG has been added to the medium
x0 = zeros(1,numElements); %Initiating vector containing starting concentrations of all elements, most are zero
%x0(14) = 2.6; %Starting concentration of plasmid in nM
x0(14) = 1; %Copies of plasmid

%%Model Calculation
%ODEsolver that evaluates Initial_Replication_function over time specified by 'tspan'
%and with initial conditions listed in 'x0'. Additionally, the presence of
%IPTG is an input variable for the function. 
[t, y] = ode15s(@(t,x) ODE(x,IPTG), tspan, x0);

%%
%Plot concentrations over time
figure('Name','Orthogonal Replication');
for i = 1:size(y,2)
    subplot(4,4,i)
    plot(t, y(:,i), 'LineWidth', 2);
    title(names(i));
    ylabel('Number of copies')
    xlabel('Time (min)')
end

%%
function dxdt = ODE(x, IPTG)

%%Constants

%Constant concentrations
% eRNAP = 4930; %Endogenous RNAP concentration in nM
% Pr_DNAP = 2.6; %DNAP promoter concentration in nM
% Pr_GFP = 2.6; %GFP promoter concentration in nM
eRNAP = 1929; %Number of endogenous RNAP copies
Pr_DNAP = 1; %Number of DNAP promoter copies
Pr_GFP = 1; %Number of GFP promoter copies

%Presence of IPTG is evaluated and active T7RNAP promoter (Lac promoter)
%concentration is updated accordingly
if IPTG == false
    %If there is no IPTG added, there is no active T7RNAP promoter
    Pr_T7RNAP = 0; 
elseif IPTG == true
    %If IPTG is added, the T7RNAP promoter is activated, leading to a 
    %concentration of 2.6 nM
    %Pr_T7RNAP = 2.6; 
    Pr_T7RNAP = 1;
end


%%Variables

%Concentration variables of all elements with changing concentrations
C1 = x(1); %Transcription initiation complex of T7RNAP
mRNA_T7RNAP = x(2); %T7RNAP mRNA 
Ribosome_T7RNAP = x(3); %Ribosomes translating T7RNAP mRNA
T7RNAP = x(4); %T7 RNA polymerase
C2 = x(5); %Transcription initiation complex of DNAP
mRNA_DNAP = x(6); %Phi29 DNAP mRNA
Ribosome_DNAP = x(7); %Ribosomes translating Phi29 DNAP mRNA
DNAP = x(8); %Phi29 DNA polymerase
C4 = x(9); %Transcription initiation complex of GFP
mRNA_GFP = x(10); %GFP mRNA
Ribosome_GFP = x(11); %Ribosomes translating GFP mRNA
GFP = x(12); %GFP 
C3 = x(13); %Replication complex
Plasmid = x(14); %Linear plasmid holding GFP

%Mole balances to calculate concentration of unoccupied elements when
%applicable
eRNAP_free = eRNAP - C1; %Endogenous RNAP that is not transcribing
T7RNAP_free = T7RNAP - C2; %T7RNAP that is not transcribing
DNAP_free = DNAP - C3; %DNAP that is not replicating
Pr_T7RNAPfree = Pr_T7RNAP - C1; %Unoccupied T7RNAP promoter
Pr_DNAPfree = Pr_DNAP - C2; %Unoccupied DNAP promoter
Pr_GFPfree = Pr_GFP - C4; %Unoccupied GFP promoter
Plasmid_free = Plasmid - C3; %Plasmid that is not being replicated

%Rate constants
k2 = 1042; %Transcription initiation rate by endogenous RNAP (min^-1)
k4 = 1.16; %Transcription elongation rate of T7RNAP by endogenous RNAP (min^-1)
k5 = 1.33; %Translation rate of T7RNAP mRNA (min^-1)
k6 = 1.6909776 * 10^12; %T7RNAP binding rate to T7 promoter (nM^-1 min^-1)
%k6 = 4320*10^3; %nM^-1 min^-1 
k7 = 240; %T7RNAP dissociation rate from T7 promoter (min^-1)
k8 = 1.50; %Transcription elongation rate of DNAP by T7RNAP (min^-1)
k9 = 1.56572*10^12; %Phi29 DNAP binding rate to oriL or oriR (nM^-1 min^-1)
%k9 = 4000*10^3; %nM^-1 min^-1 
k10 = 200; %Phi29 DNAP dissociation rate from oriL or oriR (min^-1)
k11 = 0.46; %Plasmid replication rate by Phi29 DNAP (min^-1)
k12 = 0.023; %Protein degradation rate (min^-1)
k13 = 0.139; %mRNA degradation rate (min^-1)
k14 = 2.09; %Translation rate of Phi29 DNAP (min^-1)
k15 = 3.81; %Transcription elongation rate of GFP by T7RNAP (min^-1)
k16 = 5.31; %Translation rate of GFP (min^-1)


%%Calculation of concentration changes

%Rate laws

%Pr_T7RNAPfree + eRNAP_free -> C1
r2 = k2 * Pr_T7RNAPfree * eRNAP_free;
%C1 -> Pr_T7RNAPfree + eRNAP_free + mRNA_T7RNAP
r4 = k4 * C1;
%mRNA_T7RNAP + Ribosome_T7RNAP -> mRNA_T7RNAP + Ribosome_T7RNAP + T7RNAP
r5 = k5 * mRNA_T7RNAP * Ribosome_T7RNAP;
%Pr_DNAPfree + T7RNAP_free <--> C2
r6 = k6 * Pr_DNAPfree * T7RNAP_free;
r7 = k7 * C2;
%C2 -> Pr_DNAPfree + T7RNAP_free + mRNA_DNAP
r8 = k8 * C2;
%mRNA_DNAP + Ribosome_DNAP -> mRNA_DNAP + Ribosome_DNAP + DNAP
r9 = k14 * mRNA_DNAP * Ribosome_DNAP;
%Plasmid_free + 2 DNAP_free <--> C3
r10 = k9 * Plasmid_free * DNAP_free^2;
r11 = k10 * C3;
%C3 -> 2 Plasmid_free + 2 DNAP_free + Pr_GFP
r12 = k11 * C3;
%T7RNAP -> Ø
r13 = k12 * T7RNAP;
%mRNA_T7RNAP -> Ø
r14 = k13 * mRNA_T7RNAP;
%DNAP -> Ø
r15 = k12 * DNAP;
%mRNA_DNAP -> Ø
r16 = k13 * mRNA_DNAP;
%Plasmid -> Ø
r17 = k12 * Plasmid;
%GFP -> Ø
r18 = k12 * GFP;
%mRNA_GFP -> Ø
r19 = k13 * mRNA_GFP;
%Pr_GFPfree + T7RNAP_free <--> C4
r20 = k6 * Pr_GFPfree * T7RNAP_free;
r21 = k7 * C4;
%C4 -> Pr_GFPfree + T7RNAP_free + mRNA_GFP
r22 = k15 * C4;
%mRNA_GFP + Ribosome_GFP -> mRNA_GFP + Ribosome_GFP + GFP
r23 = k16 * mRNA_GFP * Ribosome_GFP;

%ODEs for each element
dmRNA_T7RNAP = r4 - r14; %Change in T7RNAP mRNA
dT7RNAP = r5 - r13; %Change in T7RNAP 
dC1 = r2 - r4; %Change in transcription initiation complex of T7RNAP
dmRNA_DNAP = r8 - r16; %Change in Phi29 DNAP mRNA
dDNAP = r9 - r15; %Change in Phi29 DNAP
dC2 = r6 - r7 - r8; %Change in transcription initiation complex of Phi29 DNAP
dmRNA_GFP = r22 - r19; %Change in GFP mRNA
dGFP = r23 - r18; %Change in GFP
dC4 = r20 - r21 - r22; %Change in transcription initiation complex of GFP
dRibosome_T7RNAP = 70*dmRNA_T7RNAP; %Change in ribosomes translating T7RNAP mRNA
dRibosome_DNAP = 70*dmRNA_DNAP; %Change in ribosomes translating Phi29 DNAP mRNA
dRibosome_GFP = 70*dmRNA_GFP; %Change in ribosomes translating GFP mRNA
dPlasmid = r12 - r17; %Change in linear plasmid holding GFP
dC3 = r10 - r11 - r12; %Change in replication complex


%%Assigning output variables
dxdt = [dC1; dmRNA_T7RNAP; dRibosome_T7RNAP; dT7RNAP; dC2; dmRNA_DNAP; ...
    dRibosome_DNAP; dDNAP; dC4; dmRNA_GFP; dRibosome_GFP; dGFP; dC3; dPlasmid];
end