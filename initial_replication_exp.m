function dC = Initial_Replication_function(C, IPTG)

%%Constants

%Constant concentrations
Ecoli_volume = 0.65e-15;
Avogadro = 6.022e23;
eRNAP = (1929/Avogadro)/Ecoli_volume *10^9; %Endogenous RNAP concentration in nM
Pr_DNAP = (1/Avogadro)/Ecoli_volume *10^9; %DNAP promoter concentration in nM
%eRNAP = 1929; %Number of endogenous RNAP copies
%Pr_DNAP = 1; %Number of DNAP promoter copies
%Total_Ribosomes = 0.8*45.1*10^3; %Total number of ribosomes per cell * activity (80%)
Total_Ribosomes = (((0.8*45.1*10^3)/Avogadro)/Ecoli_volume)*10^9; %Active ribosome concentration

%Presence of IPTG is evaluated and active T7RNAP promoter (Lac promoter)
%concentration is updated accordingly
if IPTG == false
    %If there is no IPTG added, there are no active T7 promoters
    Pr_T7RNAP = 0;
    Pr_p5 = 0;
    Pr_p6 = 0;
elseif IPTG == true
    %If IPTG is added, the T7 promoters are activated, leading to a
    %concentration of 2.6 nM
    Pr_T7RNAP = (1/Avogadro)/Ecoli_volume *10^9;
    Pr_p5 = (1/Avogadro)/Ecoli_volume *10^9;
    Pr_p6 = (1/Avogadro)/Ecoli_volume *10^9;
    %Pr_T7RNAP = 1;
    %Pr_p5 = 1;
    %Pr_p6 = 1;
end


%%Variables

%Concentration variables of all elements with changing concentrations
C1 = C(1); %Transcription initiation complex of T7RNAP
mRNA_T7RNAP = C(2); %T7RNAP mRNA
Ribosome_T7RNAP = C(3); %Ribosomes translating T7RNAP mRNA
T7RNAP = C(4); %T7 RNA polymerase
C2 = C(5); %Transcription initiation complex of DNAP
mRNA_DNAP = C(6); %Phi29 DNAP mRNA
Ribosome_DNAP = C(7); %Ribosomes translating Phi29 DNAP mRNA
DNAP = C(8); %Phi29 DNA polymerase
C4 = C(9); %Transcription initiation complex of GFP
mRNA_GFP = C(10); %GFP mRNA
Ribosome_GFP = C(11); %Ribosomes translating GFP mRNA
GFP = C(12); %GFP
C3 = C(13); %Replication complex
Plasmid = C(14); %GFP promoter / plasmid copy number
C5 = C(15); %Transcription complex for P5
C6 = C(16); %Transcription complex for P6

%Mole balances to calculate concentration of unoccupied elements when
%applicable
eRNAP_free = eRNAP - C1; %Endogenous RNAP that is not transcribing
T7RNAP_free = T7RNAP - C2; %T7RNAP that is not transcribing
DNAP_free = DNAP - C3; %DNAP that is not replicating
Pr_T7RNAPfree = Pr_T7RNAP - C1; %Unoccupied T7RNAP promoter
Pr_DNAPfree = Pr_DNAP - C2; %Unoccupied DNAP promoter
Pr_GFPfree = Plasmid - C4; %Unoccupied GFP promoter
Pr_p5free = Pr_p5 - C5; %Unoccupied P5 promoter
Pr_p6free = Pr_p6 - C6; %Unoccupied P6 promoter
Plasmid_free = Plasmid - C3; %Plasmid that is not being replicated

%Rate constants
k1 = 3120; %Transcription elongation rate by eRNAP (min^-1)
k2 = 1042; %Transcription initiation rate by endogenous RNAP (min^-1)
k3 = 1200; %Translation rate per ribosome (min^-1)
%k4 = 1.1568; %Transcription elongation rate of T7RNAP by endogenous RNAP (min^-1)
%k5 = 1.33; %Translation rate of T7RNAP (min^-1)
k6 = 4.320; %T7RNAP binding rate to T7 promoter (nM^-1 min^-1)
k7 = 240; %T7RNAP dissociation rate from T7 promoter (min^-1)
k8 = 1.50; %Transcription elongation rate of DNAP by T7RNAP (min^-1)
k9 = 1.08; %Phi29 DNAP binding rate to oriL or oriR (nM^-1 min^-1)
k10 = 200; %Phi29 DNAP dissociation rate from oriL or oriR (min^-1)
%k11 = 6000; %Replication rate by Phi29 DNAP (min^-1)
k11 = 0.46; %Plasmid replication rate by Phi29 DNAP (min^-1)
k12 = 0.023; %Protein degradation rate (min^-1)
k13 = 0.139; %mRNA degradation rate (min^-1)
%k14 = 2.09; %Translation rate of Phi29 DNAP (min^-1)
k15 = 3.81; %Transcription elongation rate of GFP by T7RNAP (min^-1)
%k16 = 5.31; %Translation rate of GFP (min^-1)
k17 = 2580; %Transcription elongation by T7RNAP (min^-1)

%Gene lengths(bp)
L_T7RNAP = 2697;
L_DNAP = 1719;
L_GFP = 678;
L_P5 = 375;
L_P6 = 315;
L_Plasmid = 10000;

%%Calculation of concentration changes

%Rate laws

%Pr_T7RNAPfree + eRNAP_free <--> C1
r2 = k2 * Pr_T7RNAPfree * eRNAP_free;
%C1 -> Pr_T7RNAPfree + eRNAP_free + mRNA_T7RNAP
r4 = k1/L_T7RNAP * C1;
%mRNA_T7RNAP + Ribosome_T7RNAP -> mRNA_T7RNAP + Ribosome_T7RNAP + T7RNAP
r5 = k3/(L_T7RNAP/3) * mRNA_T7RNAP * Ribosome_T7RNAP;
%Pr_DNAPfree + T7RNAP_free <--> C2
r6 = k6 * Pr_DNAPfree * T7RNAP_free;
r7 = k7 * C2;
%C2 -> Pr_DNAPfree + T7RNAP_free + mRNA_DNAP
r8 = k8 * C2;
%mRNA_DNAP + Ribosome_DNAP -> mRNA_DNAP + Ribosome_DNAP + DNAP
r9 = k3/(L_DNAP/3) * mRNA_DNAP * Ribosome_DNAP;
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
r23 = k3/(L_GFP/3) * mRNA_GFP * Ribosome_GFP;
%Pr_p5free + T7RNAP_free <--> C5
r24 = k6*Pr_p5free*T7RNAP_free;
r25 = k7*C5;
%Pr_p6free + T7RNAP_free <--> C6
r26 = k6*Pr_p6free*T7RNAP_free;
r27 = k7*C6;

%Mass balances for each element
dmRNA_T7RNAP = r4 - r14; %Change in T7RNAP mRNA
dT7RNAP = r5 - r13; %Change in T7RNAP
dC1 = r2 - r4; %Change in transcription initiation complex of T7RNAP
dmRNA_DNAP = r8 - r16; %Change in Phi29 DNAP mRNA
dDNAP = r9 - r15; %Change in Phi29 DNAP
dC2 = r6 - r7 - r8; %Change in transcription initiation complex of Phi29 DNAP
dmRNA_GFP = r22 - r19; %Change in GFP mRNA
dGFP = r23 - r18; %Change in GFP
dC4 = r20 - r21 - r22; %Change in transcription initiation complex of GFP
if Ribosome_T7RNAP + Ribosome_DNAP + Ribosome_GFP ...
        + 3*dmRNA_T7RNAP + 3*dmRNA_DNAP + 3*dmRNA_GFP <= Total_Ribosomes
    dRibosome_T7RNAP = 3*dmRNA_T7RNAP; %Change in ribosomes translating T7RNAP mRNA
    dRibosome_DNAP = 3*dmRNA_DNAP; %Change in ribosomes translating Phi29 DNAP mRNA
    dRibosome_GFP = 3*dmRNA_GFP; %Change in ribosomes translating GFP mRNA
else
    dRibosome_T7RNAP = 0; %Change in ribosomes translating T7RNAP mRNA
    dRibosome_DNAP = 0; %Change in ribosomes translating Phi29 DNAP mRNA
    dRibosome_GFP = 0; %Change in ribosomes translating GFP mRNA
end
%dPlasmid = r12 - r17; %Change in linear plasmid holding GFP
dPlasmid = 0;
%dC3 = r10 - r11 - r12; %Change in replication complex
dC3 = 0;
dC5 = r24 - r25;
dC6 = r26 - r27;


%%Assigning output variables
dC(1,:) = dC1;
dC(2,:) = dmRNA_T7RNAP;
dC(3,:) = dRibosome_T7RNAP;
dC(4,:) = dT7RNAP;
dC(5,:) = dC2;
dC(6,:) = dmRNA_DNAP;
dC(7,:) = dRibosome_DNAP;
dC(8,:) = dDNAP;
dC(9,:) = dC4;
dC(10,:) = dmRNA_GFP;
dC(11,:) = dRibosome_GFP;
dC(12,:) = dGFP;
dC(13,:) = dC3;
dC(14,:) = dPlasmid;
dC(15,:) = dC5;
dC(16,:) = dC6;
end