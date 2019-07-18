%%
%Started: 11-07-2019
%Last edit: 11-07-2019

%%
clear all
close all

%%Parameters
%Basics
Ecoli_volume = 0.65e-15;
Avogadro = 6.022e23;
tspan = [0 500]; %Time span in minutes
numElements = 2; %Number of elements in the model
names = {'DNAP', 'Plasmid'}; %Names of elements in the model with complexes

%Initial concentrations
C0 = zeros(1,numElements)+0.001; %Initiating vector containing starting concentrations of all elements, most are zero
C0(1) = (10/Avogadro)/Ecoli_volume *10^9; %Starting concentration of DNAP in nM
%C0(1) = 10^4; %Copies of DNAP
C0(2) = (1/Avogadro)/Ecoli_volume *10^9; %Starting concentration of plasmid in nM
%C0(2) = 1; %Copies of plasmid

%%Model Calculation
%ODEsolver that evaluates Initial_Replication_function over time specified by 'tspan'
%and with initial conditions listed in 'C0'.
opts = odeset('Stats', 'on');
disp('ode15s stats for replication model:')
tic, [t, y] = ode15s(@(t,C) Replication_function(C), tspan, C0, opts);
toc
y(:,:) = y(:,:)*10^-9*0.65*10^-15*6.022*10^23; %Convert from concentration to copy number

%%
%Plot concentrations over time
figure('Name','Orthogonal Replication');
for i = 1:size(y,2)
    subplot(1,2,i)
    plot(t, y(:,i), 'LineWidth', 2);
    title(names(i));
    ylabel('Copy number')
    xlabel('Time (min)')
    grid on
end


function dC = Replication_function(C)
%%Variables

%Concentration variables of all elements with changing concentrations
DNAP = C(1); %Phi29 DNA polymerase
Plasmid = C(2); %GFP promoter / plasmid copy number

%Rate constants
k9 = 709; %nM^-1 min^-1
k12 = 0.023; %min^-1


%%Calculation of concentration changes

%Rate laws
%Plasmid_free + 2 DNAP_free -> 2 Plasmid + 2 DNAP
r10 = k9 * Plasmid * DNAP^2;
%Plasmid -> Ø
r17 = k12 * Plasmid;

%Mass balances for each element
dDNAP = 0; %Change in Phi29 DNAP
dPlasmid = r10 - r17; %Change in linear plasmid holding GFP


%%Assigning output variables
dC(1,:) = dDNAP;
dC(2,:) = dPlasmid;
end