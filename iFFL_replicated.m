%iFFL system replicated
%Started: 18-07-2019

%%
%Figure 1c
clear all
close all

%Parameters
copyN = 1;
R = 1:1:100;
nvector = 0.5:0.25:1.5;

figure(1)
for ii = 1:length(nvector)
    n = nvector(ii);
    G = copyN./(R.^n);
    ylim([10^-2, 10^0]);
    loglog(R, G)
    legendentries(ii) = {['n = ', num2str(nvector(ii))]};
    hold on
end
title('Figure 1c')
legend(legendentries)
xlabel('Repressor')
ylabel('GOI expression')

%%
%Figure 1d
clear all

aG = 0.685*10^-3;
bG = 5.31;
aR = 0.685*10^-3;
bR = 5.31;
yR = log(2)/30;
ym = log(2)/5;
yG = log(2)/30;
Kd = 240/4320;
alpha = (bG*(aG/Kd))/(yG*ym);
beta = (bR*aR)/(yR*ym);

copyN = 1:1:100;
nvector = 0.5:0.25:1.5;

figure(2)
for ii = 1:length(nvector)
    n = nvector(ii);
    scaling = (aG*bG*yR^n*ym^n)/(aR^n*bR^n*yG*ym);
    G_14 = scaling.*copyN.^(1-n);
    G_15 = copyN.^(1-n);
    G_16 = alpha.*((copyN.*Kd^n)./(Kd^n+beta^n.*copyN.^n));
    %ylim([10^-1, 10^1])
    
    subplot(1,3,1)
    loglog(copyN, G_15)
    title('G = c^{1-n}')
    xlabel('Copy number')
    ylabel('GOI expression')
    hold on
    
    subplot(1,3,2)
    loglog(copyN, G_14)
    title('G ~ c^{1-n}')
    hold on
    
    subplot(1,3,3)
    loglog(copyN, G_16)
    title('G = a(cKd/(Kd+Bc))')
    hold on
    
    legendentries(ii) = {['n = ', num2str(nvector(ii))]};
end
legend(legendentries)

%%
%Figure 1e

clear all

E = 0.01:0.01:1;

figure(3)
S = E./(E+1);
ylim([10^-2, 10^0])
loglog(E, S)

title('Figure 1e')
xlabel('Stabilization error')
ylabel('Stabilized promoter strength')