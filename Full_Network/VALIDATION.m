%% VALIDATION: Code for vizualition of the full-network result
clear all; close all; clc;

%% Open File and Ploat

%  Open Files
WU = load('UNIWU.txt');
SIM = load('UNIFULL.txt');
EIGHT_CH = load('UNI8CHAIN.txt');
THREE_CH = load('UNI3CHAIN.txt');

% Associate 
LWU = WU(:,1); NOMSWU = WU(:,2);
L8 = EIGHT_CH(:,1); NOMS8 = EIGHT_CH(:,2);
L3 = THREE_CH(:,1); NOMS3 = THREE_CH(:,2);
LSIM = SIM(:,1); NOMSSIM = SIM(:,2);

%% Plot

figure(1)
hold on
plot(LWU, NOMSWU,'bo',LSIM, NOMSSIM,'r-',L8,NOMS8,'k-',L3,NOMS3,'g-')
axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nomimal Stress $\left[ \mathrm{MPa} \right] $','FontSize',16,'Interpreter','Latex')
ylim([0,6])
legend('Wu Result','Implementation','8-Chain Implemantation','3-Chain Implemantation')
hold off