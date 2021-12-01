%% VALIDATION: Code for vizualition of the full-network result
clear all; close all; clc;

%% Open File and Ploat

%  Open Files
VSIM = load('VERRON.txt');
MONO = load('UNIPOLYFULL.txt');

% Associate 
L = VSIM(:,1); NOMS = VSIM(:,2);
LSIM = MONO(:,1); NOMSSIM = MONO(:,2);

%% Plot

figure(1)
hold on
plot(L, NOMS,'bo',LSIM, NOMSSIM)
axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nomimal Stress $\left[ \mathrm{MPa} \right] $','FontSize',16,'Interpreter','Latex')
ylim([0,2.5])
legend('Verron Result','Monodisperse')
hold off