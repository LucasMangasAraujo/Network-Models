%% VALIDATION: Code for vizualition of the full-network result
clear all; close all; clc;

%% Open File and Ploat

%  Open Files
VSIM = load('VERRON.txt');
SIM = load('UNIPOLYFULL_DIRAC.txt');
UNIFORM = load('UNIPOLYFULL_UNIFORM.txt');
EXPO = load('UNIPOLYFULL_EXPONENTIAL.txt');

% Associate 
L = VSIM(:,1); NOMS = VSIM(:,2);
LSIM = SIM(:,1); NOMSSIM = SIM(:,2);
LUNIFORM = UNIFORM(:,1); NOMSUNIFORM = UNIFORM(:,2);
LEXPO = EXPO(:,1); NOMSEXPO = EXPO(:,2);

%% Open data related to Itskov 

% Axial:

% Itskov data
DATA = load('COMPRESSUNI_DATA.txt');
MONO = load('COMPRESSUNIFULL.txt');
EXPO2 = load('COMPRESSUNIPOLYFULL_EXPONENTIAL.txt');
ITS = load('ITSKOV.txt');

%Associate
LEXP = DATA(:,1); NOMSEXP = DATA(:,2);
LITS = ITS(:,1);  NOMSITS = ITS(:,2); 
LMONO = MONO(:,1); NOMSMONO = MONO(:,2);
LEXPO2 = EXPO2(:,1); NOMSEXPO2 = EXPO2(:,2);
% ========================================================================

% Shear
DATASHEAR = load('COMPSHEAR_DATA.txt');
ITSHEAR = load('ITSKOV_SHEAR.txt');
MONOSHEAR = load('COMPSHEARFULL.txt');
EXPOSHEAR = load('COMPSHEARPOLYFULL_EXPONENTIAL.txt');
% EXPOSHEAR =

% Associate
LEXPS = DATASHEAR(:,1); NOMSEXPS = DATASHEAR(:,2);
LITSHEAR = ITSHEAR(:,1);  NOMSITSHEAR = ITSHEAR(:,2);
LMONOS = MONOSHEAR(:,1); NOMSMONOS = MONOSHEAR(:,2);
LEXPOS = EXPOSHEAR(:,1); NOMSEXPOS = EXPOSHEAR(:,2);

%% Plot

figure(1)
hold on
plot(L, NOMS,'bo')
plot(LSIM, NOMSSIM,'LineWidth',1.3);
plot(LUNIFORM,NOMSUNIFORM,'LineWidth',1.3)
plot(LEXPO,NOMSEXPO,'LineWidth',1.3)
axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nomimal Stress $\left[ \mathrm{MPa} \right] $','FontSize',16,'Interpreter','Latex')
ylim([0,2.5])
legend('Verron Result','Monodisperse','Uniform','Exponential')
hold off

figure(2)
hold on
plot(LEXP,abs(NOMSEXP),'bo')
plot(LITS,abs(NOMSITS),'--')
plot(LMONO,abs(NOMSMONO),'LineWidth',1.3)
plot(LEXPO2,abs(NOMSEXPO2),'LineWidth',1.3)
axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nomimal Stress $\left[ \mathrm{MPa} \right] $','FontSize',16,'Interpreter','Latex')
ylim([0,30])
set(gca,'XDir','reverse'); % How to inverte the ticks of the plot
legend('Uniaxial Compression Data','Itskov Model','Monodisperse','Exponential')
hold off


figure(3)
hold on
plot(LEXPS,abs(NOMSEXPS),'bv')
plot(LITSHEAR,abs(NOMSITSHEAR),'--')
plot(LMONOS,abs(NOMSMONOS),'LineWidth',1.3)
plot(LEXPOS,abs(NOMSEXPOS),'LineWidth',1.3)
axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nomimal Stress $\left[ \mathrm{MPa} \right] $','FontSize',16,'Interpreter','Latex')
ylim([0,30])
set(gca,'XDir','reverse'); % How to inverte the ticks of the plot
legend('Plane Strain Compression Data','Itskov Model','Monodisperse','Exponential')
hold off



