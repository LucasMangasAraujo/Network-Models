%% POSTPROCESSING CODE: PLOTING NOMINAL STRESS VERSUS STRETCH CURVES
clc; close all; close all;

%% OPEN FILES
% Data and Simulation from Boyce

% Data
DATAUNI=load('TRELOARUNI.txt');
DATABI=load('TRELOARBI.txt');
DATASHEAR=load('TRELOARSHEAR.txt');

%  Boyce and Arruda Simulations
BOYCEUNIGAUSS=load('BOYCEUNI_GAUSS.txt');
BOYCEBIGAUSS=load('BOYCEBI_GAUSS.txt');
BOYCESHEARGAUSS=load('BOYCESHEAR_GAUSS.txt');


% Simulations from Pytho
UNIGAUUS=load('UNIGAUSS.txt');
BIGAUSS=load('BIGAUSS.txt');
SHEARGAUSS=load('SHEARGAUSS.txt');

%% ASSEMLE UNIAXIAL DATA

% Experiment
UNISTRETCHEXP=DATAUNI(:,1); UNINOMSTRESEXP=DATAUNI(:,2);
BISTRETCHEXP=DATABI(:,1); BINOMSTRESEXP=DATABI(:,2);
SHEARSTRETCHEXP=DATASHEAR(:,1); SHEARNOMSTRESEXP=DATASHEAR(:,2);

%  Boyce's simulation
UNISTRETCHBOYCE=BOYCEUNIGAUSS(:,1); UNINOMSTRESBOYCE=BOYCEUNIGAUSS(:,2);
BISTRETCHBOYCE=BOYCEBIGAUSS(:,1);   BINOMSTRESBOYCE=BOYCEBIGAUSS(:,2);
SHEARSTRETCHBOYCE=BOYCESHEARGAUSS(:,1); SHEARNOMSTRESBOYCE=BOYCESHEARGAUSS(:,2);

% Python
UNISTRETCHGAUSS=UNIGAUUS(:,1); UNINOMSTRESGAUSS=UNIGAUUS(:,2);
BISTRETCHGAUSS=BIGAUSS(:,1);   BINOMSTRESGAUSS=BIGAUSS(:,2);
SHEARSTRETCHGAUSS=SHEARGAUSS(:,1); SHEARNOMSTRESGAUSS=SHEARGAUSS(:,2);

% Plot Results From Gaussian Modeol
figure(1)
plot(UNISTRETCHEXP,UNINOMSTRESEXP,'ko')
hold on
plot(BISTRETCHEXP,BINOMSTRESEXP,'kd')
plot(SHEARSTRETCHEXP,SHEARNOMSTRESEXP,'k+')
% Plot Boyce's Simulation
plot(UNISTRETCHBOYCE,UNINOMSTRESBOYCE,'ro')
plot(BISTRETCHBOYCE,BINOMSTRESBOYCE,'gd')
plot(SHEARSTRETCHBOYCE,SHEARNOMSTRESBOYCE,'c+')
% Our implemantation
plot(UNISTRETCHGAUSS,UNINOMSTRESGAUSS,'b.')
plot(BISTRETCHGAUSS,BINOMSTRESGAUSS,'g.')
plot(SHEARSTRETCHGAUSS,SHEARNOMSTRESGAUSS,'c.')

axis square
xlabel('Stretch','FontSize',16,'Interpreter','Latex')
ylabel('Nominal Stress$\left(MPa\right)$','FontSize',16,'Interpreter','Latex')
ylim([0 Inf])
legend('Treloar Uniaxial','Treloar Biaxial','Treloar Pure Shear',...
    'Boyce and Arruda Simulation Uniaxial Gauss','Boyce and Arruda Simulation Biaxial','Boyce and Arruda Simulation Shear',...
    'Implementation Uniaxial Gauss','Implementation Biaxial Gauss')
set(legend, 'Interpreter','Latex')