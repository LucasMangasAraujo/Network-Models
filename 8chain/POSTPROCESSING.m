%% POSTPROCESSING CODE: PLOTING NOMINAL STRESS VERSUS STRETCH CURVES
clc; close all; close all;

%% OPEN FILES FROM PYTHON SIMULATIONS 

%  Prompt for model nale
MODEL = input('Type model:\n','s');

% Get Files ID
EXT = '.txt'; % File Extension

PYUNI = strcat('UNI',MODEL,EXT); % Uniaxial
PYBI = strcat('BI',MODEL,EXT);  % Biaxial
PYSHEAR = strcat('SHEAR',MODEL,EXT); % Shear

UNI = load(PYUNI);
BI = load(PYBI);
SHEAR = load(PYSHEAR);

%% Experimental Data
% Data
DATAUNI=load('TRELOARUNI.txt');
DATABI=load('TRELOARBI.txt');
DATASHEAR=load('TRELOARSHEAR.txt');

%% Boyce & Arruda Simulations
%  Files names
BOYCEUNI_ID = strcat('BOYCEUNI_', MODEL, EXT);
BOYCEBI_ID = strcat('BOYCEBI_', MODEL, EXT);
BOYCESHEAR_ID = strcat('BOYCESHEAR_', MODEL, EXT);

%  Boyce and Arruda Simulations
BOYCEUNI = load(BOYCEUNI_ID);
BOYCEBI = load(BOYCEBI_ID);
BOYCESHEAR = load(BOYCESHEAR_ID);


%% ASSEMLE UNIAXIAL DATA

% Experiment
UNISTRETCHEXP=DATAUNI(:,1); UNINOMSTRESEXP=DATAUNI(:,2);
BISTRETCHEXP=DATABI(:,1); BINOMSTRESEXP=DATABI(:,2);
SHEARSTRETCHEXP=DATASHEAR(:,1); SHEARNOMSTRESEXP=DATASHEAR(:,2);

%  Boyce's simulation
UNISTRETCHBOYCE=BOYCEUNI(:,1); UNINOMSTRESBOYCE=BOYCEUNI(:,2);
BISTRETCHBOYCE=BOYCEBI(:,1);   BINOMSTRESBOYCE=BOYCEBI(:,2);
SHEARSTRETCHBOYCE=BOYCESHEAR(:,1); SHEARNOMSTRESBOYCE=BOYCESHEAR(:,2);

% Python
UNISTRETCH = UNI(:,1); UNINOMSTRES = UNI(:,2);
BISTRETCH = BI(:,1);   BINOMSTRES = BI(:,2);
SHEARSTRETCH = SHEAR(:,1); SHEARNOMSTRES = SHEAR(:,2);

% Plot Results From Gaussian Modeol
figure(1)
plot(UNISTRETCHEXP,UNINOMSTRESEXP,'ko')
hold on
plot(BISTRETCHEXP,BINOMSTRESEXP,'kd')
plot(SHEARSTRETCHEXP,SHEARNOMSTRESEXP,'k+')
% Plot Boyce's Simulation
plot(UNISTRETCHBOYCE,UNINOMSTRESBOYCE,'bo')
plot(BISTRETCHBOYCE,BINOMSTRESBOYCE,'gd')
plot(SHEARSTRETCHBOYCE,SHEARNOMSTRESBOYCE,'r+')
% Our implemantation
plot(UNISTRETCH,UNINOMSTRES,'b-')
plot(BISTRETCH,BINOMSTRES,'g-')
plot(SHEARSTRETCH,SHEARNOMSTRES,'r-')

axis square
xlabel('Stretch','FontSize',16,'Interpreter','Latex')
ylabel('Nominal Stress$\left(MPa\right)$','FontSize',16,'Interpreter','Latex')
ylim([0 7])
if ( strcmp(MODEL,'GAUSS') == 1 ) 
    legend('Treloar Uniaxial','Treloar Biaxial','Treloar Pure Shear',...
        'Boyce and Arruda Simulation Uniaxial Gauss','Boyce and Arruda Simulation Biaxial','Boyce and Arruda Simulation Shear',...
        'Implementation Uniaxial Gauss','Implementation Biaxial Gauss')
elseif(strcmp(MODEL,'3CHAIN') == 1)
        legend('Treloar Uniaxial','Treloar Biaxial','Treloar Pure Shear',...
        'Boyce and Arruda Simulation Uniaxial','Boyce and Arruda Simulation Biaxial','Boyce and Arruda Simulation Shear',...
        'Implementation Uniaxial 3-chain','Implementation Biaxial 3-chain','Implementation Shear 3-chain')
elseif(strcmp(MODEL,'8CHAIN') == 1)
        legend('Treloar Uniaxial','Treloar Biaxial','Treloar Pure Shear',...
        'Boyce and Arruda Simulation Uniaxial','Boyce and Arruda Simulation Biaxial','Boyce and Arruda Simulation Shear',...
        'Implementation Uniaxial 8-chain','Implementation Biaxial 8-chain','Implementation Shear 8-chain')
end

set(legend, 'Interpreter','Latex')