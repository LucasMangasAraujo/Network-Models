%% POSTPROCESSING CODE: PLOTING NOMINAL STRESS VERSUS STRETCH CURVES
clc; close all; close all;

%% OPEN FILES FROM PYTHON SIMULATIONS 

%  Prompt for model nale
MODEL = input('Type model:\n','s');

% Get Files ID
EXT = '.txt'; % File Extension

PYUNI = strcat('UNI',MODEL,EXT); % Uniaxial
PYBI = strcat('BI',MODEL,EXT);  % Biaxial
% PYSHEAR = strcat('SHEAR',MODEL,EXT); % Shear

UNI = load(PYUNI);
BI = load(PYBI);

%% Experimental Data
% Data
DATAUNI = load('UNIEXP.txt');
DATABI = load('BIEXP.txt');
% DATASHEAR=load('TRELOARSHEAR.txt');

%% Boyce & Arruda Simulations
%  Files names

%  Boyce and Arruda Simulations
UNIWU = load('UNIWU2.txt');
WUBI = load('BIWU.txt');


%% ASSEMLE UNIAXIAL DATA

% Experiment
UNISTRETCHEXP=DATAUNI(:,1); UNINOMSTRESEXP=DATAUNI(:,2);
BISTRETCHEXP=DATABI(:,1); BINOMSTRESEXP=DATABI(:,2);

%  Boyce's simulation
UNISTRETCHWU = UNIWU(:,1); UNINOMSTRESWU = UNIWU(:,2);
BISTRETCHWU=WUBI(:,1);   BINOMSTRESBOYCE=WUBI(:,2);


% Python
UNISTRETCH = UNI(:,1); UNINOMSTRES = UNI(:,2);
BISTRETCH = BI(:,1);   BINOMSTRES = BI(:,2);

% Plot Results From Gaussian Modeol
figure(1)
plot(UNISTRETCHEXP,UNINOMSTRESEXP,'ko')
hold on
plot(BISTRETCHEXP,BINOMSTRESEXP,'kv')
% Plot Boyce's Simulation
plot(UNISTRETCHWU,UNINOMSTRESWU,'bo')
plot(BISTRETCHWU,BINOMSTRESBOYCE,'go')

% Our implemantation
plot(UNISTRETCH,UNINOMSTRES,'b-')
plot(BISTRETCH,BINOMSTRES,'g-')

axis square
xlabel('$\lambda$','FontSize',16,'Interpreter','Latex')
ylabel('Nominal Stress$\left(MPa\right)$','FontSize',16,'Interpreter','Latex')
ylim([0 3])
legend('Experinemt Uniaxial','Experinemt Biaxial',...
        'Wu Uniaxial Simulation','Wu Simulation Biaxial',...
        'Python 21-Points Scheme Uniaxial','Python 21-Points Scheme Biaxial')


set(legend, 'Interpreter','Latex')