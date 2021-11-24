%% CODF: CODE associated with the chain oriendtation distribuition function
clc; close all; close all;

%% Variables Declaration

% Stretch
L1 = 2.5; L2 = sqrt(L1) ; L3 = sqrt(L1);

% Uniform Distribution
C0= 1.0/(4.*pi);

% Angles
THETA = [0:pi/100:pi];
VARPHI = [0:pi/100:2*pi]; NPHI =  length(VARPHI);

% Mesh
[MTHETA,MPHI] = meshgrid(THETA,VARPHI);


% Spherical Coordinates 
M1 = sin(MTHETA).*cos(MPHI); 
M2 = sin(MTHETA).*sin(MPHI);
M3 = cos(MTHETA);

% Chain stretch
L = ( ((M1./L1).^2) + ((M2./L2).^2) + ((M3./L3).^2)   ).^(-0.5);

% Distrubution
C = C0*L.^3.0; 


%% Normalized COFD
surfc (MPHI./pi, MTHETA./pi, C./C0)
hold on
%  az = 45;
%  el = 45;
%  view(az, el);
xlabel('$\frac{x}{a}$','FontSize',22,'Interpreter','Latex')
ylabel('$\frac{y}{a}$','FontSize',22,'Interpreter','Latex')
zlabel('$\frac{\sigma_x^n}{p_0}$','FontSize',22,'Interpreter','Latex')
hold off