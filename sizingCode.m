%% AME 415 Turbine Sizing Code
% by Adam 
%% Cleaning
clear; clc;
%% Constants
g0 = 9.807; % [m/s^2]

m2ft = 3.2801;
hp2w = 745.699872;

%% Inputs
T01 = 500;  % [K]
P01 = 13842074;   % [Pa]
P3 = 11011992;    % [Pa]
mdot = 13.61;   % [kg/s]

c_p = 15636.11;   % [J/kg/K]
gam = 1.3994;
R = 4123.311;
z = 1;

power = 6215;   % [HP]
N = 50000;   % [RPM]

% T01 = 753;  % [K]
% P01 = 182385;   % [Pa]
% P3 = 101325;    % [Pa]
% mdot = 9;   % [kg/s]
% 
% c_p = 1156.9;   % [J/kg/K]
% gam = 1.33;
% R = c_p*(1 - gam^-1);
% z = 1;
%
% power = 1200;   % [HP]
% N = 8000;   % [RPM]

eta_TT = 0.9;

Nb_stator = 80; Nb_rotor = 100;
alpha1 = 90;

BL = 0.9; Zwiff = 0.8;

phi = 0.4;  psi = 1.6;
Rc = 0.4;

inletParamString = strcat('------------------\n Inlet parameters \n------------------\n', ...
    'T01 = %3.0fK, P01 = %1.3E Pa, P3 = %1.3E Pa\n', ...
    'Cp = %1.3E J/kg/K, gamma = %1.3f, R = %4.3f J/kg/K\n', ...
    'Required power = %d W, blade speed = %d RPM \n', ...
    'Reaction Rc = %1.1f, eta_TT = %1.1f, flow coefficient phi = %1.1f \n');
fprintf(inletParamString, T01, P01, P3, c_p, gam, R, power, N, Rc, eta_TT, phi);

%% Calculations

UoverC0 = sqrt(1 / (phi^2 - 4*(Rc - 1)/eta_TT));

T3is = T01*(P3/P01)^((gam-1)/gam);
C0 = sqrt(2*c_p*(T01 - T3is));
U = C0*UoverC0;

Ca = U*phi;

dh0is = (C0^2 - Ca^2)/2;
W = eta_TT*dh0is;
power_est = mdot*W/hp2w;

alpha2 = atand(U*Ca/W);

C2 = Ca/sind(alpha2);
C2u = Ca/tand(alpha2);
C2a = Ca;
W2u = C2u - U;
W2a = Ca;
W2 = sqrt(W2u^2 + W2a^2);
alpha2p = atand(W2a/W2u);

W3a = Ca;
W3u = -U;
W3 = sqrt(W3a^2 + W3u^2);
C3a = Ca;
C3u = 0;
C3 = Ca;
alpha3 = 90;
alpha3p = atand(W3a/W3u);

T2 = T01 - (C2^2 / 2 / c_p);
a2 = sqrt(gam*R*T2);
M2 = C2/a2;

Mw2 = W2/a2;
Tw2 = T2*(1 + 0.5*(gam-1)*Mw2^2);
T3 = Tw2 - (W3^2 / 2 / c_p);
a3 = sqrt(gam*R*T3);
Mw3 = W3/a3;

REF = [gam, P01, M2, Mw2, Mw3, P3];
lossguess = 0.05;
f = @(y) lossCalcFunc(y, REF);
[loss,fval] = fzero(f, lossguess);
% loss = 0.02;

P02 = P01*(1 - loss);
P2 = P02/(1 + 0.5*(gam-1)*M2^2)^(gam/(gam-1));
Pw2 = P2*(1 + 0.5*(gam-1)*Mw2^2)^(gam/(gam-1));
Pw3 = Pw2*(1 - loss);

rho2 = P2/z/R/T2;
rho3 = P3/z/R/T3;

A2 = mdot/rho2/C2a;
A3 = mdot/rho3/C3a;

rmean = 30*U/pi/N;

L2 = A2/2/pi/rmean/BL;
L3 = A3/2/pi/rmean/BL;

%% Stator

sc0_stator = 0.427 + alpha2/58 - (alpha2/93)^2;
sBz_stator = Zwiff/(2*(sind(alpha2)^2)*(cotd(alpha2) - cotd(alpha1)));

s_stator = 2*pi*rmean/Nb_stator;

Bz_stator = s_stator/sBz_stator; c_stator = s_stator/sc0_stator;
betaS_stator = asind(Bz_stator/c_stator);

%% Rotor
sc0_rotor = 0.427 + abs(alpha3p)/58 - (alpha3p/93)^2;
sc1_rotor = 0.224 + (1.575 - abs(alpha3p)/90)*(abs(alpha3p)/90);

zeta = (90 - alpha2p)/(90 - abs(alpha3p));

sc_opt = sc0_rotor + (sc1_rotor - sc0_rotor)*abs(zeta)*zeta;

sBz_rotor = Zwiff/(2*(sind(alpha3p)^2)*(cotd(alpha2p) - cotd(alpha3p)));

s_rotor = 2*pi*rmean/Nb_rotor;

Bz_rotor = s_rotor/sBz_rotor; c_rotor = s_rotor/sc_opt;
betaS_rotor = asind(Bz_rotor/c_rotor);

%% 3D Design

r_tip_stator = rmean + 0.5*L2;  r_tip_rotor = rmean + 0.5*L3;
r_hub_stator = rmean - 0.5*L2;  r_hub_rotor = rmean - 0.5*L3;

U_tip = r_tip_stator*pi*N/30;
U_hub = r_hub_stator*pi*N/30;

C2u_tip = rmean*C2u/r_tip_stator;   C2a_tip = Ca;
C2u_hub = rmean*C2u/r_hub_stator;   C2a_hub = Ca;
C2_tip = sqrt(C2u_tip^2 + C2a_tip^2);
C2_hub = sqrt(C2u_hub^2 + C2a_hub^2);

W2u_tip = C2u_tip - U_tip;  W2a_tip = C2a_tip;
W2_tip = sqrt(W2u_tip^2 + W2a_tip^2);
W2u_hub = C2u_hub - U_hub;  W2a_hub = C2a_hub;
W2_hub = sqrt(W2u_hub^2 + W2a_hub^2);

alpha2_tip = asind(C2a_tip/C2_tip);  alpha2p_tip = atand(W2a_tip/W2u_tip);
alpha2_hub = asind(C2a_hub/C2_hub);  alpha2p_hub = atand(W2a_hub/W2u_hub);

C3u_tip = 0;    C3a_tip = Ca;
C3u_hub = 0;    C3a_hub = Ca;
C3_tip = sqrt(C3u_tip^2 + C3a_tip^2);
C3_hub = sqrt(C3u_hub^2 + C3a_hub^2);

W3u_tip = C3u_tip - U_tip;  W3a_tip = C3a_tip;
W3_tip = sqrt(W3u_tip^2 + W3a_tip^2);
W3u_hub = C3u_hub - U_hub;  W3a_hub = C3a_hub;
W3_hub = sqrt(W3u_hub^2 + W3a_hub^2);

alpha3_tip = asind(C3a_tip/C3_tip);  alpha3p_tip = atand(W3a_tip/W3u_tip);
alpha3_hub = asind(C3a_hub/C3_hub);  alpha3p_hub = atand(W3a_hub/W3u_hub);

tipVals = [C2a_tip, C2u_tip, C2_tip, W2a_tip, W2u_tip, W2_tip, ...
                alpha2_tip, alpha2p_tip, C3a_tip, C3u_tip, C3_tip, ...
                W3a_tip, W3u_tip, W3_tip, alpha3_tip, alpha3p_tip];
hubVals = [C2a_hub, C2u_hub, C2_hub, W2a_hub, W2u_hub, W2_hub, ...
                alpha2_hub, alpha2p_hub, C3a_hub, C3u_hub, C3_hub, ...
                W3a_hub, W3u_hub, W3_hub, alpha3_hub, alpha3p_hub];             

meanVals = [C2a, C2u, C2, W2a, W2u, W2, alpha2, alpha2p, ...
                C3a, C3u, C3, W3a, W3u, W3, alpha3, alpha3p];
outputVals = [hubVals; meanVals; tipVals];

printHeader = {'Hub'; 'Mean'; 'Tip'};

outputParameterString = strcat('C2a = %3.3f m/s, C2u = %3.3f m/s, C2 = %3.3f m/s \n', ...
    'W2a = %3.3f m/s, W2u = %3.3f m/s, W2 = %3.3f m/s \n', ...
    'alpha2 = %2.2fº, alpha2prime = %2.2fº \n', ...
    'C3a = %3.3f m/s, C3u = %3.3f m/s, C3 = %3.3f m/s \n', ...
    'W3a = %3.3f m/s, W3u = %3.3f m/s, W3 = %3.3f m/s \n', ...
    'alpha3 = %2.2fº, alpha3prime = %2.2fº \n');

for ii = 1:3
    headerStr = strcat('------------------\n',printHeader{ii}, ...
        {' Values\n'}, '------------------\n');
    fprintf(headerStr{:});
    fprintf(outputParameterString, outputVals(ii,:));
end