%%
clear all; close all; clc;
format long g
% Parameters SIP:
Freq_Center = 193.54;
CouplCoeff_k1 = 0.49832327234602;
Radius = 10e-6;
Alpha = 1.15240339832846;
Alpha_2 = 0.980630321583591;

EffRefrIndex = 2.362;
SpeedLight = 3e8; % m/s
Imag = -1e-3; % Imaginary part of the Effective refractive index. Negative if gain.
n = EffRefrIndex + j*Imag;
c = SpeedLight;
k0 = 2*pi*Freq_Center*1e12/c;
R = Radius;
k = CouplCoeff_k1;
t = sqrt(1-k^2);

% 1. Find Ksd
A = -t^4/k^4 +2*(t^6-2*t^4+t^2)/k^6;

Ksd = (1/2)*acos((1/6)*(-A-9));

Ksd_pi = Ksd/pi;
% Checked. Works well.

% 2. Find \beta = \alpha - \alpha'
B = (k^2/t^2)*3*cos(Ksd);

beta = (1/(2*R*n*k0))*acos(B);
beta = 0.171773076744869; 
B_2 = cos((2*R*n*k0)*beta);
Ksd_2 = acos(B_2*t^2/(3*k^2));
Ksd_2_pi = Ksd_2/pi;
% 3. Find \gamma = \alpha + \alpha'
C = (k^4*(cos(3*Ksd)+9*cos(Ksd)));

gamma = ((1/(2*k0*n*R))*acos(C-2*(t^2-t^4)*B))-pi;

% 4. Find \alpha & \alpha':

Alpha_a = 0.5*(beta + gamma)
Alpha_b = 0.5*(gamma-beta)


%% B & betas
B = 0.0942757341393594;
beta = 0.00771008250509886;
beta_2 = 0.171773076744869;

beta_test = (1/(2*R*n*k0))*(-i*log(B+i*sqrt(1-B^2)))
beta_test_2 = (1/(2*R*n*k0))*(-i*log(B+i*sqrt(1-B^2)))

Delta_beta = beta_2 - beta;

ShouldBeZero = exp(i*2*R*n*k0*Delta_beta) + exp(i*3*R*n*k0*Delta_beta)


