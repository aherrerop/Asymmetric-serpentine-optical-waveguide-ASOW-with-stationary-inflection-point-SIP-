%% Dispersion Diagram Analytic Expression Using Symbolic Maths
clear all; close all; clc
%% Old Notation %%%%%%%%%%%%%%%%

% Parameters
EffRefrIndex = 2.362; % Effective refractive index
Radius = 5e-6; % Radius of the "ring"
AngleSecB = 58.5*pi/180; % Angle of section B from the center of the "ring"
Freq = 2*pi*189e12;

LengthSectionB = 2*AngleSecB*Radius; % Length of section B
TotalLengthUnitCell = 2*pi*Radius + 2*LengthSectionB; % Total length of the unit cell
LightSpeed_c= 3e8; % Speed of light
Wavenumber_Vacuum = Freq/LightSpeed_c; % Wavenumber in vacuum
 
PhaseSectionA = Wavenumber_Vacuum*EffRefrIndex*pi*Radius/2; % Phase accumulated in the section A
syms PhaseSectionA
syms PhaseSectionB 
PhaseSectionB = Wavenumber_Vacuum*EffRefrIndex*LengthSectionB; % Phase accumulated in the section B
syms CouplCoeff_k1 
CouplCoeff_k1 = sqrt(0.24);
syms CouplCoeff_k2 
CouplCoeff_k2 = CouplCoeff_k1;
% Tp1: Transfer Matrix from Z0+ to Zc- (from begining of the unit cell to first coupling point)
% Tp2: Transfer Matrix from Zc+ to Zc'- (from first to second coupling
% point)
% Tp3: Transfer Matrix from Zc'+ to Z1- (from second coupling point to
% end of the unit cell)
syms Tp1(PhaseSectionA, PhaseSectionB)

Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(1j*PhaseSectionA);
Tp1(2,2)= exp(1j*PhaseSectionB);
Tp1(3,3)= exp(1j*PhaseSectionA);
Tp1(4,4)= exp(-1j*PhaseSectionA);
Tp1(5,5)= exp(-1j*PhaseSectionB);
Tp1(6,6)= exp(-1j*PhaseSectionA);

Tp2 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp2(i,j) = 0;
        end
    end
end
Tp2(1,1)= exp(1j*PhaseSectionA);
Tp2(2,2)= exp(1j*PhaseSectionB);
Tp2(3,3)= exp(1j*PhaseSectionA);
Tp2(4,4)= exp(-1j*PhaseSectionA);
Tp2(5,5)= exp(-1j*PhaseSectionB);
Tp2(6,6)= exp(-1j*PhaseSectionA);
syms T1(CouplCoeff_k1)
T1 = [1, 0, 0, 0, 0, 0;
    0, 0, 1j/CouplCoeff_k1, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0;
    0, 1j/CouplCoeff_k1, 0, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
    0, 0, 0, 1, 0, 0;
    0, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, 0, -1j/CouplCoeff_k1;
    0, 0, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, -1j/CouplCoeff_k1, 0];
syms T2(CouplCoeff_k2)
T2 = [0, 1j/CouplCoeff_k2, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0, 0;
    1j/CouplCoeff_k2, 0, 0, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0;
    0, 0, 1, 0, 0, 0;
    1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0, 0, 0, -1j/CouplCoeff_k2, 0;
    0, 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0, -1j/CouplCoeff_k2, 0, 0;
    0, 0, 0, 0, 0, 1];

syms TransferMatrix(CouplCoeff_k1, CouplCoeff_k2, PhaseSectionA, PhaseSectionB)
T1Tp1 = T1*Tp1;
Tp2T1Tp1 = Tp2*T1*Tp1;
TransferMatrix = T2*Tp2*T1*Tp1;

syms x
Dispersiondiagram = charpoly(TransferMatrix,x);

x^6 
- (2*x^5*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2) 
- (x^3*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(CouplCoeff_k1*CouplCoeff_k2 + CouplCoeff_k1*CouplCoeff_k2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*CouplCoeff_k2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) - 2*CouplCoeff_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2) - 2*CouplCoeff_k2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3) 
+ (x^2*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(3*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2^3 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2^3 + exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2))/(CouplCoeff_k1^3*CouplCoeff_k2^3) 
+ (x^4*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(3*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2^3 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2^3 + exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2))/(CouplCoeff_k1^3*CouplCoeff_k2^3) 
- (x*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2) - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3) 
+ 1

%% Filippo asked me to change the Notation E1 = E3 from the previous notation.
%% New Notation ! %%%%%% Only two phases, a & b.
% Also, Calculate the 2x2 matrix E1(Z_3.2) = [a b ; c d]*E1(Z_0);
clear all; close all; clc

% Transfer Function Calculation
%%% We use the New Notation: 
%%% [E1+(0),E1-(0),E2+(0),E2-(0),E3+(0),E3-(0),E1+(L),E1-(L),E2+(L),E2-(L),E3+(L),E3-(L)]
%%% Careful

% Parameters
EffRefrIndex = 2.362; % Effective refractive index
Radius = 5e-6; % Radius of the "ring"
AngleSecB = 58.5*pi/180; % Angle of section B from the center of the "ring"
Freq = 2*pi*189e12;

LengthSectionB = 2*AngleSecB*Radius; % Length of section B
TotalLengthUnitCell = 2*pi*Radius + 2*LengthSectionB; % Total length of the unit cell
LightSpeed_c= 3e8; % Speed of light
Wavenumber_Vacuum = Freq/LightSpeed_c; % Wavenumber in vacuum

PhaseSectionA = Wavenumber_Vacuum*EffRefrIndex*pi*Radius/2; % Phase accumulated in the section A
PhaseSectionB = Wavenumber_Vacuum*EffRefrIndex*LengthSectionB; % Phase accumulated in the section B

% Matrix that accumulates travels from the beginning of the cell to the
% First Coupling Point (on the bottom half). Adapted to the new
% notation:

syms Tp1(PhaseSectionA, PhaseSectionB)
Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(1j*PhaseSectionA);
Tp1(2,2)= exp(-1j*PhaseSectionA);
Tp1(3,3)= exp(1j*PhaseSectionB);
Tp1(4,4)= exp(-1j*PhaseSectionB);
Tp1(5,5)= exp(1j*PhaseSectionA);
Tp1(6,6)= exp(-1j*PhaseSectionA);
    
syms Tp2(PhaseSectionA, PhaseSectionB)
Tp2 = Tp1;
    
syms T2(CouplCoeff_k2)
% T2 = zeros(6,6);
% T2(1,1) = 1;
% T2(2,2) = 1;
% T2(3,4) = -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;
% T2(3,5) = 1j/CouplCoeff_k2;
% T2(4,3) = 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;
% T2(4,6) = -1j/CouplCoeff_k2;
% T2(5,3) = 1j/CouplCoeff_k2;
% T2(5,6) = -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;
% T2(6,4) = -1j/CouplCoeff_k2;
% T2(6,5) = 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;

T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 1j/CouplCoeff_k2, 0;
    0, 0, 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0, 0, -1j/CouplCoeff_k2;
    0, 0, 1j/CouplCoeff_k2, 0, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;
    0, 0, 0, -1j/CouplCoeff_k2, 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0];

syms T1(CouplCoeff_k1)
% T1 = zeros(6,6);
% T1(5,5) = 1;
% T1(6,6) = 1;
% T1(1,2) = -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
% T1(1,3) = 1j/CouplCoeff_k1;
% T1(2,1) = 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
% T1(2,4) = -1j/CouplCoeff_k1;
% T1(3,1) = 1j/CouplCoeff_k1;
% T1(3,4) = -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
% T1(4,2) = -1j/CouplCoeff_k1;
% T1(4,3) = 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;

T1 = [ 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

% Transfer Matrix of the Unit Cell (Two coupling points: Bottom half &
% Top Half)

% TransferMatrix of the Finite Structure with a finite NumberOfUnitCells
syms TransferMatrix(CouplCoeff_k1, CouplCoeff_k2, PhaseSectionA, PhaseSectionB)
TransferMatrix = T2*Tp2*T1*Tp1;

%{
TransferMatrix =
 
[                                                                                                        0,                                                            -((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1,                          (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                                       0,                                                                    0,                                                                      0]
[                                                           ((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1,                                                                                                          0,                                                                                       0,                       -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                    0,                                                                      0]
[                                                                                                        0, -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2), ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2),                                                                                       0,                             (exp(PhaseSectionA*2i)*1i)/CouplCoeff_k2,                                                                      0]
[ -(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2),                                                                                                          0,                                                                                       0, ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2),                                                                    0,                             -(exp(-PhaseSectionA*2i)*1i)/CouplCoeff_k2]
[                             -(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i))/(CouplCoeff_k1*CouplCoeff_k2),                                                                                                          0,                                                                                       0,                               (1 - CouplCoeff_k1^2)^(1/2)/(CouplCoeff_k1*CouplCoeff_k2),                                                                    0, -(exp(-PhaseSectionA*2i)*(1 - CouplCoeff_k2^2)^(1/2)*1i)/CouplCoeff_k2]
[                                                                                                        0,                             -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i))/(CouplCoeff_k1*CouplCoeff_k2),                               (1 - CouplCoeff_k1^2)^(1/2)/(CouplCoeff_k1*CouplCoeff_k2),                                                                                       0, (exp(PhaseSectionA*2i)*(1 - CouplCoeff_k2^2)^(1/2)*1i)/CouplCoeff_k2,                                                                      0]
%}

% Now we get the Dispersion Diagram
syms x
DispersionRelation = charpoly(TransferMatrix,x);

%{
DispersionRelation = x^6 - (2*x^5*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2) - (x^3*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(CouplCoeff_k1*CouplCoeff_k2 + CouplCoeff_k1*CouplCoeff_k2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*CouplCoeff_k2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) - 2*CouplCoeff_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2) - 2*CouplCoeff_k2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3) + (x^2*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(3*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2^3 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2^3 + exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2))/(CouplCoeff_k1^3*CouplCoeff_k2^3) + (x^4*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(3*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2^3 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1^3*CouplCoeff_k2 - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2^3 + exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*CouplCoeff_k1*CouplCoeff_k2))/(CouplCoeff_k1^3*CouplCoeff_k2^3) - (x*exp(-PhaseSectionB*2i)*exp(-PhaseSectionA*4i)*(2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2) - 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3) + 1

DispersionRelationWrittenInMyNotation = 
x^6 
- (2*x^5*t1*t2)/(k1*k2) 
- (x^3*exp(-phib*2i)*exp(-phia*4i)*
(k1*k2 + k1*k2*exp(phib*4i)*exp(phia*8i) + 2*k1^2*exp(phib*2i)*exp(phia*4i)*t1*t2 + 2*k2^2*exp(phib*2i)*exp(phia*4i)*t1*t2 - 2*k1^2*exp(phib*2i)*exp(phia*4i)*t1*t2^3 - 2*k2^2*exp(phib*2i)*exp(phia*4i)*t1^3*t2))/(k1^3*k2^3) 
+ (x^2*exp(-phib*2i)*exp(-phia*4i)*(3*exp(phib*2i)*exp(phia*4i)*k1^3*k2^3 - 2*exp(phib*2i)*exp(phia*4i)*k1^3*k2 - 2*exp(phib*2i)*exp(phia*4i)*k1*k2^3 + exp(phib*2i)*exp(phia*4i)*k1*k2))/(k1^3*k2^3) 
+ (x^4*exp(-phib*2i)*exp(-phia*4i)*(3*exp(phib*2i)*exp(phia*4i)*k1^3*k2^3 - 2*exp(phib*2i)*exp(phia*4i)*k1^3*k2 - 2*exp(phib*2i)*exp(phia*4i)*k1*k2^3 + exp(phib*2i)*exp(phia*4i)*k1*k2))/(k1^3*k2^3) 
- (x*exp(-phib*2i)*exp(-phia*4i)*(2*exp(phib*2i)*exp(phia*4i)*t1*t2 - 2*exp(phib*2i)*exp(phia*4i)*t1*t2^3 - 2*exp(phib*2i)*exp(phia*4i)*t1^3*t2 + 2*exp(phib*2i)*exp(phia*4i)*t1^3*t2^3))/(k1^3*k2^3) 
+ 1
 
%}

% Now we get the Taux
syms Taux(CouplCoeff_k1, CouplCoeff_k2, PhaseSectionA, PhaseSectionB)
Taux = Tp2*T1*Tp1;
%{

Taux =
 
[                                                              0,                   -((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1, (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                 0,                     0,                      0]
[                 ((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1,                                                                 0,                                                              0, -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                     0,                      0]
[ (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                 0,                                                              0,                   -((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1,                     0,                      0]
[                                                              0, -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                 ((1 - CouplCoeff_k1^2)^(1/2)*1i)/CouplCoeff_k1,                                                                 0,                     0,                      0]
[                                                              0,                                                                 0,                                                              0,                                                                 0, exp(PhaseSectionA*2i),                      0]
[                                                              0,                                                                 0,                                                              0,                                                                 0,                     0, exp(-PhaseSectionA*2i)]
 
%}

%% New Notation! Four Coupling Coefficients!! %%%%%%%%%%%%%%%%%%%%%%%
% Also, Calculate the 2x2 matrix E1(Z_3.2) = [a b ; c d]*E1(Z_0);
clear all; close all; clc

% Transfer Function Calculation
%%% We use the New Notation: 
%%% [E1+(0),E1-(0),E2+(0),E2-(0),E3+(0),E3-(0),E1+(L),E1-(L),E2+(L),E2-(L),E3+(L),E3-(L)]
%%% Careful

% Parameters
EffRefrIndex = 2.362; % Effective refractive index
Radius = 5e-6; % Radius of the "ring"
AngleSecB = 58.5*pi/180; % Angle of section B from the center of the "ring"
Freq = 2*pi*189e12;

LengthSectionB = 2*AngleSecB*Radius; % Length of section B
TotalLengthUnitCell = 2*pi*Radius + 2*LengthSectionB; % Total length of the unit cell
LightSpeed_c= 3e8; % Speed of light
Wavenumber_Vacuum = Freq/LightSpeed_c; % Wavenumber in vacuum

PhaseSectionA = Wavenumber_Vacuum*EffRefrIndex*pi*Radius/2; % Phase accumulated in the section A
PhaseSectionB = Wavenumber_Vacuum*EffRefrIndex*LengthSectionB; % Phase accumulated in the section B

% Matrix that accumulates travels from the beginning of the cell to the
% First Coupling Point (on the bottom half). Adapted to the new
% notation:

syms Tp1(PhaseSectionA, PhaseSectionB)
Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(1j*PhaseSectionA);
Tp1(2,2)= exp(-1j*PhaseSectionA);
Tp1(3,3)= exp(1j*PhaseSectionB);
Tp1(4,4)= exp(-1j*PhaseSectionB);
Tp1(5,5)= exp(1j*PhaseSectionA);
Tp1(6,6)= exp(-1j*PhaseSectionA);
    
syms Tp2(PhaseSectionA, PhaseSectionB)
Tp2 = Tp1;
    
% The First pair of Coupling Coefficients is expressed in the following
% Transfer Matrices:
syms T2(CouplCoeff_k1)

T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 1j/CouplCoeff_k2, 0;
    0, 0, 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0, 0, -1j/CouplCoeff_k2;
    0, 0, 1j/CouplCoeff_k2, 0, 0, -1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2;
    0, 0, 0, -1j/CouplCoeff_k2, 1j*sqrt(1-CouplCoeff_k2^2)/CouplCoeff_k2, 0];

syms T1(CouplCoeff_k1)

T1 = [ 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

% The Second pair of Coupling Coefficients has the same matrices with a new
% Coefficient:
syms T2_2(CouplCoeff_k2_2)

T2_2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*sqrt(1-CouplCoeff_k2_2^2)/CouplCoeff_k2_2, 1j/CouplCoeff_k2_2, 0;
    0, 0, 1j*sqrt(1-CouplCoeff_k2_2^2)/CouplCoeff_k2_2, 0, 0, -1j/CouplCoeff_k2_2;
    0, 0, 1j/CouplCoeff_k2_2, 0, 0, -1j*sqrt(1-CouplCoeff_k2_2^2)/CouplCoeff_k2_2;
    0, 0, 0, -1j/CouplCoeff_k2_2, 1j*sqrt(1-CouplCoeff_k2_2^2)/CouplCoeff_k2_2, 0];

syms T1_2(CouplCoeff_k1_2)

T1_2 = [ 0, -1j*sqrt(1-CouplCoeff_k1_2^2)/CouplCoeff_k1_2, 1j/CouplCoeff_k1_2, 0, 0, 0;
    1j*sqrt(1-CouplCoeff_k1_2^2)/CouplCoeff_k1_2, 0, 0, -1j/CouplCoeff_k1_2, 0, 0;
    1j/CouplCoeff_k1_2, 0, 0, -1j*sqrt(1-CouplCoeff_k1_2^2)/CouplCoeff_k1_2, 0, 0;
    0, -1j/CouplCoeff_k1_2, 1j*sqrt(1-CouplCoeff_k1_2^2)/CouplCoeff_k1_2, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

% Transfer Matrix of the Unit Cell (Two coupling points: Bottom half &
% Top Half)

% TransferMatrix of the Finite Structure with a finite NumberOfUnitCells
syms TransferMatrix(CouplCoeff_k1, CouplCoeff_k2, PhaseSectionA, PhaseSectionB)
TransferMatrix_4CouplCoeff = T2_2*Tp2*T1_2*Tp1*T2*Tp2*T1*Tp1;

%{
TransferMatrix_4CouplCoeff = ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k1_2),                                                                                                                                                                                                        -((1 - CouplCoeff_k2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2),                                                                                                                                 (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2),                                                                                                                                                                         -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k1_2),                                   -(exp(PhaseSectionA*1i)*exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i))/(CouplCoeff_k2*CouplCoeff_k1_2),                                                                                                                                        0;

((1 - CouplCoeff_k2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2),                                                                                                                                                                                            ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k1_2),                                                                                                                                                                           -(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k1_2),                                                                                                                              -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2),                                                                                                                                      0,                                  -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i))/(CouplCoeff_k2*CouplCoeff_k1_2):

-exp(PhaseSectionA*1i)*((exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2) + (exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2)),                                                                                 -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                                                                                    ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2), exp(-PhaseSectionB*1i)*((exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2) + (exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2)), (exp(PhaseSectionA*2i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                                                              (1 - CouplCoeff_k2^2)^(1/2)/(CouplCoeff_k2*CouplCoeff_k2_2);

-(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2), exp(-PhaseSectionA*1i)*((exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2) + (exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2)), -exp(PhaseSectionB*1i)*((exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2) + (exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2)),                                                                                                    ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                                                            (1 - CouplCoeff_k2^2)^(1/2)/(CouplCoeff_k2*CouplCoeff_k2_2), -(exp(-PhaseSectionA*2i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2);

-(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2), exp(-PhaseSectionA*1i)*((exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2) + (exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2)), -exp(PhaseSectionB*1i)*((exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2) + (exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2)),                                                                                                                                  ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                            ((1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k2*CouplCoeff_k2_2),                               -(exp(-PhaseSectionA*2i)*(1 - CouplCoeff_k1_2^2)^(1/2)*1i)/(CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2);
 -exp(PhaseSectionA*1i)*((exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*(1
- CouplCoeff_k1^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2) + (exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2)),                                                                                                               -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                                                                                                                  ((1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2))/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2), exp(-PhaseSectionB*1i)*((exp(-PhaseSectionA*2i)*exp(-PhaseSectionB*1i)*1i)/(CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2) + (exp(PhaseSectionA*2i)*exp(PhaseSectionB*1i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)*1i)/(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k2_2)),                               (exp(PhaseSectionA*2i)*(1 - CouplCoeff_k1_2^2)^(1/2)*1i)/(CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2),                                              ((1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2))/(CouplCoeff_k2*CouplCoeff_k2_2)];

%}

% Now we get the Dispersion Diagram
syms Eigen
DispersionDiagram_4CouplCoeff = charpoly(TransferMatrix_4CouplCoeff,Eigen);

%{
The Dispersion Relation has been found to be:
DispersionDiagram_4CouplCoeff =
 
Eigen^6 - (Eigen*exp(-PhaseSectionB*4i)*exp(-PhaseSectionA*8i)*(2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3) - (Eigen^5*exp(-PhaseSectionB*4i)*exp(-PhaseSectionA*8i)*(2*CouplCoeff_k1^2*CouplCoeff_k2^3*CouplCoeff_k1_2^2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k2^2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3) - (Eigen^4*exp(-PhaseSectionB*4i)*exp(-PhaseSectionA*8i)*(- 3*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2^3*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*CouplCoeff_k1^2*CouplCoeff_k2^3*CouplCoeff_k1_2^2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 4*CouplCoeff_k1^3*CouplCoeff_k2^2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3) - (Eigen^3*exp(-PhaseSectionB*4i)*exp(-PhaseSectionA*8i)*(CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2 + 8*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*8i)*exp(PhaseSectionA*16i) - 4*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1^2*CouplCoeff_k2^3*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k1^3*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2^3*CouplCoeff_k1_2^2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k2^2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k1_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k1^2*CouplCoeff_k1_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1^2*CouplCoeff_k1_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k1_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3) - (Eigen^2*exp(-PhaseSectionB*4i)*exp(-PhaseSectionA*8i)*(- 3*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 4*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 4*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 4*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 4*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 4*CouplCoeff_k1^3*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 4*CouplCoeff_k2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 4*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(3/2) - 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(1/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(3/2) - 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + 2*CouplCoeff_k1*CouplCoeff_k1_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k1^3*CouplCoeff_k1_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k2^2)^(3/2)*(1 - CouplCoeff_k2_2^2)^(3/2) + 2*CouplCoeff_k2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i)*(1 - CouplCoeff_k1^2)^(3/2)*(1 - CouplCoeff_k1_2^2)^(3/2) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^3*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^3*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1*CouplCoeff_k2^3*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^3*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*4i)*exp(PhaseSectionA*8i) + CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2*CouplCoeff_k1_2^2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2*CouplCoeff_k1_2*CouplCoeff_k2_2^2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k2^2)^(1/2)*(1 - CouplCoeff_k1_2^2)^(1/2) + CouplCoeff_k1*CouplCoeff_k2^2*CouplCoeff_k1_2^2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2) + CouplCoeff_k1^2*CouplCoeff_k2^2*CouplCoeff_k1_2*CouplCoeff_k2_2*exp(PhaseSectionB*6i)*exp(PhaseSectionA*12i)*(1 - CouplCoeff_k1_2^2)^(1/2)*(1 - CouplCoeff_k2_2^2)^(1/2)))/(CouplCoeff_k1^3*CouplCoeff_k2^3*CouplCoeff_k1_2^3*CouplCoeff_k2_2^3) + 1
 
%% In a shorter notation:

DispersionDiagram_4CouplCoeff =
 
Eigen^6 + (exp(-phiB*4i)*exp(-phiA*8i)*(exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t1_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^4*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^4*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^2*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^2*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t1_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t1_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^4*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^4*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^2*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^2*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t1_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^4*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^4*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^2*t2_2^4 + 16*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^2*t2_2^2 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t1_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^4*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^4*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^2*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^2*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t1_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1^2 + exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^4*t1_2^2 + exp(phiB*4i)*exp(phiA*8i)*t2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^4*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^4*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^2*t2_2^4 - 8*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^2*t2_2^2 + 4*exp(phiB*4i)*exp(phiA*8i)*t2^2*t1_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t2^2 + exp(phiB*4i)*exp(phiA*8i)*t1_2^4*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1_2^4*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)*t1_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t1_2^2*t2_2^4 + 4*exp(phiB*4i)*exp(phiA*8i)*t1_2^2*t2_2^2 - 2*exp(phiB*4i)*exp(phiA*8i)*t1_2^2 + exp(phiB*4i)*exp(phiA*8i)*t2_2^4 - 2*exp(phiB*4i)*exp(phiA*8i)*t2_2^2 + exp(phiB*4i)*exp(phiA*8i)))/(k1^4*k2^4*k1_2^4*k2_2^4) - (Eigen^3*exp(-phiB*4i)*exp(-phiA*8i)*(k1^2*k2^2*k1_2^2*k2_2^2 + k1^2*k2^2*k1_2^2*k2_2^2*exp(phiB*8i)*exp(phiA*16i) + k1^2*k2^2*k1_2*k2_2^3*t1*t2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^2*k1_2^3*k2_2*t1*t2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^2*k1_2*k2_2^3*t1*t2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2^2*k1_2^3*k2_2*t1*t2*exp(phiB*6i)*exp(phiA*12i) + k1*k2^2*k1_2^2*k2_2^3*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) + k1^3*k2^2*k1_2^2*k2_2*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) + 2*k1*k2^2*k1_2^3*k2_2^2*t1*t1_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^3*k2^2*k1_2*k2_2^2*t1*t1_2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2*k1_2^3*k2_2^2*t1*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^3*k1_2*k2_2^2*t1*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1*k2^2*k1_2^2*k2_2^3*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) + k1^3*k2^2*k1_2^2*k2_2*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) + 2*k1^2*k2*k1_2^2*k2_2^3*t2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^3*k1_2^2*k2_2*t2*t2_2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2*k1_2^3*k2_2^2*t1*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2^3*k1_2*k2_2^2*t1*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1*k2^3*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1^3*k2*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1*k2^3*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1^3*k2*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) + 2*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^2*k1_2^2*k2_2^2*t2^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^2*k1_2^2*k2_2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2^3*k2_2^2*t1^3*t1_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^3*k2^2*k1_2*k2_2^2*t1*t1_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2^3*t2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2^3*k1_2^2*k2_2*t2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - k1*k2^2*k1_2^2*k2_2^3*t1^2*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2^2*k1_2*k2_2^3*t1*t2*t1_2^2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2*k1_2^3*k2_2^2*t1*t2^2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2^2*k1_2^3*k2_2*t1*t2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) - k1*k2^2*k1_2^2*k2_2^3*t1^2*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2^2*k1_2*k2_2^3*t1*t2*t1_2^2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2*k1_2^3*k2_2^2*t1*t2^2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2^2*k1_2^3*k2_2*t1*t2*t2_2^2*exp(phiB*6i)*exp(phiA*12i) - k1*k2^3*k1_2^2*k2_2^2*t1^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2^3*k1_2*k2_2^2*t1*t1_2^2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^3*k2*k1_2^2*k2_2^2*t2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^3*k2^2*k1_2^2*k2_2*t2*t1_2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) - k1*k2^3*k1_2^2*k2_2^2*t1^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2^3*k1_2*k2_2^2*t1*t1_2^2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^3*k2*k1_2^2*k2_2^2*t2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^3*k2^2*k1_2^2*k2_2*t2*t1_2*t2_2^2*exp(phiB*6i)*exp(phiA*12i) - 4*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2^2*k1_2^2*k2_2^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^3*k1_2*k2_2^3*t1*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^3*k2*k1_2^3*k2_2*t1*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 8*k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^3*k1_2*k2_2^3*t1*t2*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^3*k1_2*k2_2^3*t1^3*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^3*k2*k1_2^3*k2_2*t1*t2*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1^3*k2*k1_2^3*k2_2*t1*t2^3*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^3*k1_2^2*k2_2^3*t1^2*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^3*k2^2*k1_2^3*k2_2^2*t1*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^3*k1_2*k2_2^3*t1^3*t2*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^3*k2*k1_2^3*k2_2*t1*t2^3*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i)))/(k1^4*k2^4*k1_2^4*k2_2^4) - (Eigen*exp(-phiB*4i)*exp(-phiA*8i)*(2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t1_2^4*t2_2^3 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t1_2^4*t2_2 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t1_2^2*t2_2^3 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t1_2^2*t2_2 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t2_2^3 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2^3*t2_2 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t1_2^4*t2_2^3 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t1_2^4*t2_2 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t1_2^2*t2_2^3 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t1_2^2*t2_2 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t2_2^3 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^4*t2*t2_2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2^3*t2_2^4 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2^3*t2_2^2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2^3 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2*t2_2^4 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2*t2_2^2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^4*t1_2 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^3*t1_2^3*t2_2^3 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^3*t1_2^3*t2_2 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^3*t1_2*t2_2^3 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^3*t1_2*t2_2 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2^3*t2_2^4 + 8*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2^3*t2_2^2 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2^3 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2*t2_2^4 - 8*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2*t2_2^2 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2^2*t1_2 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2*t1_2^3*t2_2^3 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2*t1_2^3*t2_2 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2*t1_2*t2_2^3 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t2*t1_2*t2_2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2^3*t2_2^4 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2^3*t2_2^2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2^3 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2*t2_2^4 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2*t2_2^2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1^3*t1_2 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t1_2^4*t2_2^3 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t1_2^4*t2_2 + 8*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t1_2^2*t2_2^3 - 8*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t1_2^2*t2_2 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t2_2^3 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2^3*t2_2 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t1_2^4*t2_2^3 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t1_2^4*t2_2 - 8*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t1_2^2*t2_2^3 + 8*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t1_2^2*t2_2 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t2_2^3 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1^2*t2*t2_2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2^3*t2_2^4 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2^3*t2_2^2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2^3 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2*t2_2^4 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2*t2_2^2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^4*t1_2 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^3*t1_2^3*t2_2^3 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^3*t1_2^3*t2_2 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^3*t1_2*t2_2^3 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^3*t1_2*t2_2 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2^3*t2_2^4 - 8*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2^3*t2_2^2 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2^3 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2*t2_2^4 + 8*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2*t2_2^2 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2^2*t1_2 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2*t1_2^3*t2_2^3 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2*t1_2^3*t2_2 - 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2*t1_2*t2_2^3 + 2*k1*k2*k1_2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t1*t2*t1_2*t2_2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2^3*t2_2^4 + 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2^3*t2_2^2 - 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2^3 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2*t2_2^4 - 4*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2*t2_2^2 + 2*k1*k1_2*exp(phiB*4i)*exp(phiA*8i)*t1*t1_2 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t1_2^4*t2_2^3 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t1_2^4*t2_2 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t1_2^2*t2_2^3 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t1_2^2*t2_2 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t2_2^3 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2^3*t2_2 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t1_2^4*t2_2^3 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t1_2^4*t2_2 + 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t1_2^2*t2_2^3 - 4*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t1_2^2*t2_2 - 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t2_2^3 + 2*k2*k2_2*exp(phiB*4i)*exp(phiA*8i)*t2*t2_2))/(k1^4*k2^4*k1_2^4*k2_2^4) - (Eigen^2*exp(-phiB*4i)*exp(-phiA*8i)*(- k1^2*k1_2^2*t1^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - k2^2*k2_2^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k1_2^2*t2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k2_2^2*t1^2*exp(phiB*4i)*exp(phiA*8i) + k2^2*k1_2^2*k2_2^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k1_2^2*k2_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2^2*k2_2^2*t1^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2^2*k1_2^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k2_2^2*t1^2*t1_2^4*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k1_2^2*t2^2*t2_2^4*exp(phiB*4i)*exp(phiA*8i) - 2*k2^2*k1_2^2*k2_2^2*t1^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) + k2^2*k1_2^2*k2_2^2*t1^4*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k1_2^2*k2_2^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k1_2^2*k2_2^2*t2^4*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k1_2^2*t1^2*t2^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k1_2^2*t1^2*t2^4*t1_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k2^2*k2_2^2*t1^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k2^2*k2_2^2*t1^4*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k1_2^2*t1^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k1_2^2*t1^2*t1_2^2*t2_2^4*exp(phiB*4i)*exp(phiA*8i) + 2*k2^2*k2_2^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k2^2*k2_2^2*t2^2*t1_2^4*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k1_2*k2_2*t1*t2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^2*k1_2*k2_2*t1*t2*exp(phiB*6i)*exp(phiA*12i) + k1*k2^2*k1_2^2*k2_2*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2*k1_2*k2_2^2*t1*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1*k2^2*k1_2^2*k2_2*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2*k1_2*k2_2^2*t1*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1*k2*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1*k2*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) - 4*k1^2*k1_2^2*t1^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k1_2^2*t1^2*t2^2*t1_2^2*t2_2^4*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k1_2^2*t1^2*t2^4*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k1_2^2*t1^2*t2^4*t1_2^2*t2_2^4*exp(phiB*4i)*exp(phiA*8i) - 4*k2^2*k2_2^2*t1^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k2^2*k2_2^2*t1^2*t2^2*t1_2^4*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k2^2*k2_2^2*t1^4*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k2^2*k2_2^2*t1^4*t2^2*t1_2^4*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1*k2^2*k1_2*k2_2^2*t1^3*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t2*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 2*k1^2*k2*k1_2^2*k2_2*t2^3*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2^3*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 2*k1*k2^2*k1_2*k2_2^2*t1^3*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2*k1_2^2*k2_2*t2^3*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - k1*k2^2*k1_2^2*k2_2*t1^2*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2^2*k1_2*k2_2*t1*t2*t1_2^2*exp(phiB*2i)*exp(phiA*4i) + 2*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k2*k1_2*k2_2^2*t1*t2^2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2^2*k1_2*k2_2*t1*t2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) - k1*k2^2*k1_2^2*k2_2*t1^2*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2^2*k1_2*k2_2*t1*t2*t1_2^2*exp(phiB*6i)*exp(phiA*12i) + 2*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t2_2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k2*k1_2*k2_2^2*t1*t2^2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2^2*k1_2*k2_2*t1*t2*t2_2^2*exp(phiB*6i)*exp(phiA*12i) - k1*k2*k1_2^2*k2_2^2*t1^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1^2*k2*k1_2*k2_2^2*t1*t1_2^2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1*k2*k1_2^2*k2_2^2*t2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) - k1*k2^2*k1_2^2*k2_2*t2*t1_2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) + 2*k1*k2^2*k1_2*k2_2^2*t1*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2*k1_2^2*k2_2*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) - k1*k2*k1_2^2*k2_2^2*t1^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1^2*k2*k1_2*k2_2^2*t1*t1_2^2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1*k2*k1_2^2*k2_2^2*t2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) - k1*k2^2*k1_2^2*k2_2*t2*t1_2*t2_2^2*exp(phiB*6i)*exp(phiA*12i) - 4*k1*k2*k1_2*k2_2*t1*t2*t1_2^3*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1*t2^3*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1*t2^3*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1^3*t2*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1^3*t2*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1^3*t2^3*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 4*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1*t2^3*t1_2^3*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1^3*t2*t1_2^3*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1^3*t2^3*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1^3*t2^3*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2^2*k1_2*k2_2^2*t1^3*t2^2*t1_2^3*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2*k1_2^2*k2_2*t1^2*t2^3*t1_2^2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1*t2*t1_2*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1*t2*t1_2^3*t2_2*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1*t2^3*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 4*k1*k2*k1_2*k2_2*t1^3*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i) + k1*k2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1*k2^2*k1_2^2*k2_2*t1^2*t2*t1_2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2*k1_2*k2_2^2*t1*t2^2*t1_2^2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^2*k1_2*k2_2*t1*t2*t1_2^2*t2_2^2*exp(phiB*2i)*exp(phiA*4i) - 4*k1*k2^2*k1_2*k2_2^2*t1*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2*k1_2^2*k2_2*t1^2*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) - 4*k1*k2*k1_2*k2_2*t1^3*t2^3*t1_2^3*t2_2^3*exp(phiB*4i)*exp(phiA*8i) + k1*k2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1*k2^2*k1_2^2*k2_2*t1^2*t2*t1_2*t2_2^2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2*k1_2*k2_2^2*t1*t2^2*t1_2^2*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2^2*k1_2*k2_2*t1*t2*t1_2^2*t2_2^2*exp(phiB*6i)*exp(phiA*12i)))/(k1^4*k2^4*k1_2^4*k2_2^4) - (Eigen^5*exp(-phiB*4i)*exp(-phiA*8i)*(2*t2*t2_2*exp(phiB*4i)*exp(phiA*8i)*k1^4*k2^3*k1_2^4*k2_2^3 + 2*t1*t1_2*exp(phiB*4i)*exp(phiA*8i)*k1^3*k2^4*k1_2^3*k2_2^4 + 2*t1*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i)*k1^3*k2^3*k1_2^3*k2_2^3))/(k1^4*k2^4*k1_2^4*k2_2^4) - (Eigen^4*exp(-phiB*4i)*exp(-phiA*8i)*(k1^2*k2^2*k1_2^4*k2_2^2*t1^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k1_2^2*k2_2^4*t2^2*exp(phiB*4i)*exp(phiA*8i) + k1^4*k2^2*k1_2^2*k2_2^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^4*k1_2^2*k2_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k2^4*k1_2^2*k2_2^4*t1^2*t1_2^2*exp(phiB*4i)*exp(phiA*8i) - k1^4*k2^2*k1_2^4*k2_2^2*t2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + k1^2*k2^2*k1_2^3*k2_2^3*t1*t2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^2*k1_2^3*k2_2^3*t1*t2*exp(phiB*6i)*exp(phiA*12i) + k1^3*k2^2*k1_2^2*k2_2^3*t2*t1_2*exp(phiB*2i)*exp(phiA*4i) + k1^2*k2^3*k1_2^3*k2_2^2*t1*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1^3*k2^2*k1_2^2*k2_2^3*t2*t1_2*exp(phiB*6i)*exp(phiA*12i) + k1^2*k2^3*k1_2^3*k2_2^2*t1*t2_2*exp(phiB*6i)*exp(phiA*12i) + k1^3*k2^3*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*2i)*exp(phiA*4i) + k1^3*k2^3*k1_2^2*k2_2^2*t1_2*t2_2*exp(phiB*6i)*exp(phiA*12i) + 2*k1^3*k2^2*k1_2^3*k2_2^2*t1*t2^2*t1_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^3*k1_2^2*k2_2^3*t1^2*t2*t2_2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^3*k2^2*k1_2^3*k2_2^2*t1*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) + 2*k1^2*k2^3*k1_2^2*k2_2^3*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) - k1^2*k2^2*k1_2^2*k2_2^2*t1^2*t2^2*t1_2^2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^2*k2^3*k1_2^2*k2_2^3*t1^2*t2*t1_2^2*t2_2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^3*k2^2*k1_2^3*k2_2^2*t1*t2^2*t1_2*t2_2^2*exp(phiB*4i)*exp(phiA*8i) - 4*k1^3*k2^3*k1_2^3*k2_2^3*t1*t2*t1_2*t2_2*exp(phiB*4i)*exp(phiA*8i)))/(k1^4*k2^4*k1_2^4*k2_2^4)
 
%}
%% New notation!!! Three Phase Delays and 1 Coupling Coefficient.
% % 
clear all; close  all; clc;
syms Tp1(PhaseSectionA, PhaseSectionB)
Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(1j*PhaseSectionA);
Tp1(2,2)= exp(-1j*PhaseSectionA);
Tp1(3,3)= exp(1j*PhaseSectionB);
Tp1(4,4)= exp(-1j*PhaseSectionB);
Tp1(5,5)= exp(1j*PhaseSectionA);
Tp1(6,6)= exp(-1j*PhaseSectionA);

syms Tp2(PhaseSectionA, PhaseSectionD)
Tp2 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp2(i,j) = 0;
        end
    end
end
Tp2(1,1)= exp(1j*PhaseSectionA);
Tp2(2,2)= exp(-1j*PhaseSectionA);
Tp2(3,3)= exp(1j*PhaseSectionD);
Tp2(4,4)= exp(-1j*PhaseSectionD);
Tp2(5,5)= exp(1j*PhaseSectionA);
Tp2(6,6)= exp(-1j*PhaseSectionA);

syms T2(CouplCoeff_k1, Trans_k1)

T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0;
    0, 0, 1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1;
    0, 0, 1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1;
    0, 0, 0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0];

syms T1(CouplCoeff_k1, Trans_k1)

T1 = [ 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

syms TransferMatrix(CouplCoeff_k1, Trans_k1, PhaseSectionA, PhaseSectionB, PhaseSectionD)
TransferMatrix = T2*Tp2*T1*Tp1;
%{
TransferMatrix
 
TransferMatrix =
 
[                                                                       0,                                              -(Trans_k1*1i)/CouplCoeff_k1,            (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                         0,                                                 0,                                                   0]
[                                             (Trans_k1*1i)/CouplCoeff_k1,                                                                         0,                                                                         0,         -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                 0,                                                   0]
[                                                                       0, -(Trans_k1*exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2, (Trans_k1^2*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0,          (exp(PhaseSectionA*2i)*1i)/CouplCoeff_k1,                                                   0]
[ -(Trans_k1*exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0,                                                                         0, (Trans_k1^2*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                 0,          -(exp(-PhaseSectionA*2i)*1i)/CouplCoeff_k1]
[          -(exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0,                                                                         0,   (Trans_k1*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                 0, -(Trans_k1*exp(-PhaseSectionA*2i)*1i)/CouplCoeff_k1]
[                                                                       0,          -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,   (Trans_k1*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0, (Trans_k1*exp(PhaseSectionA*2i)*1i)/CouplCoeff_k1,                                                   0]
 
%}
syms Eigen
DispersionDiagram = charpoly(TransferMatrix,Eigen);

%{
DispersionDiagram =
 
Eigen^6 + 
(exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^8 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^6 + 6*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^4 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^2 + exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)))/CouplCoeff_k1^8 
- (Eigen^5*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^6*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^6*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8 
- (Eigen^2*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(- exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^4*Trans_k1^4 + 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^6 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^4 + 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^2))/CouplCoeff_k1^8 
- (Eigen*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^2*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^2*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) - 2*CouplCoeff_k1^2*Trans_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) - 2*CouplCoeff_k1^2*Trans_k1^4*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) + CouplCoeff_k1^2*Trans_k1^6*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^2*Trans_k1^6*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8 
+ (Eigen^4*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(- 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^6*Trans_k1^2 + exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^4*Trans_k1^4))/CouplCoeff_k1^8 
- (Eigen^3*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^4 + CouplCoeff_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionD*2i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^4*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + 2*CouplCoeff_k1^4*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) - 2*CouplCoeff_k1^4*Trans_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) - 2*CouplCoeff_k1^4*Trans_k1^4*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8

%}
syms TransferAux(CouplCoeff_k1, Trans_k1, PhaseSectionA, PhaseSectionB, PhaseSectionD)
TransferAux = Tp2*T1*Tp1;

%{
 
TransferAux =
 
[                                                              0,                                      -(Trans_k1*1i)/CouplCoeff_k1,           (exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                         0,                     0,                      0]
[                                    (Trans_k1*1i)/CouplCoeff_k1,                                                                 0,                                                                        0,         -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                     0,                      0]
[ (exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i)*1i)/CouplCoeff_k1,                                                                 0,                                                                        0, -(Trans_k1*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i)*1i)/CouplCoeff_k1,                     0,                      0]
[                                                              0, -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i)*1i)/CouplCoeff_k1, (Trans_k1*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i)*1i)/CouplCoeff_k1,                                                                         0,                     0,                      0]
[                                                              0,                                                                 0,                                                                        0,                                                                         0, exp(PhaseSectionA*2i),                      0]
[                                                              0,                                                                 0,                                                                        0,                                                                         0,                     0, exp(-PhaseSectionA*2i)]
 
%}

%% New Notation!! 3 phase delays & 1 couplcoeff. Changes in sign
%% New notation!!! Three Phase Delays and 1 Coupling Coefficient.
% % 
clear all; close  all; clc;
syms Tp1(PhaseSectionA, PhaseSectionB)
Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(-1j*PhaseSectionA);
Tp1(2,2)= exp(1j*PhaseSectionA);
Tp1(3,3)= exp(-1j*PhaseSectionB);
Tp1(4,4)= exp(1j*PhaseSectionB);
Tp1(5,5)= exp(-1j*PhaseSectionA);
Tp1(6,6)= exp(1j*PhaseSectionA);

syms Tp2(PhaseSectionA, PhaseSectionD)
Tp2 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp2(i,j) = 0;
        end
    end
end
Tp2(1,1)= exp(-1j*PhaseSectionA);
Tp2(2,2)= exp(1j*PhaseSectionA);
Tp2(3,3)= exp(-1j*PhaseSectionD);
Tp2(4,4)= exp(1j*PhaseSectionD);
Tp2(5,5)= exp(-1j*PhaseSectionA);
Tp2(6,6)= exp(1j*PhaseSectionA);

syms T2(CouplCoeff_k1, Trans_k1)

T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0;
    0, 0, 1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1;
    0, 0, 1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1;
    0, 0, 0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0];

syms T1(CouplCoeff_k1, Trans_k1)

T1 = [ 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

syms TransferMatrix(CouplCoeff_k1, Trans_k1, PhaseSectionA, PhaseSectionB, PhaseSectionD)
TransferMatrix = T2*Tp2*T1*Tp1;
%{
TransferMatrix
[                                                                         0,                                            -(Trans_k1*1i)/CouplCoeff_k1,          (exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                         0,                                                  0,                                                  0]
[                                               (Trans_k1*1i)/CouplCoeff_k1,                                                                       0,                                                                         0,           -(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                  0,                                                  0]
[                                                                         0, -(Trans_k1*exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2, (Trans_k1^2*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0,          (exp(-PhaseSectionA*2i)*1i)/CouplCoeff_k1,                                                  0]
[ -(Trans_k1*exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                       0,                                                                         0, (Trans_k1^2*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                  0,          -(exp(PhaseSectionA*2i)*1i)/CouplCoeff_k1]
[          -(exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                       0,                                                                         0,   (Trans_k1*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i))/CouplCoeff_k1^2,                                                  0, -(Trans_k1*exp(PhaseSectionA*2i)*1i)/CouplCoeff_k1]
[                                                                         0,          -(exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,   (Trans_k1*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i))/CouplCoeff_k1^2,                                                                         0, (Trans_k1*exp(-PhaseSectionA*2i)*1i)/CouplCoeff_k1,                                                  0]
 

%}
syms Eigen
DispersionDiagram = charpoly(TransferMatrix,Eigen);

%{
DispersionDiagram =
 
Eigen^6 
+ (exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^8 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^6 + 6*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^4 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*Trans_k1^2 + exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)))/CouplCoeff_k1^8 
- (Eigen^5*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^6*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^6*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8 
- (Eigen^2*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(- exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^4*Trans_k1^4 + 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^6 - 4*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^4 + 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^2*Trans_k1^2))/CouplCoeff_k1^8 
- (Eigen*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^2*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^2*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) - 2*CouplCoeff_k1^2*Trans_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) - 2*CouplCoeff_k1^2*Trans_k1^4*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) + CouplCoeff_k1^2*Trans_k1^6*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + CouplCoeff_k1^2*Trans_k1^6*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8 
+ (Eigen^4*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(- 2*exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^6*Trans_k1^2 + exp(PhaseSectionB*1i)*exp(PhaseSectionA*4i)*exp(PhaseSectionD*1i)*CouplCoeff_k1^4*Trans_k1^4))/CouplCoeff_k1^8 
- (Eigen^3*exp(-PhaseSectionB*1i)*exp(-PhaseSectionA*4i)*exp(-PhaseSectionD*1i)*(CouplCoeff_k1^4 + CouplCoeff_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionD*2i)*exp(PhaseSectionA*8i) + 2*CouplCoeff_k1^4*Trans_k1^2*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) + 2*CouplCoeff_k1^4*Trans_k1^2*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i) - 2*CouplCoeff_k1^4*Trans_k1^4*exp(PhaseSectionB*2i)*exp(PhaseSectionA*4i) - 2*CouplCoeff_k1^4*Trans_k1^4*exp(PhaseSectionA*4i)*exp(PhaseSectionD*2i)))/CouplCoeff_k1^8
 

%}
syms TransferAux(CouplCoeff_k1, Trans_k1, PhaseSectionA, PhaseSectionB, PhaseSectionD)
TransferAux = Tp2*T1*Tp1;

%{
 
TransferAux =
 
[                                                                0,                                    -(Trans_k1*1i)/CouplCoeff_k1,         (exp(-PhaseSectionA*1i)*exp(-PhaseSectionB*1i)*1i)/CouplCoeff_k1,                                                                         0,                      0,                     0]
[                                      (Trans_k1*1i)/CouplCoeff_k1,                                                               0,                                                                        0,           -(exp(PhaseSectionA*1i)*exp(PhaseSectionB*1i)*1i)/CouplCoeff_k1,                      0,                     0]
[ (exp(-PhaseSectionA*1i)*exp(-PhaseSectionD*1i)*1i)/CouplCoeff_k1,                                                               0,                                                                        0, -(Trans_k1*exp(PhaseSectionB*1i)*exp(-PhaseSectionD*1i)*1i)/CouplCoeff_k1,                      0,                     0]
[                                                                0, -(exp(PhaseSectionA*1i)*exp(PhaseSectionD*1i)*1i)/CouplCoeff_k1, (Trans_k1*exp(-PhaseSectionB*1i)*exp(PhaseSectionD*1i)*1i)/CouplCoeff_k1,                                                                         0,                      0,                     0]
[                                                                0,                                                               0,                                                                        0,                                                                         0, exp(-PhaseSectionA*2i),                     0]
[                                                                0,                                                               0,                                                                        0,                                                                         0,                      0, exp(PhaseSectionA*2i)]
 
%}



%%
clear all; close  all; clc;
syms Tp1(PhaseSectionA, PhaseSectionB)
Tp1 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp1(i,j) = 0;
        end
    end
end
Tp1(1,1)= exp(1j*PhaseSectionA);
Tp1(2,2)= exp(-1j*PhaseSectionA);
Tp1(3,3)= exp(1j*PhaseSectionB);
Tp1(4,4)= exp(-1j*PhaseSectionB);
Tp1(5,5)= exp(1j*PhaseSectionA);
Tp1(6,6)= exp(-1j*PhaseSectionA);

syms Tp2(PhaseSectionA, PhaseSectionD)
Tp2 = sym('a%d%d', [6 6]);
for i=1:6
    for j=1:6
        if i~=j
            Tp2(i,j) = 0;
        end
    end
end
Tp2(1,1)= exp(1j*PhaseSectionA);
Tp2(2,2)= exp(-1j*PhaseSectionA);
Tp2(3,3)= exp(1j*PhaseSectionD);
Tp2(4,4)= exp(-1j*PhaseSectionD);
Tp2(5,5)= exp(1j*PhaseSectionA);
Tp2(6,6)= exp(-1j*PhaseSectionA);

syms T2(CouplCoeff_k1, Trans_k1)

T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0;
    0, 0, 1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1;
    0, 0, 1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1;
    0, 0, 0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0];

syms T1(CouplCoeff_k1, Trans_k1)

T1 = [ 0, -1j*Trans_k1/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*Trans_k1/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*Trans_k1/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*Trans_k1/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];

syms TransferMatrix(CouplCoeff_k1, Trans_k1, PhaseSectionA, PhaseSectionB, PhaseSectionD)
TransferMatrix = T2*Tp2*T1*Tp1;

%% Let's calculate the scattering matrix:
T_11 = TransferMatrix(1:3,1:3);
T_12 = TransferMatrix(1:3,4:6);
T_21 = TransferMatrix(4:6,1:3);
T_22 = TransferMatrix(4:6,4:6);

S_11 = -T_22^-1*T_21;
S_12 = T_22^-1;
S_21 = T_11 - T_12*T_22^-1*T_21;
S_22 = T_12*T_22^-1;

ScatteringMatrix(1:3,1:3) = S_11;
ScatteringMatrix(1:3,4:6) = S_12;
ScatteringMatrix(4:6,1:3) = S_21;
ScatteringMatrix(4:6,4:6) = S_22;


%%
syms DispRel(a, b, c)
DispRel = (a-b)^3*(a-c)^3;

