function [TransferMatrix] = UnitCell_ASOW_TransferMatrix (Freq, CouplCoeff_k1, Radius, Alpha, Alpha_2)
    % Returns the Transfer Function of Scheuer's SOW when given the
    % parameters
    EffRefrIndex = 2.362;
    SpeedLight = 3e8; % m/s
    % Freq is in THz
    pA = (pi^2*EffRefrIndex*Radius*1e12)/SpeedLight;
    pB = (4*pi*EffRefrIndex*Radius*Alpha*1e12)/SpeedLight;
    pD = (4*pi*EffRefrIndex*Radius*Alpha_2*1e12)/SpeedLight;
    
    PhaseSectionA = Freq * pA;
    PhaseSectionB = Freq * pB;
    PhaseSectionD = Freq * pD;
    
   
    % Tp1: Transfer Matrix from Z0+ to Zc- (from begining of the unit cell to first coupling point)
    % Tp2: Transfer Matrix from Zc+ to Zc'- (from first to second coupling
    % point)

    Tp1(1,1)= exp(1j*PhaseSectionA);
    Tp1(2,2)= exp(-1j*PhaseSectionA);
    Tp1(3,3)= exp(1j*PhaseSectionB);
    Tp1(4,4)= exp(-1j*PhaseSectionB);
    Tp1(5,5)= exp(1j*PhaseSectionA);
    Tp1(6,6)= exp(-1j*PhaseSectionA);
    
    Tp2(1,1)= exp(1j*PhaseSectionA);
    Tp2(2,2)= exp(-1j*PhaseSectionA);
    Tp2(3,3)= exp(1j*PhaseSectionD);
    Tp2(4,4)= exp(-1j*PhaseSectionD);
    Tp2(5,5)= exp(1j*PhaseSectionA);
    Tp2(6,6)= exp(-1j*PhaseSectionA);
     
T1 = [ 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
    1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
    1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0;
    0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];
 
T2 = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0;
    0, 0, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1;
    0, 0, 1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
    0, 0, 0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0]; 
       
TransferMatrix = T2*Tp2*T1*Tp1;
%
 % The TransferMatrix below has been calculated using Symbolic Math.

end