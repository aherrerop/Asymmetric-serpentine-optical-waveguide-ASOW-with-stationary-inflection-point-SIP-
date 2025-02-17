function [PortValues, TransferMatrix] = get_MyBrokenDownSOW_PortValues_ASOW_SIP_Parameters (NumberOfUnitCells, Freq, CouplCoeff_k1, Radius, Alpha, Alpha_2)
    % Transfer Function Calculation
    %%% We use the New Notation:
    %%% [E1+(0),E1-(0),E2+(0),E2-(0),E3+(0),E3-(0),E1+(L),E1-(L),E2+(L),E2-(L),E3+(L),E3-(L)]
    %%% Careful
    % 12 Unknowns. 12 Equations. One Solution.
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

    Tp1(1,1)= exp(-1j*PhaseSectionA);
    Tp1(2,2)= exp(1j*PhaseSectionA);
    Tp1(3,3)= exp(-1j*PhaseSectionB);
    Tp1(4,4)= exp(1j*PhaseSectionB);
    Tp1(5,5)= exp(-1j*PhaseSectionA);
    Tp1(6,6)= exp(1j*PhaseSectionA);
    
    Tp2(1,1)= exp(-1j*PhaseSectionA);
    Tp2(2,2)= exp(1j*PhaseSectionA);
    Tp2(3,3)= exp(-1j*PhaseSectionD);
    Tp2(4,4)= exp(1j*PhaseSectionD);
    Tp2(5,5)= exp(-1j*PhaseSectionA);
    Tp2(6,6)= exp(1j*PhaseSectionA);
     
    T1 = -1*[ 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0, 0, 0;
        1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1, 0, 0;
        1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0;
        0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, 0;
        0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 1];
 
    T2 = -1*[1, 0, 0, 0, 0, 0;
        0, 1, 0, 0, 0, 0;
        0, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 1j/CouplCoeff_k1, 0;
        0, 0, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0, 0, -1j/CouplCoeff_k1;
        0, 0, 1j/CouplCoeff_k1, 0, 0, -1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1;
        0, 0, 0, -1j/CouplCoeff_k1, 1j*sqrt(1-CouplCoeff_k1^2)/CouplCoeff_k1, 0];
       
    TransferMatrix_UnitCell = T2*Tp2*T1*Tp1;
    T_aux = Tp2 * T1 * Tp1;
    % TransferMatrix of the Finite Structure with a finite NumberOfUnitCells                                              
    TransferMatrix = T_aux * ((TransferMatrix_UnitCell)^(NumberOfUnitCells-1));

    
    %  SystemMatrix*PortValues = b = [1, 0, ..., 0]
    % The Ports are ordered like: 
    % [E1+(0),E1-(0),E2+(0),E2-(0),E3+(0),E3-(0),E1+(L),E1-(L),E2+(L),E2-(L),E3+(L),E3-(L)]
    SystemMatrix = zeros (12,12);
    % E1+(0) = 1;
    SystemMatrix(1,1) = 1;
    % 6 Equations that come from the TransferMatrix
    SystemMatrix(2:7,1:6) = TransferMatrix;
    SystemMatrix(2:7,7:12) = -eye(6);
    % E1-(L) = 0. The Waveguide is perfectly matched
    SystemMatrix(8,8) = 1;
    % Two Equations at Z = 0
    SystemMatrix(9,3) = 1;
    SystemMatrix(9,6) = -1;
    SystemMatrix(10,4) = -1;
    SystemMatrix(10,5) = 1;
    % Two Equations at Z = L
    SystemMatrix(11,10) = 1;
    SystemMatrix(11,11) = -1;
    SystemMatrix(12,9) = -1;
    SystemMatrix(12,12) = 1;
    
    % E1+(0) = 1.  
    b = zeros(12,1);
    b(1,1) = 1;
    PortValues = zeros(12,1);
    
    % We solve the system for the value of the field amplitudes at the ports.
    PortValues = SystemMatrix\b;
 

 
 
                                                             
 
end
