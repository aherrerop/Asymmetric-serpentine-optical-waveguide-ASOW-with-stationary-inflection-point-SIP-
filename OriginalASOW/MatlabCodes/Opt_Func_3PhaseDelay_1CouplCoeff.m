function [sigma] = Opt_Func_3PhaseDelay_1CouplCoeff(X)
% Sigma is the eigenvector coalescence parameter (the more similar the
% eigenvalues, the smaller it is (tends to zero)).
% X is a vector of the Parameters to tune:
Freq = 193.54;
CouplCoeff_k1 = X(1);
pA = X(2);
pB = X(3);
pD = X(4);


% We get the Transfer Matrix with these parameters:
[TransferMatrix] = get_MyBrokenDownSOW_UnitCell_3_PhaseDelay_1CouplCoeff (Freq, CouplCoeff_k1, pA, pB, pD);
% sigma = SIP_Check_Vectors(TransferMatrix);
sigma = SIP_Check_Eigenvectors(TransferMatrix);
% sigma = Sixth_Check_Eigenvectors(TransferMatrix);

% Enforce Solutions to be strictly physical:

% For the Coupling Coefficients:
if CouplCoeff_k1 > 0.95 && CouplCoeff_k1 <0.05 % && CouplCoeff_k2 > 0.95 && CouplCoeff_k2 < 0.05 && CouplCoeff_k3 > 0.95 && CouplCoeff_k3 < 0.05 && CouplCoeff_k4 > 0.95 && CouplCoeff_k4 < 0.05
    sigma = 1e3;
end
X)
% For the Phases
if Freq*pA > pi && Freq*pA < -pi && Freq*pB > pi && Freq*pB <-pi
% if pA < 0.86589 && pB > 1.14063 && pB < 1.120129
    sigma = 1e3; 
end

if Freq*pD > pi && Freq*pD < -pi 
% if pD > 1.14063 && pD < 1.120129
    sigma = 1e3; 
end

end




