function [sigma] = Opt_Func_3PhaseDelay_1CouplCoeff_Parameters(X)
% Sigma is the eigenvector coalescence parameter (the more similar the
% eigenvalues, the smaller it is (tends to zero)).
% X is a vector of the Parameters to tune:
Freq = 193.54;
CouplCoeff_k1 = X(1);
Alpha = X(2);
Alpha_2 = X(3);
Radius = 6e-6;


% We get the Transfer Matrix with these parameters:
[TransferMatrix] = get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq, CouplCoeff_k1, Radius, Alpha, Alpha_2);
sigma = SIP_Check_Vectors(TransferMatrix);
% sigma = SIP_Check_Eigenvectors(TransferMatrix);
% sigma = Sixth_Check_Eigenvectors(TransferMatrix);

% Enforce Solutions to be strictly physical:

% For the Coupling Coefficients:
if CouplCoeff_k1 > 0.95 && CouplCoeff_k1 <0.05 % && CouplCoeff_k2 > 0.95 && CouplCoeff_k2 < 0.05 && CouplCoeff_k3 > 0.95 && CouplCoeff_k3 < 0.05 && CouplCoeff_k4 > 0.95 && CouplCoeff_k4 < 0.05
    sigma = 1e3;
end

% For the Angles:
if Alpha > 1.035 && Alpha < 1 && Alpha_2 > 1.035 && Alpha_2 < 1
    sigma = 1e3;
end

% Get SIP far from center:
[Eigenvectors, EigenValues_X] = eig(TransferMatrix);
Eigenvalues = diag(EigenValues_X);
[Eigenvalues_X,index_X] = sort(Eigenvalues);
Ksd_X = -log(Eigenvalues_X)./(1j);
Ksd_pi = Ksd_X/pi;
if abs(Ksd_pi) > 0.5
    sigma = 1e3;
end
end




