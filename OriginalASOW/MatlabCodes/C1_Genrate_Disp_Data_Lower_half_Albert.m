
clear all;

% Parameters 6DBE:
Freq_Center = 193.54;
Radius = 10e-6;
CouplCoeff_k1 = 0.408249494346565;
Alpha = 1.06438571584396;
Alpha_2 = 1.03641573751797;

% Frequency Sweep
Freqmax = Freq_Center + 0.040;
Freqmin = Freq_Center - 0.040;
Freqsteps = 5001;

Freq_low = linspace(Freqmin, Freq_Center, Freqsteps);
Freq_high = linspace(Freq_Center,Freqmax, Freqsteps);

Freq_low_plt = Freq_low./Freq_Center;
Freq_high_plt = Freq_high./Freq_Center;
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
for ii = 1:Freqsteps    
    TransferMatrices_Vec(:,:,ii) =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq_low(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    Sigma(ii) = SIP_Check_Vectors(TransferMatrices_Vec(:,:,ii));
end

for ii = 1:6
    kdlow1 = Eigenvalues(1,:);
    kdlow2 = Eigenvalues(2,:);
    kdlow3 = Eigenvalues(3,:);
    kdlow4 = Eigenvalues(4,:);
    kdlow5 = Eigenvalues(5,:);
    kdlow6 = Eigenvalues(6,:);
end