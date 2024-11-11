% DispersionDiagramDrawing
clear all; close all; clc;

% Parameters SIP:
Freq_Center = 193.54;
CouplCoeff = 0.49832327234602;
Radius = 10e-6;
Alpha = 1.15240339832846;
Alpha_2 = 0.980630321583591;


Loss = 0;%4.7e-6; % Loss is positive!        
Gain = 0; % Gain is negative!          
ImagRefrIndex = Gain + Loss; % Imaginary part of the Effective refractive index. Negative if gain.
% Constants
EffRefrIndex = 2.362  -1j*ImagRefrIndex; %  0.01 is a 5nm width change. No SIP distortion, only displacement.
SpeedLight = 3e8; % m/s


% Frequency Sweep
Freqmax = Freq_Center + 0.05;
Freqmin = Freq_Center - 0.05;
Freqsteps = 5000;

Freq = linspace(Freqmin, Freqmax, Freqsteps);
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
NumberOfUnitCells = 9;
for jj = 1:Freqsteps    
%     TransferMatrices_Vec(:,:,jj) =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq(jj), CouplCoeff, Radius, Alpha, Alpha_2);
    [~, TransferMatrices_Vec(:,:,jj)] = ASOW_Gain_PortValues (NumberOfUnitCells, Freq(jj), CouplCoeff, Radius, Alpha, Alpha_2);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
%     Sigma(ii) = SIP_Finding_Hyperdistance(TransferMatrices_Vec(:,:,ii));
end
% 
% Sigma_peaks = -findpeaks(-Sigma);
% 
% index_Top = find(Sigma == Sigma_peaks(3));
% index_Bottom = find(Sigma == Sigma_peaks(1));
% Freq_RBE_Top = Freq(index_Top);
% Freq_RBE_Bottom = Freq(index_Bottom);
% % Freq_RBE_Top = Freq_Center;
% % Freq_RBE_Bottom = Freq_Center;

% plot dispersion
    
Eigenvalues_plt = Eigenvalues'/pi;
Eigenvalues_plt_real = real(Eigenvalues_plt);
Eigenvalues_plt_imag = imag(Eigenvalues_plt);
Freq_plt = (ones(6,1)*Freq./Freq_Center)';

%
figure(1)
subplot(2,2,1)
hold on
plot(Eigenvalues_plt_real,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
% line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];  
for k = 1:1
    plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
%     plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.k')
end
% yline(1,'b','Linewidth',1.5)
% yline(Freq_RBE_Top/Freq_Center,'r','Linewidth',1.5)
% yline(Freq_RBE_Bottom/Freq_Center,'r','Linewidth',1.5)
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Re[kd/\pi]')
ylabel('Normalized Frequency')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');
%
subplot(2,2,2)
hold on
plot(Eigenvalues_plt_imag,Freq_plt,'.')
plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
for k = 1:1
    plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
%     plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.k')
end
% yline(1,'b','Linewidth',1.5)
% yline(Freq_RBE_Top/Freq_Center,'r','Linewidth',1.5)
% yline(Freq_RBE_Bottom/Freq_Center,'r','Linewidth',1.5)
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Im[kd/\pi]')
ylabel('Normalized Frequency')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');


subplot(2,2,3)
hold on
plot(Eigenvalues_plt_real,Eigenvalues_plt_imag,'.')
for k = 1:6
    plot(Eigenvalues_plt_real(:,k),Eigenvalues_plt_real(:,k),'.','Color',line_color(k,:))
%     plot(Eigenvalues_plt_real(:,k),Eigenvalues_plt_real(:,k),'.k')
end
axis([-1,1, -1,1])
xlabel('Re[kd/\pi]')
ylabel('Im[kd/\pi]')
grid on

% subplot(2,2,4)
% plot(Freq_plt(:,k),Sigma,'.')
% % axis([-1,1, -1,1])
% % xline(1,'b','Linewidth',1.5)
% % xline(Freq_RBE_Top/Freq_Center,'r','Linewidth',1.5)
% % xline(Freq_RBE_Bottom/Freq_Center,'r','Linewidth',1.5)
% ylabel('Coalescence parameter, \sigma')
% xlabel('Normalized Frequency')
% grid on
% ylabel('Sigma')

% Sigma_min = min(Sigma)
