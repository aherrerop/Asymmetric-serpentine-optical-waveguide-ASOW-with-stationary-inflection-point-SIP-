% DispersionDiagramDrawing
clear all; close all; clc;

% Parameters SIP:
% Freq_Center = 193.54;
% CouplCoeff_k1 = 0.49832327234602;
% Radius = 10e-6;
% Alpha = 1.15240339832846;
% Alpha_2 = 0.980630321583591;
Freq_Center = 193.54;
CouplCoeff_k1 = 0.49832327234602;
A = pi / 191.486975362555;
m = 4;
Radius = 10e-6;
Alpha = 1.15240339832846 - m*A;
Alpha_2 = 0.980630321583591 + m*A;

% Constants
EffRefrIndex = 2.362;
SpeedLight = 3e8; % m/s

FSR = (SpeedLight/EffRefrIndex)/(2*Radius*(pi+Alpha + Alpha_2));
FSR_Freq_Center = FSR/(Freq_Center*1e12);

pA = (pi^2*EffRefrIndex*Radius*1e12)/SpeedLight;% +(4*pi*EffRefrIndex*Radius*(2*A)*1e12)/SpeedLight ;
pB = (4*pi*EffRefrIndex*Radius*(Alpha)*1e12)/SpeedLight;
pD = (4*pi*EffRefrIndex*Radius*(Alpha_2)*1e12)/SpeedLight;

% Frequency Sweep
Freqmax = Freq_Center + 0.04;
Freqmin = Freq_Center - 0.04;
Freqsteps = 5000;

Freq = linspace(Freqmin, Freqmax, Freqsteps);
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);
NumberOfUnitCells = 1;
for ii = 1:Freqsteps    
    TransferMatrices_Vec(:,:,ii) =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    Sigma(ii) = SIP_Check_Eigenvectors(TransferMatrices_Vec(:,:,ii));
end

% plot dispersion
    
Eigenvalues_plt = Eigenvalues'/pi;
Eigenvalues_plt_real = real(Eigenvalues_plt);
Eigenvalues_plt_imag = imag(Eigenvalues_plt);
Freq_plt = (ones(6,1)*Freq./Freq_Center)';

Eigenvalues_plt(400,:)';
Eigenvalues_plt(4600,:)';
Eigenvalues_plt(2500,:)';

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
end
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Re[kd/\pi]')
ylabel('Normalized Frequency')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

subplot(2,2,2)
hold on
plot(Eigenvalues_plt_imag,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
for k = 1:1
    plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
end
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
end
axis([-1,1, -1,1])
xlabel('Re[kd/\pi]')
ylabel('Im[kd/\pi]')
grid on

subplot(2,2,4)
plot(Freq_plt(:,k),Sigma,'.')
% axis([-1,1, -1,1])
ylabel('Normalized Frequency')
grid on
ylabel('Sigma')

% Sigma_min = min(Sigma)
