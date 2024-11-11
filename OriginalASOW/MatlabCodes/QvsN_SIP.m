%% Calculate Q for N @ SIP

clear all; close all; clc
% Parameters SIP:
Freq_Center = 193.54;
CouplCoeff_k1 = 0.49832327234602;
Radius = 10e-6;
Alpha = 1.15240339832846;
Alpha_2 = 0.980630321583591;

% Constants
EffRefrIndex = 2.362;
SpeedLight = 3e8; % m/s

pA = (pi^2*EffRefrIndex*Radius*1e12)/SpeedLight;
pB = (4*pi*EffRefrIndex*Radius*Alpha*1e12)/SpeedLight;
pD = (4*pi*EffRefrIndex*Radius*Alpha_2*1e12)/SpeedLight;
% Frequency Sweep

Freqsteps = 5000;

% Diff_Freq = (Freq(end)-Freq(1))/Freqsteps; % As the Freq. vector is uniform, it's derivative is constant
NumberOfUnitCells =5:1:40; % Number of Unit Cells Used

PortValues = zeros(12,Freqsteps,length(NumberOfUnitCells));
for ii=1:length(NumberOfUnitCells)
    NumberOfUnitCells(ii)
    Freqmax = Freq_Center + 3e-3*(20/NumberOfUnitCells(ii))^3;%1e-3;
    Freqmin = Freq_Center - 3e-3*(20/NumberOfUnitCells(ii))^3;%8e-5*(40/NumberOfUnitCells)^2;
   Freq = linspace(Freqmin, Freqmax, Freqsteps);
    Freq = Freq';
    for jj = 1: Freqsteps
        [PortValues(:,jj,ii), TransferMatrices_Vec(:,:,jj)] = get_MyBrokenDownSOW_PortValues_ASOW_SIP_Parameters (NumberOfUnitCells(ii), Freq(jj), CouplCoeff_k1, Radius, Alpha, Alpha_2);
    end
end

%% Transfer Function & Q


TransferFunction = zeros(Freqsteps,length(NumberOfUnitCells));
PhaseTransferFunction = TransferFunction;
GroupDelay = zeros(Freqsteps-1,length(NumberOfUnitCells));
index = {};
Q = zeros(Freqsteps-1,length(NumberOfUnitCells));
% Q_0 = EffRefrIndex * 2*Radius .* NumberOfUnitCells;
GroupDelayBaseline = EffRefrIndex*NumberOfUnitCells*(2*pi*Radius + 2*(Alpha + Alpha_2)*Radius) / (SpeedLight);

for ii = 1:length(NumberOfUnitCells)
    NumberOfUnitCells(ii)
    Freqmax = Freq_Center + 3e-3*(20/NumberOfUnitCells(ii))^3;%1e-3;
    Freqmin = Freq_Center - 3e-3*(20/NumberOfUnitCells(ii))^3;%8e-5*(40/NumberOfUnitCells)^2;
    Freq = linspace(Freqmin, Freqmax, Freqsteps);
    Freq = 1e12*Freq';
    TransferFunction(:,ii) = PortValues(7,:,ii); % Input is the unity.
    TransferFunctiondB(:,ii) = 20*log10(abs(TransferFunction(:,ii)));
    PhaseTransferFunction(:,ii) = angle(TransferFunction(:,ii));
    GroupDelay(:,ii) = -diff(PhaseTransferFunction(:,ii))./((2*pi).*diff(Freq));
    [GroupDelay_min, index_min] = max(GroupDelay(:,ii));
    Freq_res(ii) = Freq(index_min);
    Q_unnormalized(:,ii) = pi*Freq_res(ii).*GroupDelay(:,ii);
end

Freq_plt = Freq./Freq_Center; % Frequency normalized to the SIP frequency
%% Check the peaks

% for ii = 1:length(NumberOfUnitCells)
%     figure(ii)
%     plot(Freq(2:end)./Freq_Center,Q_unnormalized(:,ii))
%     grid on
%     title('Q vs Freq')
% end

%% Plot of Q vs N
for ii = 1:length(NumberOfUnitCells)
    Q_peak_vec(ii) = max(Q_unnormalized(:,ii));
end

%% %% %%
close all
figure(4)
w_line = 2;
semilogy(NumberOfUnitCells, Q_peak_vec,'ob')
% loglog(NumberOfUnitCells, Q_peak_vec,'ob',NumberOfUnitCells, fit,'-k','Linewidth',1.5)
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
ylabel('$Q$','FontSize', 20,'Interpreter','latex');
xlabel('Odd number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
xlim([5 40]);
% ylim([5e5 2e7]);

%% Now we plot

% NumberOfUnitCells = 1:2:50;
% fit = 128.93*NumberOfUnitCells.^3 - 2.677e5; 
% figure(7)
% hold on
% % plot(NumberOfUnitCells, Q_peak_vec,'ob','Linewidth',3)
% plot(NumberOfUnitCells, 128.9.*NumberOfUnitCells.^3/(5354.*NumberOfUnitCells),'-k','Linewidth',1.5)
% % plot(NumberOfUnitCells, 115.8.*NumberOfUnitCells.^3/(7.729e4),'-b','Linewidth',1.5)
% 
% 
% hold off
% grid on
% legend('aN^3/b','aN^3/(+bN)')
% xlabel('Number of unit cells')
% ylabel('Fittings')

