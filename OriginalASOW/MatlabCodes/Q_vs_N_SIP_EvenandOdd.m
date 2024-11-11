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
Freqmax = Freq_Center + 1e-4;
Freqmin = Freq_Center - 1e-4;
Freqsteps = 33;

Freq = linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';
Diff_Freq = (Freq(end)-Freq(1))/Freqsteps; % As the Freq. vector is uniform, it's derivative is constant
NumberOfUnitCells = 2:2:50; % Number of Unit Cells Used
oNumberOfUnitCells = 1:2:50; % Odd

PortValues = zeros(12,Freqsteps,length(NumberOfUnitCells));
for ii=1:length(NumberOfUnitCells)
    for jj = 1: Freqsteps
        [PortValues(:,jj,ii), SystemMatrix] = get_MyBrokenDownSOW_PortValues_3_PhaseDelay_1CouplCoeff (NumberOfUnitCells(ii),Freq(jj),CouplCoeff_k1, pA, pB, pD);
        [oPortValues(:,jj,ii), SystemMatrix] = get_MyBrokenDownSOW_PortValues_3_PhaseDelay_1CouplCoeff (oNumberOfUnitCells(ii),Freq(jj),CouplCoeff_k1, pA, pB, pD);
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
    TransferFunction(:,ii) = PortValues(7,:,ii); % Input is the unity.
    TransferFunctiondB(:,ii) = 20*log10(abs(TransferFunction(:,ii)));
    PhaseTransferFunction(:,ii) = angle(TransferFunction(:,ii));
    GroupDelay(:,ii) = -diff(PhaseTransferFunction(:,ii))./((2*pi).*diff(Freq));
    [GroupDelay_min, index_min] = max(GroupDelay(:,ii));
    Freq_res(ii) = Freq(index_min);
    Q_unnormalized(:,ii) = pi*Freq_res(ii).*GroupDelay(:,ii);
 
    oTransferFunction(:,ii) = oPortValues(7,:,ii); % Input is the unity.
    oTransferFunctiondB(:,ii) = 20*log10(abs(oTransferFunction(:,ii)));
    oPhaseTransferFunction(:,ii) = angle(oTransferFunction(:,ii));
    oGroupDelay(:,ii) = -diff(oPhaseTransferFunction(:,ii))./((2*pi).*diff(Freq));
    [oGroupDelay_min, oindex_min] = max(oGroupDelay(:,ii));
    oFreq_res(ii) = Freq(oindex_min);
    oQ_unnormalized(:,ii) = pi*oFreq_res(ii).*oGroupDelay(:,ii);

    Q_normalization(ii) = 2*pi*Freq_res(ii)*GroupDelayBaseline(ii);
    Q_normalized(:,ii) = Q_unnormalized(:,ii) ./ Q_normalization(ii);
end

Freq_plt = Freq./Freq_Center; % Frequency normalized to the SIP frequency

% %% Check the peaks
% figure(1)
% hold on
% for ii = 1:length(NumberOfUnitCells)
%     plot(Freq(2:end)./Freq_Center,Q_unnormalized(:,ii))
% end
% hold off
% grid on
% title('Q vs Freq')

%% Plot of Q vs N
for ii = 1:length(NumberOfUnitCells)
    Q_peak_vec(ii) = max(Q_unnormalized(:,ii));
    oQ_peak_vec(ii) = max(oQ_unnormalized(:,ii));
end

%% Normalized Q plot

figure(4)
hold on
w_line = 2;
hold on
plot(NumberOfUnitCells, Q_peak_vec,'ob','linewidth',w_line)
plot(oNumberOfUnitCells, oQ_peak_vec,'ok','linewidth',w_line)
plot(NumberOfUnitCells, Q_peak_vec,'-b','linewidth',w_line)
plot(oNumberOfUnitCells, oQ_peak_vec,'-k','linewidth',w_line)
hold off
grid on
set(gca,'FontSize',28,'FontName', 'Times New Roman');
ylabel('$Q$','FontSize', 28,'Interpreter','latex');
xlabel('Number of unit cells, $N$','FontSize', 28,'Interpreter','latex')
legend('Even N','Odd N');
% pbaspect([1.3 1 1]);
xlim([0 max(NumberOfUnitCells)]);
% ylim([0 3e7])

%% In dB
figure(5)
w_line = 2;
hold on
semilogx(NumberOfUnitCells, Q_peak_vec,'ob','linewidth',w_line)
semilogx(oNumberOfUnitCells, oQ_peak_vec,'ok','linewidth',w_line)
semilogx(NumberOfUnitCells, Q_peak_vec,'-b','linewidth',w_line)
semilogx(oNumberOfUnitCells, oQ_peak_vec,'-k','linewidth',w_line)
hold off
grid on
set(gca,'FontSize',28,'FontName', 'Times New Roman');
ylabel('$Q$','FontSize', 28,'Interpreter','latex');
xlabel('Number of unit cells, $N$','FontSize', 28,'Interpreter','latex')
legend('Even N','Odd N')
% pbaspect([1.3 1 1]);
xlim([0 max(NumberOfUnitCells)]);
% ylim([0 15e6])