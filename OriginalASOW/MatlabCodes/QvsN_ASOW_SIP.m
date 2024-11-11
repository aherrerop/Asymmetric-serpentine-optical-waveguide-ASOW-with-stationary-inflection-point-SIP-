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
Freqsteps = 5000;

Freq = linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';
Diff_Freq = (Freq(end)-Freq(1))/Freqsteps; % As the Freq. vector is uniform, it's derivative is constant
NumberOfUnitCells = 5:1:40; % Number of Unit Cells Used

PortValues = zeros(12,Freqsteps,length(NumberOfUnitCells));
for ii=1:length(NumberOfUnitCells)
    for jj = 1: Freqsteps
        [PortValues(:,jj,ii), TransferMatrices_Vec(:,:,jj)] = ASOW_Gain_PortValues (NumberOfUnitCells(ii),Freq(jj), CouplCoeff_k1, Radius, Alpha, Alpha_2, EffRefrIndex);
    end
end

%% Transfer Function & Q
Freq_Center = 1e12*Freq_Center;
Freq = 1e12*linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';

TransferFunction = zeros(Freqsteps,length(NumberOfUnitCells));
PhaseTransferFunction = zeros(Freqsteps,length(NumberOfUnitCells));
% The group delay has one less element bc 
GroupDelay = zeros(Freqsteps-1,length(NumberOfUnitCells));
Q = zeros(Freqsteps-1,length(NumberOfUnitCells));
% Define the vector tau_0 (group delay going through a straight waveguide
% of the same length as each finite length structure).
GroupDelayBaseline = EffRefrIndex*(NumberOfUnitCells)*(2*pi*Radius + 2*(Alpha + Alpha_2)*Radius) / (SpeedLight);

for ii = 1:length(NumberOfUnitCells)
    TransferFunction(:,ii) = PortValues(7,:,ii); % The incident field amplitude is the unity.
    TransferFunctiondB(:,ii) = 20*log10(abs(TransferFunction(:,ii))); % Take the magnitude in dB
    PhaseTransferFunction(:,ii) = angle(TransferFunction(:,ii)); % Take the phase of the transfer function
    % The group delay is the derivate of the phase of the transfer function
    % wrt the angular frequency (hence the 2pi dividing). We will add the
    % sign later.
    GroupDelay(:,ii) = diff(PhaseTransferFunction(:,ii))./((2*pi).*diff(Freq));
    [GroupDelay_min, index_min] = min(GroupDelay(:,ii)); % Because the group delay is negative now. The sign will be added later.
    GroupDelay_norm(:,ii) = GroupDelay(:,ii) / GroupDelayBaseline(ii);  % Normalize the group delay
    Freq_res(ii) = Freq(index_min); % We get the frequency of the highest peak of the group delay (lowest bc we have yet to change the sign)
    % We can use different definitions of Q. We use the unnormalized one.
    Q(:,ii) = pi*Freq_res(ii).*GroupDelay_norm(:,ii); % Q with normalized group delay
    Q_unnormalized(:,ii) = pi*Freq_res(ii).*GroupDelay(:,ii); % Q unnormalized
    Q_normalization(ii) = pi*Freq_res(ii)*GroupDelayBaseline(ii); % normalization Q
    Q_normalized(:,ii) = Q_unnormalized(:,ii) ./ Q_normalization(ii); % Q normalized 2.
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
    Q_peak_vec(ii) = -min(Q_unnormalized(:,ii));
%     Q_peak_norm_vec(ii) = max(Q_norm(:,ii));
end
w_line = 2;
figure(4)
hold on
plot(Freq(2:end)./Freq_Center,-GroupDelay_norm(:,1),'-k','linewidth',w_line)
plot(Freq(2:end)./Freq_Center,-GroupDelay_norm(:,2),'-b','linewidth',w_line)
plot(Freq(2:end)./Freq_Center,-GroupDelay_norm(:,3),'-','linewidth',w_line)
% plot(Freq_Q_plt_4,-GroupDelay_plt_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,-GroupDelay_plt_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,-GroupDelay_plt_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% legend('$N$=8', '$N$=12','$N$=16','FontSize', 20,'Interpreter','latex');
% title('Normalized $\tau_g$ versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('$\tau_g / \tau_{0}$','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');

grid on



w_line = 2;
figure(5)
hold on
plot(Freq(2:end)./Freq_Center,-Q_unnormalized(:,1),'-k','linewidth',w_line)
plot(Freq(2:end)./Freq_Center,-Q_unnormalized(:,2),'-b','linewidth',w_line)
plot(Freq(2:end)./Freq_Center,-Q_unnormalized(:,3),'-','linewidth',w_line)
% plot(Freq_Q_plt_4,-GroupDelay_plt_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,-GroupDelay_plt_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,-GroupDelay_plt_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% legend('$N$=8', '$N$=12','$N$=16','FontSize', 20,'Interpreter','latex');
% title('Normalized $\tau_g$ versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('$Q$','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');

grid on

%% Fit
% f=polyfit(NumberOfUnitCells(5:end),Q_peak_vec(5:end),3);
% fit = f(1).*NumberOfUnitCells.^3 + f(2).*NumberOfUnitCells.^2 + f(3).*NumberOfUnitCells + f(4);
late_start = 1; %First position from which we will do the fitting. Up to 25.
x = NumberOfUnitCells(late_start:end)';
y = Q_peak_vec(late_start:end)';
% Copy the vectors to another file and run the next comments. "f" has the
% coefficients you need to apply to get "fit".
% ft = fittype('a*x.^3+b');
% options = fitoptions(ft);
% [f,g] = fit(x,y,ft,options);
   
% fit = 128.9*NumberOfUnitCells.^3 -5354; %even fitting
fit = 99.8*NumberOfUnitCells.^3 + 3.2e4*NumberOfUnitCells-3.4e5; % odd fitting

% Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)
Rsq2 = 1 - sum((Q_peak_vec - fit).^2)/sum((Q_peak_vec - mean(Q_peak_vec)).^2);

%%
figure(10)
w_line = 2;
semilogy(NumberOfUnitCells, -Q_peak_vec,'ob',NumberOfUnitCells, fit,'-k','Linewidth',1.5)
% loglog(NumberOfUnitCells, -Q_peak_vec,'ob',NumberOfUnitCells, fit,'-k','Linewidth',1.5)
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
ylabel('$Q$','FontSize', 20,'Interpreter','latex');
xlabel('Odd number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
xlim([20 50]);
ylim([5e5 2e7])

% % %%
% figure(5)
% hold on
% w_line = 2;
% plot(NumberOfUnitCells, 20*log(Q_peak_vec),'ob','Linewidth',3)
% plot(NumberOfUnitCells, 20*log(fit),'-k','Linewidth',1.5)
% plot(NumberOfUnitCells, 20*log(fitx),'-r','Linewidth',1.5)
% plot(NumberOfUnitCells, 20*log(fitx2),'-','Linewidth',1.5)
% plot(NumberOfUnitCells, 20*log(fitx4),'-','Linewidth',1.5)
% plot(NumberOfUnitCells, 20*log(fitxall),'-','Linewidth',1.5)
% hold off
% legend('Q','Q \approx 93.49N^3+1.031 \times 10^6, R^2=0.965','Q \approx 15.49N^3+3.104 \times 10^5N-5.878 \times 10^6, R^2 = 0.9955','Q \approx -74.45N^3+9332N^2-2.58 \times 10^6, R^2 = 0.99689','Q \approx 1.791N^4+2.084 \times 10^6, R^2 = 0.925','Q \approx -355.5N^3+3.915 \times 10^4N^2-1.017 \times 10^6N+8.52 \times 10^6, R^2 = 0.9988')
% % legend('data','Q \approx aN^3+b, R^2=0.965','Q \approx aN^3+bN+c, R^2 = 0.9955','Q \approx aN^3+bN^2+d, R^2 = 0.99689','Q \approx aN^4+b, R^2 = 0.925','Q \approx aN^3+bN^2+cN+d, R^2 = 0.9988')
% grid on
% set(gca,'FontSize',28,'FontName', 'Times New Roman');
% ylabel('$20\log{(Q)}$','FontSize', 28,'Interpreter','latex');
% xlabel('Odd number of unit cells, $N$','FontSize', 28,'Interpreter','latex')
% xlim([20 50]);
%% Now we plot
% 
% figure(5)
% hold on
% w_line = 2;
% loglog(NumberOfUnitCells, Qdb,'o','linewidth',w_line)
% % loglog(NumberOfUnitCells, 0.0009673*NumberOfUnitCells.^3,'-','linewidth',w_line)
% hold off
% grid on
% set(gca,'FontSize',28,'FontName', 'Times New Roman');
% ylabel('$Q$','FontSize', 28,'Interpreter','latex');
% xlabel('Even number of unit cells, $N$','FontSize', 28,'Interpreter','latex')
% ylim([0 3e7])
% pbaspect([1.3 1 1]);
% xlim([min(NumberOfUnitCells) max(NumberOfUnitCells)]);

% %% Only fitting
% NumberOfUnitCells_2 = linspace(5, 100, 95);
% fit_2 = f(1).*NumberOfUnitCells_2.^3 + f(2).*NumberOfUnitCells_2.^2 + f(3).*NumberOfUnitCells_2 + f(4);
% 
% figure(2)
% hold on
% plot(NumberOfUnitCells_2, fit_2,'-k','Linewidth',1.5)
% plot(NumberOfUnitCells, Q_peak_vec,'ob','Linewidth',1.5)
% hold off
% grid on