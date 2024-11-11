%% Calculate Transfer Function and Q factor at a SIP
% Normalized Group Delay & Q.
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

% pA = (pi^2*EffRefrIndex*Radius*1e12)/SpeedLight;
% pB = (4*pi*EffRefrIndex*Radius*Alpha*1e12)/SpeedLight;
% pD = (4*pi*EffRefrIndex*Radius*Alpha_2*1e12)/SpeedLight;
% Frequency Sweep
Delta = 0.04;
Freqmax = Freq_Center + Delta;
Freqmin = Freq_Center - Delta;
Freqsteps = 5000;
% Frequency Vector
Freq = linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';
% Length of the finite length structure
NumberOfUnitCells = [8 12 16]; % Number of Unit Cells Used
% Vector with the field amplitudes at either side for all frequencies &
% lengths of the finite-length structure.
PortValues = zeros(12,Freqsteps,length(NumberOfUnitCells));
for ii=1:length(NumberOfUnitCells)
    for jj = 1: Freqsteps
        % This function generates the 12x12 system matrix with the finite
        % length structure transfer matrix and the boundary condition.
%         [PortValues(:,jj,ii), SystemMatrix] = get_MyBrokenDownSOW_PortValues_3_PhaseDelay_1CouplCoeff (NumberOfUnitCells(ii),Freq(jj),CouplCoeff_k1, pA, pB, pD);
        [PortValues(:,jj,ii), TransferMatrices_Vec(:,:,jj)] = get_MyBrokenDownSOW_PortValues_ASOW_SIP_Parameters (NumberOfUnitCells(ii), Freq(jj), CouplCoeff_k1, Radius, Alpha, Alpha_2);
    end
end
% Let's assign the righ units for the Freq_Center.
Freq_Center = 1e12*Freq_Center;
Freq = 1e12*linspace(Freqmin, Freqmax, Freqsteps);
Freq = Freq';

% Transfer Function & Q
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
    
    ReflectionFunction(:,ii) = PortValues(2,:,ii); % The incident field amplitude is the unity.
    ReflectionFunctiondB(:,ii) = 20*log10(abs(ReflectionFunction(:,ii))); % Take the magnitude in dB
    PhaseReflectionFunction(:,ii) = angle(ReflectionFunction(:,ii)); % Take the phase of the Reflection function
end

Freq_plt = Freq./Freq_Center; % Frequency normalized to the SIP frequency
% Lossless check:
% Want to check that |S_11|^2 + |S_21|^2 = 1
ShouldBeOne = abs(TransferFunction).^2 + abs(ReflectionFunction).^2;

%
figure(1)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, TransferFunctiondB(:,ii),'-','linewidth',w_line)
end
hold off
grid on
% pbaspect([1.3 1 1]);
legend('$N$=8', '$N$=12','$N$=16','FontSize', 24,'Interpreter','latex');
set(gca,'FontSize',24,'FontName', 'Times New Roman');
% title('Transfer Function [dB] versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('$\mid$ Transfer function $\mid$ (dB)','FontSize', 24,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 24,'Interpreter','latex');
ylim([-40 0])
xlim([1e12*Freqmin/Freq_Center 1e12*Freqmax/Freq_Center])

%%
figure(2)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, PhaseTransferFunction(:,ii),'-','linewidth',w_line)
end
hold off
pbaspect([1.3 1 1]);
legend('$N$=8', '$N$=9','$N$=10','FontSize', 20,'Interpreter','latex');
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% title('Transfer Function [dB] versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('Angle Transfer functions','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
grid on
%%

figure(3)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, ReflectionFunctiondB(:,ii),'-','linewidth',w_line)
end
hold off
% pbaspect([1.3 1 1]);
legend('$N$=8', '$N$=9','$N$=10','FontSize', 20,'Interpreter','latex');
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% title('Reflection Function [dB] versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
% ylabel('$\mid$ Reflection function $\mid$ (dB)','FontSize', 20,'Interpreter','latex');
% xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
% ylim([-15 0])
xlim([1-0.00001 1+0.00001])
% axis([1-0.00001 1+0.00001 -15 0])
grid on
%%
figure(4)
hold on
w_line = 2;
for ii = 1:length(NumberOfUnitCells)
    plot(Freq_plt, PhaseReflectionFunction(:,ii),'-','linewidth',w_line)
end
hold off
pbaspect([1.3 1 1]);
legend('$N$=8', '$N$=9','$N$=10','FontSize', 20,'Interpreter','latex');
set(gca,'FontSize',20,'FontName', 'Times New Roman');
% title('Reflection Function [dB] versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('Angle Reflection functions','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
grid on
%%
% Filter the discontinuities in the Angle Transfer Function,
% which causes discontinuities in the Group Delay:
index_1 = find(GroupDelay_norm(:,1)<0);
index_12 = find(GroupDelay_norm(:,1)<(-200));
index_1 = setdiff(index_1,index_12);
Freq_Q_plt_1 = Freq_plt(index_1);
GroupDelay_plt_1 = GroupDelay_norm(index_1,1);
Q_norm_plt_1 = Q_normalized(index_1,1);
Q_plt_1 = Q_unnormalized(index_1,1);
   
index_2 = find(GroupDelay_norm(:,2)<0);
index_21 = find(GroupDelay_norm(:,2)<(-200));
index_2 = setdiff(index_2,index_21);
Freq_Q_plt_2 = Freq_plt(index_2);
GroupDelay_plt_2 = GroupDelay_norm(index_2,2);
Q_norm_plt_2 = Q_normalized(index_2,2);
Q_plt_2 = Q_unnormalized(index_2,2);

index_3 = find(GroupDelay_norm(:,3)<0);
index_32 = find(GroupDelay_norm(:,3)<(-200));
index_3 = setdiff(index_3,index_32);
Freq_Q_plt_3 = Freq_plt(index_3);
GroupDelay_plt_3 = GroupDelay_norm(index_3,3);
Q_norm_plt_3 = Q_normalized(index_3,3);
Q_plt_3 = Q_unnormalized(index_3,3);

% index_4 = find(GroupDelay(:,4)<0);
% Freq_Q_plt_4 = Freq_plt(index_4);
% GroupDelay_plt_4 = GroupDelay_norm(index_4,4);
% Q_norm_plt_4 = Q_normalized(index_4,4);
% 
% index_5 = find(GroupDelay(:,5)<0);
% Freq_Q_plt_5 = Freq_plt(index_5);
% GroupDelay_plt_5 = GroupDelay_norm(index_5,5);
% Q_norm_plt_5 = Q_normalized(index_5,5);
% 
% index_6 = find(GroupDelay(:,6)<0);
% Freq_Q_plt_6 = Freq_plt(index_6);
% GroupDelay_plt_6 = GroupDelay_norm(index_6,6);
% Q_norm_plt_6 = Q_normalized(index_6,6);

% Q_1 = pi*Freq_res(1).*abs(GroupDelay_plt_1);
% Q_2 = pi*Freq_res(2).*abs(GroupDelay_plt_2);
% Q_3 = pi*Freq_res(3).*abs(GroupDelay_plt_3);
% Q_4 = pi*Freq_res(4).*abs(GroupDelay_plt_4);
% Q_5 = pi*Freq_res(5).*abs(GroupDelay_plt_5);
% Q_6 = pi*Freq_res(6).*abs(GroupDelay_plt_6);




%%
w_line = 2;
figure(5)
hold on
plot(Freq_Q_plt_1,-GroupDelay_plt_1,'-','linewidth',w_line)
plot(Freq_Q_plt_2,-GroupDelay_plt_2,'-','linewidth',w_line)
plot(Freq_Q_plt_3,-GroupDelay_plt_3,'-','linewidth',w_line)
% plot(Freq_Q_plt_4,-GroupDelay_plt_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,-GroupDelay_plt_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,-GroupDelay_plt_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',24,'FontName', 'Times New Roman');
legend('$N$=8', '$N$=9','$N$=10','FontSize', 24,'Interpreter','latex');
% title('Normalized $\tau_g$ versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
grid on
ylabel('$\tau_g / \tau_{0}$','FontSize', 24,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 24,'Interpreter','latex');
xlim([1e12*Freqmin/Freq_Center 1e12*Freqmax/Freq_Center])
% xlim([1-2e-4 1+2e-4])
% pbaspect([1.3 1 1]);

figure(6)
hold on
plot(Freq_Q_plt_1,-Q_plt_1,'-','linewidth',w_line)
plot(Freq_Q_plt_2,-Q_plt_2,'-','linewidth',w_line)
plot(Freq_Q_plt_3,-Q_plt_3,'-','linewidth',w_line)
% plot(Freq_Q_plt_4,Q_4,'-','linewidth',w_line)
% plot(Freq_Q_plt_5,Q_5,'-','linewidth',w_line)
% plot(Freq_Q_plt_6,Q_6,'-','linewidth',w_line)
hold off
set(gca,'FontSize',20,'FontName', 'Times New Roman');
legend('$N$=8', '$N$=9','$N$=10','FontSize', 20,'Interpreter','latex');
grid on
% title('Normalized $Q$ versus $\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
ylabel('$Q$','FontSize', 20,'Interpreter','latex');
xlabel('$\omega / \omega_s$','FontSize', 20,'Interpreter','latex');
% pbaspect([1.3 1 1]);
% axis([Freqmin*1e12/Freq_Center Freqmax*1e12/Freq_Center 0 3.2e5])
% xlim([1-2e-4 1+2e-4])
