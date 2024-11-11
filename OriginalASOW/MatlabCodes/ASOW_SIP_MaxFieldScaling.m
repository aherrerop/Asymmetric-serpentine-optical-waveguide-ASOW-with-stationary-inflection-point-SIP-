%% ASOW 6DBE maximum field amplitude scaling
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

NumberOfUnitCells = [1:1:50]; % Number of Unit Cells Used

% NumberOfUnitCells = linspace(1,128,128);
PortValues = zeros(12,1);
for jj = 1:length(NumberOfUnitCells)

    [PortValues, SystemMatrix, TransferMatrix, T_aux] = get_MyBrokenDownSOW_PortValues_3_PhaseDelay_1CouplCoeff (NumberOfUnitCells(jj),Freq_Center,CouplCoeff_k1, pA, pB, pD);
    % Values at each interface
    E_0 = PortValues(1:6);
    E_End = PortValues(7:12);
    
    % Now let's find the values in the middle
    % Forwards
    E_f = zeros(6,NumberOfUnitCells(jj)+1);
    E_f(:,1) = E_0;
    % E_f(NumberOfUnitCells,:) = E_End;
    for ii = 2:NumberOfUnitCells(jj)
        E_f(:,ii)= TransferMatrix*E_f(:,ii-1);
    end
    E_f(:,NumberOfUnitCells(jj)+1) = T_aux*E_f(:,NumberOfUnitCells(jj));
    % Backwards
    E_b = zeros(6,NumberOfUnitCells(jj)+1);
    E_b(:,end) = E_End;
    E_b(:,end-1) = T_aux^-1*E_b(:,end);
    for ii = 2:NumberOfUnitCells(jj)
        E_b(:,end-ii) = TransferMatrix^-1*E_b(:,end-ii+1);
    end
    
%     if abs(E_b-E_f)<1e-2
    E_1 = abs(E_f(1,:) + E_f(2,:));
    E_1b = abs(E_b(1,:) + E_b(2,:));
    MaxField(jj) = max(E_1);
%     end
    % Plotting for each NumberOfUnitCells(jj) to check the maximum
%     UnitCell_vec = linspace(0,NumberOfUnitCells(jj), NumberOfUnitCells(jj)+1); 
%     w_line=1.5;
%     figure(jj)
%     hold on
%     plot(UnitCell_vec, abs(E_1),'ok','linewidth',w_line)
%     plot(UnitCell_vec, abs(E_1b),'-b','linewidth',1)
%     hold off
%     legend('Calculated forward','Calculated backward')
%     xlabel('Number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
%     ylabel('$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
%     grid on

end
MaxIntensity = abs(MaxField).^2;
%
%% Fitting
UnitCell_vec = linspace(0,NumberOfUnitCells(end), NumberOfUnitCells(end)+1); 

%% Separate even and odd behavior:
Even_NumberOfUnitCells = NumberOfUnitCells(2:2:end);
Odd_NumberOfUnitCells = NumberOfUnitCells(1:2:end);
Even_MaxIntensity = MaxIntensity(2:2:end);
Odd_MaxIntensity = MaxIntensity(1:2:end);
Even_MaxField = MaxField(2:2:end);
Odd_MaxField = MaxField(1:2:end);



%% Fitting
x = Even_NumberOfUnitCells(1,:)';
y = Even_MaxField(1,:)';
ft = fittype('a7*x^7 + a6*x^6 + a5*x^5 + a4*x^(4) + a3*x^3 + a2*x^2 + a1*x + a0');
[cf, g] = fit(x,y,ft)

% 
% Even_Field_fit = 0.0004366*NumberOfUnitCells.^4-0.01071*NumberOfUnitCells.^3 +0.1458*NumberOfUnitCells.^2-0.61*NumberOfUnitCells +1.243;
% Odd_Field_fit = 0.0024*NumberOfUnitCells.^4 -0.06*NumberOfUnitCells.^3+0.6*NumberOfUnitCells.^2 -1.365*NumberOfUnitCells -0.7499;
% Even_Intensity_fit = 0.000158*NumberOfUnitCells.^4 +0.126*NumberOfUnitCells.^2 +1.094*NumberOfUnitCells -6.401;
% Odd_Intensity_fit = 0.02*NumberOfUnitCells.^4 -0.52*NumberOfUnitCells.^3 + 4.652*NumberOfUnitCells.^2-16.22*NumberOfUnitCells +15.05;
%% Plot with fitting
% w_line = 2;
% figure(100)
% hold on
% plot(Even_NumberOfUnitCells,Even_MaxIntensity,'ok','linewidth',w_line)
% plot(NumberOfUnitCells,Even_Intensity_fit,'r','linewidth',w_line-1)
% plot(Odd_NumberOfUnitCells,Odd_MaxIntensity,'ob','linewidth',w_line)
% % plot(Odd_NumberOfUnitCells,Odd_Intensity_fit,'r','linewidth',w_line)
% hold off
% legend('Even N','Even fit','Odd N')
% xlabel('Number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
% ylabel('Max$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
% grid on
%%
w_line = 2;
figure(101)
hold on
plot(Even_NumberOfUnitCells,Even_MaxField,'ok','linewidth',w_line)
plot(NumberOfUnitCells,Field_fit,'r','linewidth',w_line-1)
plot(Odd_NumberOfUnitCells,Odd_MaxField,'ob','linewidth',w_line)
% plot(Odd_NumberOfUnitCells,Odd_Intensity_fit,'r','linewidth',w_line)
hold off
legend('Even N','Even fit','Odd N')
xlabel('Number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
ylabel('Max$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
grid on

%% Plotting
% 
% 
% w_line = 2;
% figure(1)
% hold on
% plot(NumberOfUnitCells,MaxIntensity,'ok','linewidth',w_line)
% % plot(NumberOfUnitCells,MaxIntensity,'b','linewidth',w_line-1)
% plot(NumberOfUnitCells,Intensity_fit,'r','linewidth',w_line-1)
% hold off
% xlabel('Number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
% ylabel('Max$|E_1(n)|^2  /  |E_{inc}|^2$','FontSize', 20,'Interpreter','latex')
% grid on
% % 
% % 
% figure(2)
% hold on
% plot(NumberOfUnitCells,MaxField,'ok','linewidth',w_line)
% % plot(NumberOfUnitCells,MaxField,'b','linewidth',w_line-1)
% plot(NumberOfUnitCells,Field_fit,'r','linewidth',w_line-1)
% hold off
% legend('Data','fit')
% xlabel('Number of unit cells, $N$','FontSize', 20,'Interpreter','latex')
% ylabel('Max$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
% grid on
%%
% 
% % figure(3)
% % hold on
% % plot(UnitCell_vec,abs(E_f(:,1) + E_b(:,2)),'ok','linewidth',w_line)
% % plot(UnitCell_vec,abs(E_f(:,1) + E_b(:,2)),'b','linewidth',w_line-1)
% % hold off
% % grid on
% 
