%% Field Amplitudes at every cell for a finite length structure
%
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

NumberOfUnitCells = 31; % Number of Unit Cells Used

% Let's find the values at each side:
PortValues = zeros(12,1);

[PortValues, SystemMatrix, TransferMatrix, T_aux] = get_MyBrokenDownSOW_PortValues_3_PhaseDelay_1CouplCoeff (NumberOfUnitCells,Freq_Center,CouplCoeff_k1, pA, pB, pD);
% TransferMatrix2 =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq_Center, CouplCoeff_k1, Radius, Alpha, Alpha_2);

E_0 = PortValues(1:6);
E_End = PortValues(7:12);

% Now let's find the values in the middle
% Forwards
E_f = zeros(6,NumberOfUnitCells+1);
E_f(:,1) = E_0;
% E_f(NumberOfUnitCells,:) = E_End;
for ii = 2:NumberOfUnitCells
    E_f(:,ii)= TransferMatrix*E_f(:,ii-1);
end
E_f(:,NumberOfUnitCells+1) = T_aux*E_f(:,ii);

% Backwards
E_b = zeros(6,NumberOfUnitCells+1);
E_b(:,end) = E_End;
E_b(:,end-1) = T_aux^-1*E_b(:,end);
for ii = 2:NumberOfUnitCells
    E_b(:,end-ii) = TransferMatrix^-1*E_b(:,end-ii+1);
end

UnitCell_vec = linspace(0,NumberOfUnitCells, NumberOfUnitCells+1); 
%% Let's plot
w_line = 2;
w_line = 2;

figure(6)
hold on
plot(UnitCell_vec,abs(E_f(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(E_b(1,:)),'-b','linewidth',1)
hold off
% title('Normalized $|E_1^+|$ vs unit cell number, $n$','FontSize', 28,'Interpreter','latex')
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1^+(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 12])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%%
figure(7)
% plot(UnitCell_vec, abs(E_backwards(:,2)),'-','linewidth',w_line)
hold on
plot(UnitCell_vec,abs(E_f(2,:)),'ok','linewidth',w_line)
plot(UnitCell_vec,abs(E_b(2,:)),'-b','linewidth',1)
hold off
% title('Normalized $|E_1^-|$ vs unit cell number, $n$','Interpreter','latex')
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1^-(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
axis([0 NumberOfUnitCells 0 10])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
%%
figure(8)
% plot(UnitCell_vec, abs(E_backwards(:,2) + E_forwards(:,1)),'-','linewidth',w_line)
hold on
plot(UnitCell_vec, abs(E_f(2,:) + E_b(1,:)),'ok','linewidth',w_line)
plot(UnitCell_vec, abs(E_f(2,:) + E_b(1,:)),'-b','linewidth',1)
hold off
% title('Normalized $|E_1|$ vs unit cell number, $n$','Interpreter','latex')
xlabel('Unit cell number, $n$','FontSize', 20,'Interpreter','latex')
ylabel('$|E_1(n)|  /  |E_{inc}|$','FontSize', 20,'Interpreter','latex')
% axis([0 NumberOfUnitCells 0 6])
grid on
set(gca,'FontSize',20,'FontName', 'Times New Roman');
