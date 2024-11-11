%% Simplify it
clear all; close all; clc
% Start the loop:
Sigma_EPD = [];
Eigenvalues_EPD = [];
X_EPD = [];
Ksd_pi_EPD = [];
Eigenvectors_EPD = [];
NumberOfIterations = 10;
for ii = 1:NumberOfIterations
    
    % Set a CouplCoeff_min & a CouplCoeff_max.
    CouplCoeff_min = 0;
    CouplCoeff_max = 1;
    % CouplCoeff_alpha_seed between CouplCoeff_min & CouplCoeff_max.
    CouplCoeff_k1_seed= CouplCoeff_min + rand()*(CouplCoeff_max-CouplCoeff_min);
%     CouplCoeff_k2_seed = CouplCoeff_min + rand()*(CouplCoeff_max-CouplCoeff_min);
    
    % Frequency:
    Freq_Center = 193.54;
    
    % Phases:
    Alpha_min = 1;
    Alpha_max = 1.035;
    %
    Alpha_seed = Alpha_min + rand()*(Alpha_max-Alpha_min);
    Alpha_2_seed = Alpha_min + rand()*(Alpha_max-Alpha_min);

    
    % Initial guess (randomly assigned between X0_min & X0_max).
    X0 = [CouplCoeff_k1_seed, Alpha_seed, Alpha_2_seed];
    opts = optimset('MaxIter',1e8, 'Display','none');
    % Gives X at the SIP. Sigma is the eigenvector coalescence parameter (the more similar the
    % eigenvalues, the smaller it is (tends to zero)).
    [X,sigma, exitflag] = fminsearch(@Opt_Func_3PhaseDelay_1CouplCoeff_Parameters, X0, opts);
    if exitflag>0 && X(3) < 1.035 && X(3) > 1 && X(2) < 1.035 && X(2) > 1 
        Freq_Center = 193.54;
        CouplCoeff_k1 = X(1);
        Alpha = X(2);
        Alpha_2 = X(3);
        Radius = 10e-6;

        
        [TransferMatrix] = get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq_Center, CouplCoeff_k1, Radius, Alpha, Alpha_2);
        [Eigenvectors, EigenValues_X] = eig(TransferMatrix);
        Eigenvalues = diag(EigenValues_X);
        [Eigenvalues_X,index_X] = sort(Eigenvalues);
        Ksd_X = -log(Eigenvalues_X)./(1j);
        Ksd_pi = Ksd_X/pi;
        Eigenvectors_X = sort(Eigenvectors(:,index_X));
        
    end
    Sigma_EPD(ii) = sigma;
    Eigenvalues_EPD(:,ii) = Eigenvalues_X;
    X_EPD(:,ii) = X;
    Ksd_pi_EPD(:,ii) = Ksd_pi;
    Eigenvectors_EPD(:,:,ii) = Eigenvectors_X;
end
display('3PhaseDelays_1CouplCoeff_Parameters')
[Sigma_min, index] = min(Sigma_EPD);
Sigma_min
Eigenvalues_min = Eigenvalues_EPD(:,index)
Ksd_pi_min = Ksd_pi_EPD(:,index)
Eigenvectors_min = Eigenvectors_EPD(:,:,ii);
X_min = X_EPD(:,index)

%%
Freq_Center = 193.54;
CouplCoeff_k1 = X_min(1);
Alpha = X_min(2);
Alpha_2 = X_min(3);
Radius = 10e-6;


% Frequency Sweep
Freqmax = Freq_Center + 0.04;
Freqmin = Freq_Center - 0.04;
Freqsteps = 50000;

Freq = linspace(Freqmin, Freqmax, Freqsteps);
TransferMatrices_Vec = zeros(6,6,Freqsteps);
Eigenvectors = zeros(6,Freqsteps);

for ii = 1:Freqsteps    
    TransferMatrices_Vec(:,:,ii) =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2);
end
[EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
for ii = 1:Freqsteps
    Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
end

% plot dispersion
    
Eigenvalues_plt = Eigenvalues'/pi;
Eigenvalues_plt_real = real(Eigenvalues_plt);
Eigenvalues_plt_imag = imag(Eigenvalues_plt);
Freq_plt = (ones(6,1)*Freq./Freq_Center)';

figure(1)
subplot(2,2,1)
hold on
plot(Eigenvalues_plt_real,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
% line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];  
for k = 1:6
    plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
end
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Re[kd/\pi]')
ylabel('Frequency (THz)')
grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

subplot(2,2,2)
hold on
plot(Eigenvalues_plt_imag,Freq_plt,'.')
% plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
for k = 1:6
    plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
end
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
xlabel('Im[kd/\pi]')
ylabel('Frequency (THz)')
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
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

% subplot(2,2,4)
% hold on
% plot(X_EPD(1,:)/Freq_Center,Sigma_EPD,'.')
% xlabel('Frequency(THz)')
% ylabel('\sigma: Coalescence Factor')
% grid on
% set(gca,'FontSize',24,'FontName', 'Times New Roman');

%% If the results are good, save:

writematrix(Sigma_EPD,'ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Sigma_EPD_Freq_pA_-pi.csv');
writematrix(Eigenvalues_EPD,'ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Eigenvalues_EPD_Freq_pA_-pi.csv');
writematrix(X_EPD,'ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Parameters_EPD_Freq_pA_-pi.csv');






%% Load, Plot and Analyze

%
clear all; close all; clc

Sigma_EPD = load('ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Sigma_EPD_Freq_pA_-pi.csv');
Eigenvalues_EPD = load('ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Eigenvalues_EPD_Freq_pA_-pi.csv');
X_EPD = load('ListsOfEPDs_WithPotentialSIPs/4-8-2021_1e4_3P2k_1_Parameters_EPD_Freq_pA_-pi.csv');
%%
close all;
% intervale = [1:50];
j=1;
interval = [];
for ii = 1:50
    inter = [j:j+200];
    j = j+200;
    interval = [interval inter'];
end
jj = 0;
for jj = 1:50
    [Sigma_min, index] = min(Sigma_EPD(interval(:,jj)));
    Sigma_min;
    Eigenvalues_min = Eigenvalues_EPD(:,index);
    X_min = X_EPD(:,index);
    
    Freq_Center = X_min(1);
    CouplCoeff_k1 = X_min(2);
    pA = X_min(3);
    pB = X_min(4);
    pC = X_min(5);
    
    
    % Frequency Sweep
    Freqmax = Freq_Center + 0.1;
    Freqmin = Freq_Center - 0.1;
    Freqsteps = 5000;
    
    Freq_Center = linspace(Freqmin, Freqmax, Freqsteps);
    TransferMatrices_Vec = zeros(6,6,Freqsteps);
    Eigenvectors = zeros(6,Freqsteps);
    
    for ii = 1:Freqsteps
        TransferMatrices_Vec(:,:,ii) =  get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq_Center, CouplCoeff_k1, Radius, Alpha, Alpha_2);
    end
    [EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
    for ii = 1:Freqsteps
        Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    end
    
    % plot dispersion
    
    Eigenvalues_plt = Eigenvalues'/pi;
    Eigenvalues_plt_real = real(Eigenvalues_plt);
    Eigenvalues_plt_imag = imag(Eigenvalues_plt);
    Freq_plt = (ones(6,1)*Freq_Center/Freq_Center)';
    
    figure(jj)
    subplot(2,2,1)
    hold on
    plot(Eigenvalues_plt_real,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    % line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
    line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];
    for k = 1:6
        plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Re[kd/\pi]')
    ylabel('Frequency (THz)')
    grid on
    % set(gca,'FontSize',24,'FontName', 'Times New Roman');
    
    subplot(2,2,2)
    hold on
    plot(Eigenvalues_plt_imag,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    for k = 1:6
        plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Im[kd/\pi]')
    ylabel('Frequency (THz)')
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
end

%%
close all; clc

jj=0;
kk = 1;
for ii=1:length(Sigma_EPD)
    if X_EPD(2,ii) < 1.035 && X_EPD(2,ii) > 1 && X_EPD(3,ii) < 1.035 && X_EPD(3,ii) > 1 
        JJ(kk) = ii;
        kk = kk+1;
    end
end

JJ = [1 2 5 14 37];
format long g
for jj = 1:1%length(JJ)
%     jj = JJ(jj)
    jj = 1;
    Sigma_EPD(jj)
    Eigenvalues_min = Eigenvalues_EPD(:,jj)
    Ksd_min = -log(Eigenvalues_min)./(1j)/pi
    X_min = X_EPD(:,jj)
    
    Freq_Center = 193.54;
    CouplCoeff_k1 = X_min(1);
    Alpha = X_min(2);
    Alpha_2 = X_min(3);
    Radius = 10e-6;

    
    
    % Frequency Sweep
    Freqmax = Freq_Center + 0.1;
    Freqmin = Freq_Center - 0.1;
    Freqsteps = 5000;
    
    Freq = linspace(Freqmin, Freqmax, Freqsteps);
    TransferMatrices_Vec = zeros(6,6,Freqsteps);
    Eigenvectors = zeros(6,Freqsteps);
    
    for ii = 1:Freqsteps    
    TransferMatrices_Vec(:,:,ii) = get_UnitCell_3_PhaseDelay_1CouplCoeff_Parameters (Freq(ii), CouplCoeff_k1, Radius, Alpha, Alpha_2);
    end
    [EigenVectors,EigenValues] = eigenshuffle(TransferMatrices_Vec);
    for ii = 1:Freqsteps
        Eigenvalues(:,ii) = -log(EigenValues(:,ii))./(1j);
    end
    
    % plot dispersion
    
    Eigenvalues_plt = Eigenvalues'/pi;
    Eigenvalues_plt_real = real(Eigenvalues_plt);
    Eigenvalues_plt_imag = imag(Eigenvalues_plt);
    Freq_plt = (ones(6,1)*Freq./Freq_Center)';
    
    figure(jj)
    subplot(2,2,1)
    hold on
    plot(Eigenvalues_plt_real,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    % line_color = ['b' 'g' 'y' 'c' 'm' 'r'];
    line_color = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]	; 	[0.9290, 0.6940, 0.1250] ; [0.4940, 0.1840, 0.5560] ; [0.4660, 0.6740, 0.1880]; [0.3010, 0.7450, 0.9330]];
    for k = 1:6
        plot(Eigenvalues_plt_real(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Re[kd/\pi]')
    ylabel('Frequency (THz)')
    grid on
    % set(gca,'FontSize',24,'FontName', 'Times New Roman');
    
    subplot(2,2,2)
    hold on
    plot(Eigenvalues_plt_imag,Freq_plt,'.')
    % plot(linspace(-1,1,1e3),Freq_Center.*ones(1e3),'r')
    for k = 1:6
        plot(Eigenvalues_plt_imag(:,k),Freq_plt(:,k),'.','Color',line_color(k,:))
    end
    axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
    xlabel('Im[kd/\pi]')
    ylabel('Frequency (THz)')
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
end