w_line=3;
% close all
BC=[255   , 135    ,0]/255;
GC='g';%[153   , 255    ,0]/255;

% kdhigh1=kdhigh1-real(kdhigh1)+wrapTo2Pi_RM(real(kdhigh1));
% kdhigh2=kdhigh2-real(kdhigh2)+wrapTo2Pi_RM(real(kdhigh2));
% kdhigh3=kdhigh3-real(kdhigh3)+wrapTo2Pi_RM(real(kdhigh3));
% kdhigh4=kdhigh4-real(kdhigh4)+wrapTo2Pi_RM(real(kdhigh4));
% kdhigh5=kdhigh5-real(kdhigh5)+wrapTo2Pi_RM(real(kdhigh5));
% kdhigh6=kdhigh6-real(kdhigh6)+wrapTo2Pi_RM(real(kdhigh6));


figure(1);
hold on;

plot(real(kdhigh2)/pi,Freq_high_plt,'-','linewidth',w_line,'color','g');
plot(real(kdhigh3)/pi,Freq_high_plt,'-','linewidth',w_line,'color','r');
plot(real(kdhigh4)/pi,Freq_high_plt,'--','linewidth',w_line,'color',BC);
plot(real(kdhigh5)/pi,Freq_high_plt,'--','linewidth',w_line,'color','b');
plot(real(kdhigh6)/pi,Freq_high_plt,'--','linewidth',w_line,'color','k');
plot(real(kdhigh1)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');


plot(real(kdlow1)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
plot(real(kdlow2)/pi,Freq_low_plt,'-','linewidth',w_line,'color','b');
plot(real(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
plot(real(kdlow4)/pi,Freq_low_plt,'--','linewidth',w_line,'color','g');
plot(real(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(real(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');


set(gca,'FontSize',24,'FontName', 'Times New Roman');
xlabel('Re($kd / \pi$)','FontSize', 24,'Interpreter','latex');
ylabel('$\omega / \omega_d$','FontSize', 24,'Interpreter','latex');
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
grid on;
% pbaspect([1.3 1 1]);
% %

figure(2)
hold on;
plot(imag(kdhigh3)/pi,Freq_high_plt,'-','linewidth',w_line,'color','b');
plot(imag(kdhigh5)/pi,Freq_high_plt,'-','linewidth',w_line,'color','r');
plot(imag(kdhigh4)/pi,Freq_high_plt,'--','linewidth',w_line,'color','g');
plot(imag(kdhigh2)/pi,Freq_high_plt,'--','linewidth',w_line,'color',BC);


plot(imag(kdlow4)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
plot(imag(kdlow1)/pi,Freq_low_plt,'-','linewidth',w_line,'color','g');
plot(imag(kdlow2)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
plot(imag(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','b');
plot(imag(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(imag(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(imag(kdhigh6)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');
plot(imag(kdhigh1)/pi,Freq_high_plt,'-','linewidth',w_line,'color','k');

xlabel('Im($kd / \pi$)','FontSize', 24,'Interpreter','latex');
ylabel('$\omega / \omega_d$','FontSize', 24,'Interpreter','latex');
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
grid on;
% pbaspect([1.3 1 1]);


% % 
figure(3);
hold on;
plot(real(kdhigh3)/pi,imag(kdhigh3)/pi,'-','linewidth',w_line,'color','b');
plot(real(kdhigh4)/pi,imag(kdhigh4)/pi,'-','linewidth',w_line,'color','g');

plot(real(kdhigh2)/pi,imag(kdhigh2)/pi,'-','linewidth',w_line,'color',BC);
plot(real(kdhigh5)/pi,imag(kdhigh5)/pi,'-','linewidth',w_line,'color','r');

plot(real(kdlow3)/pi,imag(kdlow3)/pi,'-','linewidth',w_line,'color','b');
plot(real(kdlow1)/pi,imag(kdlow1)/pi,'-','linewidth',w_line,'color','g');
plot(real(kdlow2)/pi,imag(kdlow2)/pi,'-','linewidth',w_line,'color','r');
plot(real(kdlow4)/pi,imag(kdlow4)/pi,'-','linewidth',w_line,'color',BC);

plot(real(kdhigh6)/pi,imag(kdhigh6)/pi,'-','linewidth',w_line,'color','k');
plot(real(kdhigh1)/pi,imag(kdhigh1)/pi,'-','linewidth',w_line,'color','k');
plot(real(kdlow5)/pi,imag(kdlow5)/pi,'-','linewidth',w_line,'color','k');
plot(real(kdlow6)/pi,imag(kdlow6)/pi,'-','linewidth',w_line,'color','k');

grid on;
xlabel('Re($kd / \pi$)','FontSize', 24,'Interpreter','latex');
ylabel('Im($kd / \pi$)','FontSize', 24,'Interpreter','latex');
grid on;
% pbaspect([1.3 1 1]);
axis([-0.6,0.6, -0.5,0.5]);

figure(4);
hold on;
% plot(Freq_low_plt,Sigma,'b.');
plot(Freq_high_plt,Sigma,'b.');
hold off;
ylabel('Coalescence parameter, $\sigma$','FontSize', 24,'Interpreter','latex');
xlabel('$\omega / \omega_d$','FontSize', 24,'Interpreter','latex');
set(gca,'FontSize',24,'FontName', 'Times New Roman');
xlim([Freqmin/Freq_Center, Freqmax/Freq_Center]);
grid on;
% 
