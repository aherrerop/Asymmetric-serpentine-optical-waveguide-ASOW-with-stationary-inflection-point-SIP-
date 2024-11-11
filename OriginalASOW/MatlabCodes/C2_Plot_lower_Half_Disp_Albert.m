w_line=3;
BC=[255   , 135    ,0]/255;
Wf=14;
% fig = figure('units','inch','position',[0,0,Wf,Wf/2]);
hold on;
figure(1)
% kdlow1=kdlow1-real(kdlow1)+wrapTo2Pi_RM(real(kdlow1));
% kdlow2=kdlow2-real(kdlow2)+wrapTo2Pi_RM(real(kdlow2));
% kdlow3=kdlow3-real(kdlow3)+wrapTo2Pi_RM(real(kdlow3));
% kdlow4=kdlow4-real(kdlow4)+wrapTo2Pi_RM(real(kdlow4));
% kdlow5=kdlow5-real(kdlow5)+wrapTo2Pi_RM(real(kdlow5));
% kdlow6=kdlow6-real(kdlow6)+wrapTo2Pi_RM(real(kdlow6));
hold on

plot(real(kdlow2)/pi,Freq_low_plt,'-','linewidth',w_line,'color','b');
plot(real(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','g');
plot(real(kdlow4)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
plot(real(kdlow1)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
plot(real(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(real(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
xlabel('Re(\itkdlow\rm/\pi)','FontSize', 24);
ylabel('\omega/\omega_d','FontSize', 24);
set(gca,'FontSize',24);
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);
grid on;
set(gca,'FontSize',24,'FontName', 'Times New Roman');
pbaspect([1.3 1 1]);


figure(2)
hold on;
plot(imag(kdlow2)/pi,Freq_low_plt,'-','linewidth',w_line,'color','b');
plot(imag(kdlow4)/pi,Freq_low_plt,'-','linewidth',w_line,'color',BC);
plot(imag(kdlow1)/pi,Freq_low_plt,'--','linewidth',w_line,'color','g');
plot(imag(kdlow3)/pi,Freq_low_plt,'--','linewidth',w_line,'color','r');
plot(imag(kdlow6)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');
plot(imag(kdlow5)/pi,Freq_low_plt,'-','linewidth',w_line,'color','k');

% plot(imag(kdlow3)/pi,ff/fo);
% plot(imag(kdlow4)/pi,ff/fo);
xlabel('Im(\itkdlow\rm/\pi)','FontSize', 24);
ylabel('\omega/\omega_d','FontSize', 24);
set(gca,'FontSize',24);
% ylim([fo-df fo+df]/fo);
grid on;
set(gca,'FontSize',24,'FontName', 'Times New Roman');
pbaspect([1.3 1 1]);
axis([-1,1, Freqmin/Freq_Center, Freqmax/Freq_Center]);

figure(3);
hold on;
plot(real(kdlow1)/pi,imag(kdlow1)/pi,'-','linewidth',w_line,'color',BC);
plot(real(kdlow2)/pi,imag(kdlow2)/pi,'-','linewidth',w_line,'color','b');
plot(real(kdlow3)/pi,imag(kdlow3)/pi,'-','linewidth',w_line,'color','r');
plot(real(kdlow4)/pi,imag(kdlow4)/pi,'-','linewidth',w_line,'color','g');
plot(real(kdlow5)/pi,imag(kdlow6)/pi,'-','linewidth',w_line,'color','k');
plot(real(kdlow6)/pi,imag(kdlow5)/pi,'-','linewidth',w_line,'color','k');

grid on;
xlabel('Real (\itkd\rm/\pi)','FontSize', 24);
ylabel('Imag (\itkd\rm/\pi)','FontSize', 24);
set(gca,'FontSize',24,'FontName', 'Times New Roman');
axis([-1,1, -1, 1]);



figure(4);
hold on;
plot(Freq_low_plt,Sigma,'b.');
ylabel('Coalescence parameter, \sigma','FontSize', 28);
xlabel('\omega/\omega_d','FontSize', 28);
set(gca,'FontSize',28);
grid on;
set(gca,'FontSize',28,'FontName', 'Times New Roman');

