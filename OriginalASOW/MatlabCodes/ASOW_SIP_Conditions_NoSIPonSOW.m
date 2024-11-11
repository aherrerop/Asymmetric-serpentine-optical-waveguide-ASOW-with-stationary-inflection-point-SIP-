%% No SIP on SOW
clear all; close all; clc
k = linspace(0,1,1e2);
t = 1-k.^2;
ksd = acos((1-k.^2)./(3.*k.^2));
A = -t.^4./k.^4 +2.*(t.^6-2.*t.^4+t.^2)./k.^6;
Ksd = (1/2).*acos((1/6).*(-A-9)); 

figure(1)
hold on
plot(k,ksd,'b')
plot(k,Ksd,'r')
legend('k_sd for \beta = 0','k_sd \forall \beta')
xlabel('\kappa')
ylabel('k_sd')
grid on

