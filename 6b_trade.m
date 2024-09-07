close;
clear all£»
MSE_10=[0.73,0.65,0.52,0.46,0.32,0.21];
MSE_20=[0.68,0.61,0.54,0.42,0.30,0.19];
MSE_30=[0.60,0.57,0.50,0.36,0.27,0.16];
%mse2=cell2mat(struct2cell(load('mse2.mat')));

sum_rate_10=[21.1,20.9,20.2,19.8,18.6,17.4];
sum_rate_20=[21.31,21.0,20.7,20.0,19.1,18.0];
sum_rate_30=[21.43,21.3,21.0,20.23,19.4,18.3];
 figure(1)
 fs = 12;
 linewidth = 1.5;

%plot(NRF,eefully, '-or', 'LineWidth', linewidth),hold on;

plot(MSE_10,sum_rate_10, '-or', 'LineWidth', linewidth),hold on;
plot(MSE_20,sum_rate_20, '-ob', 'LineWidth', linewidth),hold on;
plot(MSE_30,sum_rate_30, '-ok', 'LineWidth', linewidth),hold on;
% 

grid on
hold on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('MSE_R')
ylabel('Sum-Rate (bits/s/Hz/W)')
legend('DBF', 'BEST','CNN-HBF' ,'Low-resolution','PE-MO','OMP');