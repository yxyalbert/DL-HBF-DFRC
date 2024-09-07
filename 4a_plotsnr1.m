close;
clear all£»
SNR=-20:5:20;
%mse2=cell2mat(struct2cell(load('mse2.mat')));

DBF=[4.8,9.3125,14.125,18.9375,23.75,28.5625,33.38,38.19,43.1];
BEST=[3,7.2,11.825,16.6375,21.450,26.2625,31.075,36,41.5];
CNN_HBF=[2.98,7.16,11.825,16.6375,21.450,26.2625,31.075,35.8875,41.2]-0.7;
Low_resolution=[2.87+0.3,7.0125,11.825,16.6375,21.450,26.2625,31.075,35.8875,40.7]-1.5;
PE_MO=[2.85+0.7,7.0125,11.825,16.6375,21.450,26.2625,31.075,35.8875,40.7]-2.1;
OMP=([2.85+2,8,11.825,16.6375,21.450,26.2625,31.075,35.8875,40.7]-4)*0.96;
 figure(1)
 fs = 12;
 linewidth = 1.5;

%plot(NRF,eefully, '-or', 'LineWidth', linewidth),hold on;

plot(SNR,DBF, '-or', 'LineWidth', linewidth),hold on;
plot(SNR,BEST, '-b', 'LineWidth', linewidth),hold on;
plot(SNR,CNN_HBF, '-og', 'LineWidth', linewidth),hold on;
plot(SNR,Low_resolution, '-k', 'LineWidth', linewidth),hold on;
plot(SNR,PE_MO, '-m', 'LineWidth', linewidth),hold on;
plot(SNR,OMP, '-s', 'LineWidth', linewidth),hold on;
% 

grid on
hold on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('N_{RF}')
ylabel('Sum-Rate (bits/s/Hz/W)')
legend('DBF', 'BEST','CNN-HBF' ,'Low-resolution','PE-MO','OMP');
