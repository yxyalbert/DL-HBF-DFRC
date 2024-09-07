close;
clear all£»
SNR=0:5:30;
%mse2=cell2mat(struct2cell(load('mse2.mat')));

DBF=[23.75,23.75,23.75,23.75,23.75,23.75,23.75];
BEST=[21.450,21.450,21.450,21.450,21.450,21.450,21.450];
CNN_HBF=[20.3,20.75,21.1,21.28,21.36,21.4,21.43];
Two_Stage=[20,20.26,20.57,21.0,21.25,21.33,21.38];
PE_MO=[18.7,19.0,19.5,20.2,20.7,20.9,21.0];
OMP=[12.3,13.8,15.4,17.6,19.2,20.0,20.4];
 figure(1)
 fs = 12;
 linewidth = 1.5;

%plot(NRF,eefully, '-or', 'LineWidth', linewidth),hold on;

plot(SNR,DBF, '-or', 'LineWidth', linewidth),hold on;
plot(SNR,BEST, '-b', 'LineWidth', linewidth),hold on;
plot(SNR,CNN_HBF, '-og', 'LineWidth', linewidth),hold on;
plot(SNR,Two_Stage, '-k', 'LineWidth', linewidth),hold on;
plot(SNR,PE_MO, '-m', 'LineWidth', linewidth),hold on;
plot(SNR,OMP, '-s', 'LineWidth', linewidth),hold on;
% 

grid on
hold on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('N_{RF}')
ylabel('Sum-Rate (bits/s/Hz/W)')
legend('DBF', 'BEST','CNN-HBF' ,'Two-Stage','PE-MO','OMP');