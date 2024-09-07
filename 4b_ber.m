 close;
clear all£»
SNR=[-20,-15,-10,-5,0];
%mse2=cell2mat(struct2cell(load('mse2.mat')));



DBF=[0.1,0.023,0.0017,0.00004,0.00000009];
BEST=[0.13,0.03,0.0032,0.00009,0.0000003];
CNN_HBF=[0.135,0.04,0.0042,0.00012,0.00000044];
PE_MO=[0.135,0.05,0.007,0.00024,0.0000019];
Two_stage=[0.14,0.055,0.012,0.0005,0.000009];
OMP=[0.16,0.1,0.06,0.02,0.004];
 figure(1)
 fs = 12;
 linewidth = 1.5;

%plot(NRF,eefully, '-or', 'LineWidth', linewidth),hold on;

semilogy(SNR,DBF, '-or', 'LineWidth', linewidth),hold on;
semilogy(SNR,BEST, '-b', 'LineWidth', linewidth),hold on;
semilogy(SNR,CNN_HBF, '-og', 'LineWidth', linewidth),hold on;
semilogy(SNR,PE_MO, '-k', 'LineWidth', linewidth),hold on;
semilogy(SNR,Two_stage, '-m', 'LineWidth', linewidth),hold on;
semilogy(SNR,OMP, '-s', 'LineWidth', linewidth),hold on;
% 

grid on
hold on
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel('N_{RF}')
ylabel('Sum-Rate (bits/s/Hz/W)')
legend('DBF', 'BEST','CNN-HBF' ,'Two-stage','PE-MO','OMP');