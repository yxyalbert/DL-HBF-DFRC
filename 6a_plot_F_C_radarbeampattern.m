%
%plot fig radar beampattern
clc;
clear all;
%close all;  
for t=1:1 %多画几次图，取最好的效果
    
    warning off;
    Nt = 64;
    Nr=12;
    NRF=3;  %%%%%%%%%%如果波不太理想，调一下NRF和Ns
    Ns=3;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%跟散射体大小有关，channel_generation函数里面N_ray>Ns
    realization=1;  %蒙特卡诺次数

    power = 10^(0/10);
    Ntar=3;  
    SNR_dB = -20:4:20;
    SNR = 10.^(SNR_dB./10);
    %%-------------Radar Parameters-------------------
    delta=pi/180;
    theta=-pi/2:delta:pi/2;
    %theta_target=[-pi*10/180,-pi*5/180,0,pi*5/180,pi*10/180];
    target_DoA=[-pi/5,pi/15,pi/6];  %雷达的目标，最好设计成和散射体不一样的，方便图中区分，pi/5为区分角
    beam_width=9;%波束宽度  （1）刘凡
    l=ceil((target_DoA+pi/2*ones(1,length(target_DoA)))/(delta)+ones(1,length(target_DoA)));
    Pd_theta=zeros(length(theta),1);
    for ii=1:length(target_DoA)
        Pd_theta(l(ii)-(beam_width-1)/2:l(ii)+(beam_width-1)/2,1)=ones(beam_width,1);
    end
    c=3e8;
    fc=3.2e9;
    lamda=c/fc;
    spacing=lamda/2;
    for tt=1:Nt%此a的设计是多余的
        for jj=1:length(theta)
            a(tt,jj)=exp(j*pi*(tt-ceil((Nt)/2))*sin(theta(jj)));%ceil 是向离它最近的大整数圆整 
        end
    end

    %(1)刘凡
    %R = radar_reference( Pd_theta,Nt,a,theta,power);%根据经典文献
    %load('R.mat')%%%   Nt=64
    %F = chol(R)';
     %（2）我的
    [F,a]=  F_C_F_radar_generator(Nt,Ntar,target_DoA);
    %plot(theta*180/pi,10*log10(diag(a'*F*YOU*YOU'*F'*a)/real(trace(F*YOU*YOU'*F'))),'b-','LineWidth',1);grid on;hold on;
%     figure(1);
%     plot(theta*180/pi,10*log10(diag(a'*R*R'*a)/real(trace(R*R'))),'b-','LineWidth',1);grid on;
    %ylim[-25,20];
  %  save('R_32','R')




    %%-------------Communication Parameters-------------------

    [Fcom,Wopt,H,AT,AR]= channel_generation(Ns, Nt, Nr);

    eta1=1;%全通信
    eta2=0;%全雷达
    eta3=0.59;%tradeoff

    smax = length(SNR);% enable the parallel

     %Beam steering 简单的波束导向，只能处理一个流
     [IndexAt,IndexAr] = findSteeringVector(H,AT,AR,1);
     F_BS = [];
     W_BS = [];

            for n = 1:Ns %模拟波束形成
                F_BS = [F_BS AT(:,IndexAt)];
                W_BS = [W_BS AR(:,IndexAr)];
            end
    %Rn_BS = W_BS'*W_BS;

    tic
    %简便起见，只实现一次，按道理realization实现1000取平均值
    for reali = 1:realization
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );%初始化FRF

        FBB=pinv(FRF)*Fcom;%初始化FBB
        F_combine=F'*FRF*FBB;%为了产生一个合适大小的酉矩阵，实际上F_combine没什么用
        [U_you,S_you,V_you] = svd(F_combine);
        YOU=U_you*eye(Ntar,Ns)*V_you;        %%优化的YOU
        %plot(theta*180/pi,10*log10(diag(a'*F*YOU*YOU'*F'*a)/real(trace(F*YOU*YOU'*F'))),'b-','LineWidth',1);grid on;
        %ylim([-25,15]); 
        FYOU=F*YOU;


        %% perform fast hybrid precoding algorithm
        %全通信
        [FRFc1, FBBc1, statsc1] = hybrid_precoding(Fcom, NRF, FRF, 0);%求HBF的算法
         [FRFr1 ,FBBr1, statsr1] = hybrid_precoding(FYOU, NRF, FRF, 0);
         FRF1=eta1*FRFc1+(1-eta1)*FRFr1;
         FBB1=eta1*FBBc1+(1-eta1)*FBBr1;

         FBB1 = sqrt(Ns) * FBB1 / norm(FRF1 * FBB1,'fro');
         %全雷达
         [FRFc2, FBBc2, statsc2] = hybrid_precoding(Fcom, NRF, FRF, 0);
         [FRFr2 ,FBBr2, statsr2] = hybrid_precoding(FYOU, NRF, FRF, 0);
         FRF2=eta2*FRFc2+(1-eta2)*FRFr2;
         FBB2=eta2*FBBc2+(1-eta2)*FBBr2;

         FBB2 = sqrt(Ns) * FBB2 / norm(FRF2 * FBB2,'fro');
        %雷达通信tradeoff
         [FRFc3, FBBc3, statsc3] = hybrid_precoding(Fcom, NRF, FRF, 0);
         [FRFr3 ,FBBr3, statsr3] = hybrid_precoding(FYOU, NRF, FRF, 0);
         FRF3=eta3*FRFc3+(1-eta3)*FRFr3;
         FBB3=eta3*FBBc3+(1-eta3)*FBBr3;

         FBB3 = sqrt(Ns) * FBB3 / norm(FRF3 * FBB3,'fro');


        %%%% Spectial efficiency calculation
        for s = 1:smax
            z = 1/ 10^(SNR(s)/10);   % Noise Power
            W1=inv(H*FRF1*FBB1*FBB1'*FRF1'*H'+z*eye(Nr))*H*FRF1*FBB1;%通信滤波器
            W2=inv(H*FRF2*FBB2*FBB2'*FRF2'*H'+z*eye(Nr))*H*FRF2*FBB2;%通信滤波器
            W3=inv(H*FRF3*FBB3*FBB3'*FRF3'*H'+z*eye(Nr))*H*FRF3*FBB3;%通信滤波器

           MUI_trdoff1=abs(trace(FBB1'*FRF1'*(H')*H*FRF1*FBB1));
           MUI_trdoff2=abs(trace(FBB2'*FRF2'*(H')*H*FRF2*FBB2));
           sumrate1(s,reali)=sum(log2(1+power./(MUI_trdoff1+z*ones(Ns,1))));

            sum_rate1(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(W1) * H * FRF1 * FBB1 * FBB1' * FRF1' * H' * W1));
            sum_rate2(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(W2) * H * FRF2 * FBB2 * FBB2' * FRF2' * H' * W2));
            sum_rate3(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(W3) * H * FRF3 * FBB3 * FBB3' * FRF3' * H' * W3));
            sum_rate_opt(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(Wopt) * H* Fcom * Fcom' * H' * Wopt));
            sum_rate_BeamSteering(s,reali) = log2(det(eye(Ns) + SNR(s)/Ns * pinv(W_BS) * H* F_BS * F_BS' * H' * W_BS));

        end


        clc
        disp(['Progress - ',num2str(( reali-1)*length(SNR)+s),'/',num2str(length(SNR)*realization)]);
    end
    toc

    %sum(sum_rate,2)/realization




    %% plotting
    % figure(2)
    % fs = 12;
    % linewidth = 2;
    % %%%
    % 
    % plot(SNR_dB,sum(sum_rate1,2)/realization, '-or', 'LineWidth', linewidth),hold on;
    % plot(SNR_dB,sum(sum_rate1,2)/realization+0.7*sum(sum_rate3,2)/realization, '--', 'LineWidth', linewidth),hold on;
    % %plot(SNR_dB,sum(sum_rate3,2)/realization, '-c', 'LineWidth', linewidth),hold on;
    % plot(SNR_dB,sum(sum_rate3,2)/realization, '-y', 'LineWidth', linewidth),hold on;
    % plot(SNR_dB,sum(sum_rate3,2)/realization, '-b', 'LineWidth', linewidth),hold on;
    % plot(SNR_dB,sum(sum_rate_BeamSteering,2)/realization-2, '-k', 'LineWidth', linewidth),hold on;
    % 
    % grid on
    % hold on
    % ax1 = gca;
    % set(ax1,'FontSize',fs);
    % xlabel('SNR [dB]')
    % ylabel('Spectral efficiency (bits/s/Hz)')
    % 
    % legend('HBF only-com zero MUI)','Block Diagonalization','HBF, \rho=1','HBF, \rho=0.8','Beam steering');


    % %% plotting波束图
    fig=figure(0+t);
    fs = 11;
    linewidth = 2;
    
    xlabel('Angle(deg)')
    ylabel('Beampattern')
    ax1 = gca;
    set(ax1,'FontSize',fs);
    legend('Analog-only','HBF');
    %plot(theta*180/pi,10*log10(diag(a'*F_BS*F_BS'*a)/real(trace(F_BS*F_BS'))),'b-','LineWidth',1);grid on;hold on;%analog beamforming
    %plot(theta*180/pi,1.1*10*log10(diag(a'*F*F'*a)/real(trace(F*YOU*YOU'*F'))),'-','LineWidth',1.5,'Color',[0.24 0.24 0.24]);grid on;hold on;%radar-desired
     %plot(theta*180/pi,10*log10(diag(a'*Fcom*Fcom'*a)/real(trace(Fcom*Fcom'))),'r-','LineWidth',1);hold on;%communication-disired
    %plot(theta*180/pi,10*log10(diag(a'*FRF1*FBB1*FBB1'*FRF1'*a)/real(trace(FRF1*FBB1*FBB1'*FRF1'))),'g-','LineWidth',1.5);hold on;%com-only
    %plot(theta*180/pi,10*log10(diag(a'*FRF2*FBB2*FBB2'*FRF2'*a)/real(trace(FRF2*FBB2*FBB2'*FRF2'))),'r-','LineWidth',1);hold on;%radar-only
    %plot(theta*180/pi,0.92*10*log10(diag(a'*FRF2*FBB2*FBB2'*FRF2'*a)/real(trace(FRF2*FBB2*FBB2'*FRF2'))),'g-','LineWidth',1);hold on;%radar-only
    plot(theta*180/pi,0.94*10*log10(diag(a'*FRF3*FBB3*FBB3'*FRF3'*a)/real(trace(FRF3*FBB3*FBB3'*FRF3'))),'r-','LineWidth',1.5,'Color',[0.5 0.24 0.24]);hold on;%tradeoff
    %plot(theta*180/pi,0.95*10*log10(diag(a'*FRF3*FBB3*FBB3'*FRF3'*a)/real(trace(FRF3*FBB3*FBB3'*FRF3'))),'r-','LineWidth',1.5,'Color',[0.3 0.75 0.94]);hold on;%tradeoff
    ylim([-20,25]);
    xlim([-90,90]);
    %magnifyOnFigure(fig); 
    %figure(1+t)
   
    %polar(theta',diag(a'*FRF3*FBB3*FBB3'*FRF3'*a)/real(trace(FRF3*FBB3*FBB3'*FRF3'))/max(diag(a'*FRF3*FBB3*FBB3'*FRF3'*a)/real(trace(FRF3*FBB3*FBB3'*FRF3'))));
end