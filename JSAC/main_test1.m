clear all
clc   

addpath(pwd);
cd FPS-AltMin/;
addpath(genpath(pwd));
cd ..;

SNR_dB=[-20:2:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=100; 
Nt=64;Nr=16;Ns = 6;Nc = 30;
Nt_RF=6;N_RF=6;
V_RF=eye(Nt,N_RF);
W_RF=eye(Nr,N_RF);
V_RF1=eye(Nt,N_RF);
W_RF1=eye(Nr,N_RF);
L=5; % number of rays
M=Nt/Nt_RF; % number of antennas connected to one RF chains
fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;
c = 1/sqrt(Nc)*exp(1i*[0:2*pi/Nc:2*pi-2*pi/Nc])';
C = kron(eye(N_RF),c);
for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    temp1=0;temp2=0;temp3=0;temp4=0;temp5=0;temp6=0;temp7=0;temp8=0;
    P=sqrt(Nr*Nt/L);
    sigma2=P/SNR;
%     sigma2=SNR/16;
    for i=1:N_iter
        [H,power_matrix,A_BS,~,Fopt,Wopt]=mmWave_channel(Nr,Nt,Ns,L,lamada);
%         H=randn(Nr,Nt)+1i*randn(Nr,Nt);
        %%%%%% conventional method %%%%%%%%%%
%         F1=SVD_precoding(Nt_RF,H);
%         temp1=temp1+log2(det(eye(Nr)+(SNR/16)*H*F1*F1'*H'));
%         %%%%% proposed method %%%%%%%%%%%%%%%
%         F2=propose_precoding(Nt_RF,Nt,Nr,M,H,SNR);
%         temp2=temp2+log2(det(eye(Nr)+(SNR/16)*H*F2*F2'*H'));
%         %%%%%% analog method %%%%%%%%%%%%%%%
% %         F3=full_analog(Nt_RF,Nt,Nr,M,H,SNR);
% %         temp3=temp3+log2(det(eye(Nr)+(SNR/16)*H*F3*F3'*H'));
% 
%         %%%%%% spatially_sparse_precoding %%%%%%%%%%%%%%%%%%%%%%%%%
%         [F_RF,F_BB]=spatially_sparse_precoding(Nt_RF,H,A_BS);
%         F4=F_RF*F_BB;
%         temp4=temp4+log2(det(eye(Nr)+(SNR/16)*H*F4*F4'*H'));
        %%%% SVD %%%%%%%%%%%%%%%%%%%%
%         F5=hybrid_precoding(Nt_RF,Nt,Nr,M,H,SNR);
%         temp5=temp5+log2(det(eye(Nr)+(SNR/16)*H*F5*F5'*H'));
        %%%%%5月24日例会讲的算法%%%%%
        %[1]Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays
        [V_t,V_RF,V_D] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2);
        [W_t] = hybrid_beamforming_wt_w(H,N_RF,W_RF,V_t,Nr,P,sigma2);
        temp6=temp6+log2(det(eye(Nr)+(1/sigma2)*W_t*((W_t'*W_t)^(-1))*W_t'*H*V_t*V_t'*H'));
        %%%%%5月24日例会讲的算法%%%%%
        %[1]Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays
%         [V_t1,V_RF1,V_D1] = hybrid_beamforming_vt2_w(H,N_RF,V_RF1,P,Nt,sigma2);
%         [W_t1] = hybrid_beamforming_wt2_w(H,N_RF,W_RF1,V_t1,Nr,P,sigma2);
%         temp7=temp7+log2(det(eye(Nr)+(1/sigma2)*W_t1*((W_t1'*W_t1)^(-1))*W_t1'*H*V_t1*V_t1'*H'));
        
        %%%%% FPS-AltMin算法 %%%%%
        %Hybrid Precoding in Millimeter Wave Systems:How Many Phase Shifters Are Needed?
%         [ H1,At,Ar,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns);
        [ FRF, FBB ] = my_AltMin_new( Fopt, C);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [ WRF, WBB ] = my_AltMin_new( Wopt, C);
%         
%         [V_t1,V_RF1,V_D1] = hybrid_beamforming_vt2_w(H,N_RF,V_RF1,P,Nt,sigma2);
%         [W_t1] = hybrid_beamforming_wt2_w(H,N_RF,W_RF1,V_t1,Nr,P,sigma2);
        temp8=temp8+log2(det(eye(Ns) +  (P/(sigma2*Ns))* pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
%                     ()log2(det(eye(Ns) +  * pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
    end
    C1(i_snr)= real(temp1/N_iter);
    C2(i_snr)= real(temp2/N_iter);     
    C3(i_snr)= real(temp3/N_iter);
    C4(i_snr)= real(temp4/N_iter);
    C5(i_snr)= real(temp5/N_iter);
    C6(i_snr)= real(temp6/N_iter);
    C7(i_snr)= real(temp7/N_iter);
    C8(i_snr)= real(temp8/N_iter);
end
plot(SNR_dB,C1,'r-o','Linewidth',1.5);
hold on
plot(SNR_dB,C2,'b-s','Linewidth',1.5);
hold on
% plot(SNR_dB,C3,'g-^','Linewidth',1.5);
% hold on
plot(SNR_dB,C4,'c-^','Linewidth',1.5);
hold on
% plot(SNR_dB,C5,'m-^','Linewidth',1.5);
% hold on
plot(SNR_dB,C6,'k-^','Linewidth',1.5);
% hold on
% plot(SNR_dB,C7,'m-^','Linewidth',1.5);
hold on
plot(SNR_dB,C8,'m-^','Linewidth',1.5);
xlabel('SNR (dB)')
ylabel('Capacity')
grid on 
leg1='conventional method(full-connected)';%设置图例
leg2='SIC-based hybrid precoding';
leg3='full-analog method';
leg4='spatially-sparse-precoding';
leg5='SVD hybrid precoding';
leg6='hybrid precoding(metting)';
leg7='hybrid precoding2(metting)';
leg8='111';
legend(leg6,leg8);
title(['iter = ',num2str(N_iter), ' , Nt = ',num2str(Nt),' , Nr = ',num2str(Nr),' , RF chains = ',num2str(N_RF),]);