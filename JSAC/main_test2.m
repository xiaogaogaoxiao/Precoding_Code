clear all
clc   

addpath(pwd);
cd FPS-AltMin/;
addpath(genpath(pwd));
cd ..;
%%%%%%%%%%%%
% test github
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyb_method = [1 2 5 6 8];           % 混合预编码算法：
%[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% 3:full-analog method   4:spatially_sparse_precoding
% 5:SVD hybrid precoding  6:hybrid precoding(metting)
% 7:hybrid precoding2(metting)   8:FPS-AltMin Algorithm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_dB=[-10:2:6];  SNR_linear=10.^(SNR_dB/10.);
N_iter=1; 
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
num_of_HybridMethod = length(hyb_method);
CC=zeros(length(SNR_linear),num_of_HybridMethod);
for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    temp1=0;temp2=0;temp3=0;temp4=0;temp5=0;temp6=0;temp7=0;temp8=0;
    temp=zeros(1,num_of_HybridMethod);
    P=sqrt(Nr*Nt/L);
    sigma2=P/SNR;
%     sigma2=SNR/16;
    for i=1:N_iter
        for index_type = 1:length(hyb_method)
            tmptype = hyb_method(index_type);
            [H,power_matrix,A_BS,~,Fopt,Wopt]=mmWave_channel(Nr,Nt,Ns,L,lamada);
%             H=randn(Nr,Nt)+1i*randn(Nr,Nt);
            switch (tmptype)
                case 1    %-- conventional method(full-connected) 
                    F1=SVD_precoding(Nt_RF,H);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F1*F1'*H'));
                case 2    %-- SIC-based hybrid precoding
                    F2=propose_precoding(Nt_RF,Nt,Nr,M,H,SNR);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F2*F2'*H'));
                case 3    %-- full-analog method
                    F3=full_analog(Nt_RF,Nt,Nr,M,H,SNR);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F3*F3'*H'));
                case 4    %-- spatially_sparse_precoding
                    [F_RF,F_BB]=spatially_sparse_precoding(Nt_RF,H,A_BS);
                    F4=F_RF*F_BB;
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F4*F4'*H'));
                case 5    %-- SVD hybrid precoding
                    F5=hybrid_precoding(Nt_RF,Nt,Nr,M,H,SNR);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F5*F5'*H'));
                case 6    %-- point-to-point MIMO 场景的算法（5月24日例会讲的算法）
                    %[1]Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays
                    [V_t,V_RF,V_D] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2);
                    [W_t] = hybrid_beamforming_wt_w(H,N_RF,W_RF,V_t,Nr,P,sigma2);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(1/sigma2)*W_t*((W_t'*W_t)^(-1))*W_t'*H*V_t*V_t'*H'));
                case 7 
                    [V_t1,V_RF1,V_D1] = hybrid_beamforming_vt2_w(H,N_RF,V_RF1,P,Nt,sigma2);
                    [W_t1] = hybrid_beamforming_wt2_w(H,N_RF,W_RF1,V_t1,Nr,P,sigma2);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(1/sigma2)*W_t1*((W_t1'*W_t1)^(-1))*W_t1'*H*V_t1*V_t1'*H'));
                case 8    %-- Fixed Phase Shifters(FPS-AltMin)算法 
                    %Hybrid Precoding in Millimeter Wave Systems:How Many Phase Shifters Are Needed?
                    [ FRF, FBB ] = my_AltMin_new( Fopt, C);
                    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
                    [ WRF, WBB ] = my_AltMin_new( Wopt, C);
                    temp(index_type)=temp(index_type)+log2(det(eye(Ns) +  (P/(sigma2*Ns))* pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
                otherwise
                    error('Unknown TxRx.Decoder.LDPC.Type.')
            end
        end 
    end
    CC(i_snr,:)= real(temp/N_iter);
end
marker_style = {'o-','o--','+-','+--','v-','v--','<-','<--','x-','x--','^-','^--','*-','*--','d-','d--','h-','h--','^-','hexagram--','diamond:','pentagram:','v:'};
marker_color = [...
            0.0000    0.4470    0.7410;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0.3010    0.7450    0.9330;...
            0.6350    0.0780    0.1840;...
            0.7500    0.7500    0.0000;...
            0.7500    0.0000    0.7500;...
            0.0000    0.5000    0.0000;...
            0.0000    0.0000    1.0000;...
            0.2000    0.5000    0.6600;...
            0.6500    0.4500    0.5500;...
            0.6000    0.2000    0.5000;...
            0.2000    0.8880    1.0000;...
            0.3400    0.4300    0.6600;...
            0.3660    0.2880    1.0000;...
            0.6400    0.6500    0.8800;...
        ];
MarkerSize =7;
LineWidth = 2.2;

figure(1)
for i=1:num_of_HybridMethod
    plot(SNR_dB(1:end),CC(:,i),marker_style{i}, 'color', marker_color(i,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
%    semilogy(SNR_dB(1:end),CC(:,i),marker_style{i}, 'color', marker_color(i,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
    hold on;grid on;
end

legend('1','2','5','hybrid precoding(metting)','FPS-AltMin Algorithm');
%[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% 3:full-analog method   4:spatially_sparse_precoding
% 5:SVD hybrid precoding  6:hybrid precoding(metting)
% 7:hybrid precoding2(metting)   8:FPS-AltMin Algorithm]
title(['iter = ',num2str(N_iter), ' , Nt = ',num2str(Nt),' , Nr = ',num2str(Nr),' , RF chains = ',num2str(N_RF),]);
FontSize =10;
xlabel('SNR (dB)')
ylabel('Spectral efficiency (bits/s/Hz)')
hold on