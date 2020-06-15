clear all
clc   

% addpath(pwd);
% cd FPS-AltMin/;
% addpath(genpath(pwd));
% cd ..;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyb_method = [1];           % 混合预编码算法：
%[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% 3:full-analog method   4:spatially_sparse_precoding
% 5:SVD hybrid precoding  6:hybrid precoding(metting)
% 7:hybrid precoding2(metting)   8:FPS-AltMin Algorithm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_dB=[-10:2:-8];  SNR_linear=10.^(SNR_dB/10.);
N_iter=1; 
num_of_HybridMethod = length(hyb_method);
CC=zeros(length(SNR_linear),num_of_HybridMethod);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns=32;  %卫星天线数
L=9;    %卫星散射分量
Lm=10;  %基站到用户的路径数
N=8;    %FSS的个数
Np=7;   %每个基站的天线数
M=6;    %基站个数
fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    P=sqrt(Ns/L);
    sigma2=P/SNR;
    for i=1:N_iter
        for index_type = 1:length(hyb_method)
            tmptype = hyb_method(index_type);
            [H] =channel_realization_0602(Ns,L,N,lamada); %-- H(Ns*N):the channel matrix of FSS terminal
            [G] =channel_realization_0602(Np,Lm,M,lamada); %-- G(Np*M):基站到用户的信道矩阵
            switch (tmptype)
                case 1    %-- conventional method(full-connected) 
                    tao_n=tao_compute_0602(N,H,Ns,sigma2);         %(1*N)FSS终端的信干燥比
%                     F1=SVD_precoding(Nt_RF,H);
%                     temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F1*F1'*H'));
                
                otherwise
                    error('Unknown TxRx.Decoder.LDPC.Type.')
            end
        end 
    end
%     CC(i_snr,:)= real(temp/N_iter);
end
% marker_style = {'o-','o--','+-','+--','v-','v--','<-','<--','x-','x--','^-','^--','*-','*--','d-','d--','h-','h--','^-','hexagram--','diamond:','pentagram:','v:'};
% marker_color = [...
%             0.0000    0.4470    0.7410;...
%             0.8500    0.3250    0.0980;...
%             0.9290    0.6940    0.1250;...
%             0.4940    0.1840    0.5560;...
%             0.4660    0.6740    0.1880;...
%             0.3010    0.7450    0.9330;...
%             0.6350    0.0780    0.1840;...
%             0.7500    0.7500    0.0000;...
%             0.7500    0.0000    0.7500;...
%             0.0000    0.5000    0.0000;...
%             0.0000    0.0000    1.0000;...
%             0.2000    0.5000    0.6600;...
%             0.6500    0.4500    0.5500;...
%             0.6000    0.2000    0.5000;...
%             0.2000    0.8880    1.0000;...
%             0.3400    0.4300    0.6600;...
%             0.3660    0.2880    1.0000;...
%             0.6400    0.6500    0.8800;...
%         ];
% MarkerSize =7;
% LineWidth = 2.2;
% 
% figure(1)
% for i=1:num_of_HybridMethod
%     plot(SNR_dB(1:end),CC(:,i),marker_style{i}, 'color', marker_color(i,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
% %    semilogy(SNR_dB(1:end),CC(:,i),marker_style{i}, 'color', marker_color(i,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
%     hold on;grid on;
% end
% 
% legend('1','2','5','hybrid precoding(metting)','FPS-AltMin Algorithm');
% %[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% % 3:full-analog method   4:spatially_sparse_precoding
% % 5:SVD hybrid precoding  6:hybrid precoding(metting)
% % 7:hybrid precoding2(metting)   8:FPS-AltMin Algorithm]
% title(['iter = ',num2str(N_iter), ' , Nt = ',num2str(Nt),' , Nr = ',num2str(Nr),' , RF chains = ',num2str(N_RF),]);
% FontSize =10;
% xlabel('SNR (dB)')
% ylabel('Spectral efficiency (bits/s/Hz)')
% hold on