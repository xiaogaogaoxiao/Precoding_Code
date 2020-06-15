clear all
clc   

%%  保存的结果路径
% filenameTmp = [date,'-FrameNum-',num2str(Frame_Num),'__R-',num2str(R_range),'__A-',num2str(A_range),'__BG',num2str(BG),'__',Modulation,'__AWGN_iterations-',num2str(iterations)];
filenameTmp = ['1_',date];
res.mkdir_str=strcat('.\Simulation\',filenameTmp);
mkdir(res.mkdir_str);%一运行就会在当前文件夹下创建simulation文件夹
res.mkdir_str1 =strcat(res.mkdir_str,'\');
filenameTmp1= filenameTmp;%配置文件的保存的名字
res.mkdir_str =strcat(res.mkdir_str1,filenameTmp1);
res.mkdir_str =strcat(res.mkdir_str,'.m');
res.Save_link1=Save_link_file('main.m',res.mkdir_str);%
res.result_name = 'Spectral efficiency.txt';
res.result_dir = strcat(res.mkdir_str1,res.result_name);
res.fp = fopen(res.result_dir,'w+');

% addpath('hybrid precoding(point-to-point MIMO)_5');
% addpath('FPS-AltMin_7');
% addpath('Spatially sparse precoding_9');
%把当前matlab下面的子目录全部加入到当前工作环境中
fompath = fileparts(mfilename('fullpath'));
addpath(genpath(fompath));
global Nc Nray Vrf z c epsilon_B epsilon_D eta_B eta_D theta1 theta2 theta_opt phi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyb_method = [9 10];           % 混合预编码算法：
%[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% 3:full-analog method   4:SVD hybrid precoding
% 5:hybrid precoding(metting)  6:hybrid precoding2(metting)
% 7:FPS-AltMin Algorithm]   8:Algorithms 3   9:spatially_sparse_precoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   :2:6]
SNR_dB=[-40:5:0];  SNR_linear=10.^(SNR_dB/10.);
par.mod = 'QPSK';

Nt=64;Nr=16;Nc = 8;Ns=2;%Ns_num = 2;
Nt_RF=6;N_RF=4;Nr_RF=4;
N_iter=100; 

% K = 8; % K个单天线用户
beta_K = ones(1,Nr);
P_K = eye(Nr);
V_RF=eye(Nt,N_RF);
V_RF2=eye(Nt,N_RF);
VD=eye(N_RF,Nr);
W_RF=eye(Nr,N_RF);
V_RF1=eye(Nt,N_RF);
W_RF1=eye(Nr,N_RF);
F_RF=eye(Nt,Nt_RF);
Nray=10; % number of rays
M=Nt/Nt_RF; % number of antennas connected to one RF chains
fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;
cc = 1/sqrt(Nc)*exp(1i*[0:2*pi/Nc:2*pi-2*pi/Nc])';
C = kron(eye(N_RF),cc);
num_of_HybridMethod = length(hyb_method);
CC=zeros(length(SNR_linear),num_of_HybridMethod);
z = zeros(Nt,N_RF);c = zeros(Nt,N_RF);
epsilon_B = zeros(Nt,N_RF);epsilon_D = zeros(Nt,N_RF);eta_B = zeros(Nt,N_RF);
eta_D = zeros(Nt,N_RF);theta1 = zeros(Nt,N_RF);theta2 = zeros(Nt,N_RF);
theta_opt = zeros(Nt,N_RF);phi = zeros(Nt,N_RF);Vrf=zeros(Nt,N_RF);
[par]=symbol_w(par);
% for Ns_idx=1:length(Ns_num)
%     Ns=Ns_num(Ns_idx);

% -- start simulation
% track simulation time

for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    temp=zeros(1,num_of_HybridMethod);
    R=0;
    P=sqrt(Nr*Nt/Nray);
    sigma2=P/SNR;
%     sigma2=SNR/16;
    for ii=1:N_iter
        for index_type = 1:length(hyb_method)
            tmptype = hyb_method(index_type);
%             [H,power_matrix,A_BS,~,Fopt,Wopt]=mmWave_channel(Nr,Nt,Ns,Nray,lamada);
            [H,At,Ar,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns,lamada);
            % generate random bit stream
            res.b = randi([0 1],Ns,par.bps);
            % generate transmit symbols
            res.idx = bi2de(res.b,'left-msb')+1;
            s = par.symbols(res.idx).';
            % generate noise vector
            n = sqrt(0.5)*(randn(Nr,1)+1i*randn(Nr,1));
            
            if ii==1 
                [~,V_RF2,~] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2);
            end
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
                case 4    %-- SVD hybrid precoding
                    F5=hybrid_precoding(Nt_RF,Nt,Nr,M,H,SNR);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(SNR/16)*H*F5*F5'*H'));
                case 5    %-- point-to-point MIMO 场景的算法（5月24日例会讲的算法）
                    %Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays
                    [V_t,V_RF,V_D] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2);
                    [W_t] = hybrid_beamforming_wt_w(H,N_RF,W_RF,V_t,Nr,P,sigma2);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(1/sigma2)*W_t*((W_t'*W_t)^(-1))*W_t'*H*V_t*V_t'*H'));
                case 6 
                    [V_t1,V_RF1,V_D1] = hybrid_beamforming_vt2_w(H,N_RF,V_RF1,P,Nt,sigma2);
                    [W_t1] = hybrid_beamforming_wt2_w(H,N_RF,W_RF1,V_t1,Nr,P,sigma2);
                    temp(index_type)=temp(index_type)+log2(det(eye(Nr)+(1/sigma2)*W_t1*((W_t1'*W_t1)^(-1))*W_t1'*H*V_t1*V_t1'*H'));
                case 7    %-- Fixed Phase Shifters(FPS-AltMin)算法 
                    %Hybrid Precoding in Millimeter Wave Systems:How Many Phase Shifters Are Needed?
                    %(用到和速率公式)Alternating Minimization Algorithms for Hybrid Precoding in Millimeter Wave MIMO Systems
                    [ FRF, FBB ] = my_AltMin_new( Fopt, C);
                    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
                    [ WRF, WBB ] = my_AltMin_new( Wopt, C);
                    temp(index_type)=temp(index_type)+log2(det(eye(Ns) +  (P/(sigma2*Ns))* pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
                case 8    %-- Algorithms 3
                    F = [1000];
                    while (isempty(F) || F > 1e-6) %判断是否收敛
                        %Algorithm 3的第2行到第10行
                        [ V_RF2 ] = hybrid_beamforming_vrf_w(N_RF,P_K,V_RF2,H,Nt);
%                         [~,V_RF2,~] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2);
%                         [power2]= Copy_of_water_filling(P,P_K);
                        %Algorithm 3的第12行
                        [ P_K ] = hybrid_beamforming_P_w(V_RF2,H,Nr,P,sigma2,beta_K);                           
                        H_VRF2=H*V_RF2;
                        F=F-norm(V_RF2,'fro')^2;
                    end
                    %Algorithm 3的第14行
                    VD = H_VRF2'*inv(H_VRF2*H_VRF2')*(P_K^0.5); % 求解VD
                    
                    for k=1:Nr
                        zz=0;
                        for l=1:Nr
                            if(l~=k)
                                zz=zz+abs(H(k,:)*V_RF2*VD(:,l))^2;
                            else 
                                continue;
                            end
                        end
                        %公式（20）
                        temp(index_type)=temp(index_type)+log2(1+(abs(H(k,:)*V_RF2*VD(:,k))^2)/(sigma2+zz)); 
                    end
                case 9    %-- Spatially Sparse Precoding    
                    [F_RF,F_BB]=Spatially_Sparse_Precoding_vt_w(Fopt,Nt_RF,At,Ns);
                    F_RF_BB=F_RF*F_BB;
                    y=sqrt(P)*H*F_RF_BB*s+n;
%                     W_MMSE=(mean(s*y')*inv(mean(y*y','all')))';
                    W_MMSE=((1/sqrt(P))*inv(F_RF_BB'*H'*H*F_RF_BB+(sigma2*Ns/P)*eye(Ns))*F_RF_BB'*H')';
                    [W_RF,W_BB] = Spatially_Sparse_Precoding_wt_w(W_MMSE,Nr_RF,Ar,y,Nr);
                    W_RF_BB=W_RF*W_BB;
                    R_n=sigma2*W_RF_BB'*W_RF_BB;
                    WF=W_RF_BB'*H*F_RF_BB;
                    temp(index_type)=temp(index_type)+log2(det(eye(Ns)+(P/Ns)*inv(R_n)*WF*WF'));
                case 10    %-- Spatially Sparse Precoding    
%                     [F_RF,F_BB]=Spatially_Sparse_Precoding_vt_w(Fopt,Nt_RF,At,Ns);
                    F_RF_BB=Fopt;
                    y=sqrt(P)*H*F_RF_BB*s+n;
%                     W_MMSE=(mean(s*y')*inv(mean(y*y','all')))';
                    W_MMSE=((1/sqrt(P))*inv(F_RF_BB'*H'*H*F_RF_BB+(sigma2*Ns/P)*eye(Ns))*F_RF_BB'*H')';
                    [W_RF,W_BB] = Spatially_Sparse_Precoding_wt_w(W_MMSE,Nr_RF,Ar,y,Nr);
                    W_RF_BB=W_RF*W_BB;
                    [F_RF,F_BB]=Spatially_Sparse_Precoding_vt_w(Fopt,Nt_RF,At,Ns);
                    F_RF_BB=F_RF*F_BB;
                    R_n=sigma2*W_RF_BB'*W_RF_BB;
                    WF=W_RF_BB'*H*F_RF_BB;
                    temp(index_type)=temp(index_type)+log2(det(eye(Ns)+(P/Ns)*inv(R_n)*WF*WF'));
                otherwise
                    error('Unknown TxRx.Decoder.LDPC.Type.')
            end
        end 
    end
    CC(i_snr,:)= real(temp/N_iter);
end
%     CCC{Ns_idx}=CC;
% end
% marker_style = {'o-','o--','+-','+--','v-','v--','<-','<--','x-','x--','^-','^--','*-','*--','d-','d--','h-','h--','^-','hexagram--','diamond:','pentagram:','v:'};
marker_style = {'o-','+-','v-','<-','x-','^-','*-','*--','d-','d--','h-','h--','^-','hexagram--','diamond:','pentagram:','v:'};
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

% figure(1)
% for j=1:length(Ns_num)
% for i=1:num_of_HybridMethod
%     plot(SNR_dB(1:end),CCC{j}(:,i),marker_style{j}, 'color', marker_color(j,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
% %    semilogy(SNR_dB(1:end),CC(:,i),marker_style{i}, 'color', marker_color(i,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize);
%     hold on;grid on;
% end
% end

legend('Ns=1','Ns=2','5','hybrid precoding(metting)','FPS-AltMin Algorithm');
%[1:conventional method(full-connected)  2:SIC-based hybrid precoding
% 3:full-analog method   4:SVD hybrid precoding
% 5:hybrid precoding(metting)  6:hybrid precoding2(metting)
% 7:FPS-AltMin Algorithm]   8:Algorithms 3   9:spatially_sparse_precoding
title(['iter = ',num2str(N_iter), ' , Nt = ',num2str(Nt),' , Nr = ',num2str(Nr),' , RF chains = ',num2str(N_RF),]);
FontSize =10;
xlabel('SNR (dB)')
ylabel('Spectral efficiency (bits/s/Hz)')
hold on

strsave= strcat('.\Simulation\',filenameTmp,'\'); 
filenameTmp1= strcat(filenameTmp,'.mat'); % 每次可以load最新数据
strsave= strcat(strsave,filenameTmp1); 
s=['save ' strsave];% 保持.mat 文件，以后仿真结果可以再次确认,以后一定注意可以再次画图。
eval(s);
save main.mat;
fprintf('\nsimulation has finished!\n\n');