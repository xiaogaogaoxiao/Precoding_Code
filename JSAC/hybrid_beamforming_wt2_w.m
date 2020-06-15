function [W_t] = hybrid_beamforming_wt_w(H,N_RF,W_RF,V_t,Nr,P,sigma2)
%À„∑®2 °∞[2]Hybrid Analog and Digital Beamforming for mmWave OFDM Large-Scale Antenna Arrays°±
    for j=1:N_RF
        W_RF1=W_RF;
        W_RF1(:,j)=[];
        F2=H*V_t*V_t'*H';
%         P=sqrt(Nr*Nt/L);
        r2=N_RF/Nr;
        C1=eye(N_RF-1)+(r2/sigma2)*W_RF1'*F2*W_RF1;
        G1=(r2/sigma2)*F2-(r2^2/sigma2^2)*F2*W_RF1*inv(C1)*W_RF1'*F2;
        for i=1:Nr
            n2(i,j)=0;
            for l=1:Nr
                if l~=i
                    n2(i,j)=n2(i,j)+G1(i,l)*W_RF(l,j);  
                end
            end
            if n2(i,j)==0
                W_RF(i,j)=1;
            else
                W_RF(i,j)=n2(i,j)/abs(n2(i,j));
            end
        end
    end
%     J=W_RF'*H*V_RF*V_D*V_D'*V_RF'*H'*W_RF+sigma2*W_RF'*W_RF;
%     W_D=J^(-1)*W_RF'*H*V_RF*V_D;
    J=W_RF'*H*V_t*V_t'*H'*W_RF+sigma2*W_RF'*W_RF;
    W_D=J^(-1)*W_RF'*H*V_t;
    W_t=W_RF*W_D;
end

