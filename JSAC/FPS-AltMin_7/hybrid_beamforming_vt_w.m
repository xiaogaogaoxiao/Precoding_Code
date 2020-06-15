function [V_t,V_RF,V_D] = hybrid_beamforming_vt_w(H,N_RF,V_RF,P,Nt,sigma2)
%À„∑®1 °∞[1]Hybrid Digital and Analog Beamforming Design for Large-Scale Antenna Arrays°±
    for j=1:N_RF
        V_RF1=V_RF;
        V_RF1(:,j)=[];
        F1=H'*H;
%         P=sqrt(Nr*Nt/L);
        r2=P/(Nt*N_RF);
        C=eye(N_RF-1)+(r2*sigma2)*V_RF1'*F1*V_RF1;
        G=(r2*sigma2)*F1-(r2^2*sigma2^2)*F1*V_RF1*inv(C)*V_RF1'*F1;
        for i=1:Nt
            nl(i,j)=0;
            for l=1:Nt
                if l~=i
                    nl(i,j)=nl(i,j)+G(i,l)*V_RF(l,j);  
                end
            end
            if nl(i,j)==0
                V_RF(i,j)=1;
            else
                V_RF(i,j)=nl(i,j)/abs(nl(i,j));
            end
        end
    end
    H_eff=H*V_RF;
    H_eff_Q=H_eff*(V_RF'*V_RF)^(-1/2);
    %uu=(V_RF'*V_RF)^(-1/2);
    %U_e=svd(H_eff_Q);
    [~,~,U_e]=svd(H_eff_Q);
    T_e=sqrt(P/N_RF)*eye(N_RF);
    V_D=(V_RF'*V_RF)^(-1/2)*U_e*T_e;
    V_t=V_RF*V_D;
end

