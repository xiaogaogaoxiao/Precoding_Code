function [Vrf] = hybrid_beamforming_vrf_w(N_RF,P,V_RF,H,Nt)
global Vrf z c epsilon_B epsilon_D eta_B eta_D theta1 theta2 theta_opt phi;
H_hat=inv(P^0.5)*H;
F = [1000];
while (isempty(F) || F > 1e-6)   %判断是否收敛
    for j = 1:N_RF
        Vrfj=V_RF;
        Vrfj(:,j)=[];%删除Vrf的第j列后的矩阵
        H_Vrfj=H*Vrfj;
        Aj = inv(P^0.5)*H_Vrfj*H_Vrfj'*inv(P^0.5);  %Algorithm 3的第3行
        Bj = H_hat'*inv(Aj^2)*H_hat;  %APPENDIX A中的矩阵Bj
        Dj = H_hat'*inv(Aj)*H_hat;    %APPENDIX A中的矩阵Dj
        for i =1: Nt                     
            epsilon_B(i,j) = epsilon_b(Bj,V_RF,i,j,Nt);%第5行
            epsilon_D(i,j) = epsilon_d(Dj,V_RF,i,j,Nt);%第5行
            eta_B(i,j) = eta_b(Bj,V_RF,i,j,Nt);        %第5行
            eta_D(i,j) = eta_d(Dj,V_RF,i,j,Nt);        %第5行
            %公式（28）上面的z(i,j)
            z(i,j) = imag(2*conj(eta_B(i,j))*eta_D(i,j));
            %公式（28）上面的c(i,j)
            c(i,j) = (1+epsilon_D(i,j))*eta_B(i,j)-epsilon_B(i,j)*eta_D(i,j);
            if real(c(i,j)) >= 0     %公式（28）
                phi(i,j) = asin(imag(c(i,j))/abs(c(i,j)));
            else
                phi(i,j) =pi-asin(imag(c(i,j))/abs(c(i,j)));
            end
            theta1(i,j) = -phi(i,j) + asin(z(i,j)/abs(c(i,j)));%公式（27a）
            theta2(i,j) = pi-phi(i,j) - asin(z(i,j)/abs(c(i,j)));%公式（27b）
            if theta1(i,j)<0      %将公式（27）求解的两个角度转换在[0,2pi)上
                theta1(i,j)=theta1(i,j)+2*pi;
            elseif theta1(i,j)>2*pi
                theta1(i,j)=theta1(i,j)-2*pi;
            end
            if theta2(i,j)<0
                theta2(i,j)=theta2(i,j)+2*pi;
            elseif theta2(i,j)>2*pi
                theta2(i,j)=theta2(i,j)-2*pi;
            end
            %使用公式（26）求解函数值
            f1 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(theta1(i,j))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(theta1(i,j))*eta_D(i,j)));
            f2 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(theta2(i,j))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(theta2(i,j))*eta_D(i,j)));
%             f1 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(exp(-1i*theta1(i,j)))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(exp(-1i*theta1(i,j)))*eta_D(i,j)));
%             f2 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(exp(-1i*theta2(i,j)))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(exp(-1i*theta2(i,j)))*eta_D(i,j)));
            if f1>=f2     %公式（29）
                theta_opt(i,j)= theta2(i,j);
            else
                theta_opt(i,j)= theta1(i,j);
            end
            Vrf(i,j)=exp(-1i*theta_opt(i,j));  %Algorithm 3第8行
        end
    end
    F=F-norm(Vrf,'fro')^2;   
    V_RF=Vrf;
end
end
% %     z = zeros(N_RF,Nt);c = zeros(N_RF,Nt);
%     H_hat=inv(P^0.5)*H;
%     for j = 1:N_RF
%         Vrfj = [V_RF(:,1:j-1) V_RF(:,j+1:N_RF)]; %删除Vrf的第j列后的矩阵
%         Aj = H_hat*Vrfj*Vrfj'*H'*inv(P^0.5);
%         Bj = H_hat'*inv(Aj^2)*H_hat;
%         Dj = H_hat'*inv(Aj)*H_hat;
%         for i =1: Nt                     
%             epsilon_B(i,j) = epsilon_b(Bj,V_RF,i,j,Nt);
%             epsilon_D(i,j) = epsilon_d(Dj,V_RF,i,j,Nt);
%             eta_B(i,j) = eta_b(Bj,V_RF,i,j,Nt);
%             eta_D(i,j) = eta_d(Dj,V_RF,i,j,Nt);
%             z(i,j) = imag(2*eta_B(i,j)*eta_D(i,j));
%             c(i,j) = (1+epsilon_D(i,j))*eta_B(i,j)-epsilon_B(i,j)*eta_D(i,j);
%             if real(c(i,j)) >= 0
%                 phi(i,j) = asin(imag(c(i,j))/abs(c(i,j)));
%             else
%                 phi(i,j) =pi-asin(imag(c(i,j))/abs(c(i,j)));
%             end
%             theta1(i,j) = -phi(i,j) + asin(z(i,j)/abs(c(i,j)));
%             theta2(i,j) = pi-phi(i,j) - asin(z(i,j)/abs(c(i,j)));
%             f1 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(theta1(i,j))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(theta1(i,j))*eta_D(i,j)));
%             f2 = Nt*trace(inv(Aj))-Nt*(epsilon_B(i,j)+2*real(conj(theta2(i,j))*eta_B(i,j)))/(1+epsilon_D(i,j)+2*real(conj(theta2(i,j))*eta_D(i,j)));
%             if f1>=f2
%                 theta_opt(i,j)= theta2(i,j);
%             else
%                 theta_opt(i,j)= theta1(i,j);
%             end
%             Vrf(i,j)=exp(-1i*theta_opt(i,j));
%         end
%     end
  



