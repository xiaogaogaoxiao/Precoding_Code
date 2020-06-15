function [h,power_matrix,A_BS,A_MS,Fopt,Wopt]=mmWave_channel(Nr,Nt,Ns,L,lamada)
%蜂窝网信道矩阵的实现
count = 0;
power=sqrt(Nr*Nt/L);
alpha=(randn(1,L)+1i*randn(1,L))/sqrt(2);
power_matrix=power*diag(alpha);
AoA=2*pi*rand(1,L)-pi;   %到达角
AoD=2*pi*rand(1,L)-pi;   %偏离角
% AoA=(pi/3)*rand(1,L)-pi/6;
% AoD=(pi/3)*rand(1,L)-pi/6;
d=lamada/2;
for l=1:L
    A_BS(:,l)=array_respones(AoD(l),Nt,d,lamada);    %基站的阵列相应
    A_MS(:,l)=array_respones(AoA(l),Nr,d,lamada);    %移动设备的阵列相应
end
h=A_MS*power_matrix*A_BS';    %信道矩阵
% if(rank(h(:,:))>=Ns)
%     count = count + 1;
    [U,S,V] = svd(h(:,:));
    Fopt(:,:) = V([1:Nt],[1:Ns]);
    Wopt(:,:) = U([1:Nr],[1:Ns]);
% end