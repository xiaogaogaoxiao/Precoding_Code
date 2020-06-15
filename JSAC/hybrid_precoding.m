function F=hybrid_precoding(Nt_RF,Nt,Nr,M,H,SNR)
F=zeros(Nt,Nt_RF);
for i=1:Nt_RF
    G=H'*inv(eye(Nr)+SNR*H*F(:,1:(i-1))*F(:,1:(i-1))'*H')*H;
    f=zeros(Nt,1);
    temp=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M);
    [~,~,op]=svd(temp);
    f(M*(i-1)+1:M*(i-1)+M)=op(:,1);
    F(:,i)=f;
end