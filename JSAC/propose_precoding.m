function F=propose_precoding(Nt_RF,Nt,Nr,M,H,SNR)
F=zeros(Nt,Nt_RF);
for i=1:Nt_RF
    G=H'*inv(eye(Nr)+SNR*H*F(:,1:(i-1))*F(:,1:(i-1))'*H')*H;
    f=zeros(Nt,1);
    temp=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M);
    [~,~,V]=svd(temp);
    op=V(:,1);
    phase=exp(1i*angle(op))/sqrt(M);
    a=(op'*phase+phase'*op)/(2*phase'*phase);
    f(M*(i-1)+1:M*(i-1)+M)=a*phase;
    F(:,i)=f;
end