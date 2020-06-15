function F=full_analog(Nt_RF,Nt,Nr,M,H,SNR)
F=zeros(Nt,Nr);
for i=1:Nt_RF
    G=H'*inv(eye(Nr)+SNR*H*F(:,1:(i-1))*F(:,1:(i-1))'*H')*H;
    f=zeros(Nt,1);
    for j=1:M
        temp=G([1:j-1,j+1:length(G)],[1:j-1,j+1:length(G)]);
        f(M*(i-1)+j)=(1/sqrt(Nt_RF))*exp(1i*angle(temp(j,:)*[zeros(M*(i-1),1); f([(M*(i-1)+1):(M*(i-1)+j-1),(M*(i-1)+j+1):(M*(i-1)+M)]);  zeros(M*(Nt_RF-i),1)]));
    end
    F(:,i)=f;
end