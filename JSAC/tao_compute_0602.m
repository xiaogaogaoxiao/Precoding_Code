function a=tao_compute_0602(N,H,Ns,sigma2)
tao_n=zeros(1,N);
w=zeros(Ns,N);%ÎÀÐÇµÄ²¨ÊøÐÎ³É¾ØÕó
for n = 1:N
    a=rand(Ns,1);
    b=rand(Ns,1);
    w(:,n)=(1/sigma2)*complex(a,b);
    R=H(:,n)*H(:,n)';
end
for i = 1:N
    I_int=;
end
    
    I_ext=;
    I_AN=;
    tao_n_fenzi=w'*R*w;
    tao_n_fenmu=w'*R*w;
    tao_n(n)=1;
    za=0;
end
P=H(:,1)*H(:,1)';
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1i*[0:N-1]*2*pi*d*sin(azimuth(i))/lamada)).'];
end