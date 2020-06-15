function a=array_respones_0602(azimuth,N,d,lamada)
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1i*[0:N-1]*2*pi*d*sin(azimuth(i))/lamada)).'];
end