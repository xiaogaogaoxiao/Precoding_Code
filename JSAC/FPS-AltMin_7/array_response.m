% function y = array_response(a1,a2,N)
% for m= 0:sqrt(N)-1
%     for n= 0:sqrt(N)-1
%         y(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
%     end
% end
% y = y.'/sqrt(N);
% end
function a=array_response(azimuth,N,d,lamada)
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1i*[0:N-1]*2*pi*d*sin(azimuth(i))/lamada)).'];
end