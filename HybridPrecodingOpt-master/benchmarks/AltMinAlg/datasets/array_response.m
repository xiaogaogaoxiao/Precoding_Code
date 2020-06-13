function y = array_response(a1,a2,N)
%[1] Xianghao Yu, et al.  "Alternating Minimization Algorithms for Hybrid
%Precoding in Millimeter Wave MIMO Systems " IEEE JOURNAL OF SELECTED TOPICS IN SIGNAL PROCESSING, VOL. 10, NO. 3, APRIL 2016
% uniform square planar array (USPA) with
for m= 0:sqrt(N)-1
    for n= 0:sqrt(N)-1
        y(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(a1)*sin(a2) + n*cos(a2) ) ); %(4) of [1]
    end
end
y = y.'/sqrt(N);
end