function [ H,At,Ar,Fopt,Wopt ] = channel_realization_0602(Ns,L,N,lamada)
%      power=sqrt(Nr*Nt/L);
    gamma = sqrt(Ns/L); %normalization factor
    sigma = 1; %according to the normalization condition of the H
    angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
    d=lamada/2;
    
    AoD_m = unifrnd(-pi/2,pi/2);
    AoD(1,[1:L]) = laprnd(1,L,AoD_m(1),angle_sigma);
    H(:,:) = zeros(Ns,N);
    for n = 1:N
        for L_ray = 1:L                    
            At(:,L_ray) = array_respones_0602(AoD(1,L_ray),Ns,d,lamada); %UPA array response
            alpha(L_ray) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
            H(:,n) = H(:,n) + alpha(L_ray) * At(:,L_ray);
        end
    end
    H(:,:) = gamma * H(:,:);
end