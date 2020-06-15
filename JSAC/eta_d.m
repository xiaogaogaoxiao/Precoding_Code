function eta_D = eta_d(Dj,Vrf,i,j,N)
D = 0;
for l = 1:N
    if l == i
        continue;
    end
    D = D+Dj(i,l)*Vrf(l,j);
end
eta_D = D;