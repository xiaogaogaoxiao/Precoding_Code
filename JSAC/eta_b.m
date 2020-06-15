function eta_B = eta_b(Bj,Vrf,i,j,N)
B = 0;
for l = 1:N
    if l == i
        continue;
    end
    B = B+Bj(i,l)*Vrf(l,j);
end
eta_B = B;