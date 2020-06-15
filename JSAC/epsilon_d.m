function [epsilon_D] = epsilon_d(Dj,Vrf,i,j,N)
Vrf1 = conj(Vrf);
C =0;
for m = 1:N
    if m == i
        continue;
    end
    for n = 1:N
        if n == i
            continue;
        end
        C = C+Vrf1(m,j)*Dj(m,n)*Vrf(n,j);
    end
end
epsilon_D = Dj(i,i)+2*real(C);