function F=SVD_precoding(Nt_RF,H)
[~,~,V]=svd(H);
F=V(:,1:Nt_RF);