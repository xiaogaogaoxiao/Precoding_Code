function [W,W_BB] = Spatially_Sparse_Precoding_wt_w(W_MMSE,Nr_RF,Ar,y,Nr)
W=[];
% [~,L]=size(At);
% [~,~,V]=svd(H);
% F_opt=V(:,1:Nt_RF);
W_res=W_MMSE;
for i=1:Nr_RF
    pha=Ar'*mean(y*y','all')*W_res;
    [~,k]=max(diag(pha*pha'));
    W=[W Ar(:,k)];
    W_BB=inv(W'*mean(y*y','all')*W)*W'*mean(y*y','all')*W_MMSE;
    W_res=(W_MMSE-W*W_BB)/norm(W_MMSE-W*W_BB,'fro');
end
W_BB=sqrt(Nr)*W_BB/norm(W*W_BB,'fro');