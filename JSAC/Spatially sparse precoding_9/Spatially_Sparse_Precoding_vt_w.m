function [F,F_BB] = Spatially_Sparse_Precoding_vt_w(F_opt,Nt_RF,At,Ns)
F=[];
% [~,L]=size(At);
% [~,~,V]=svd(H);
% F_opt=V(:,1:Nt_RF);
F_res=F_opt;
for i=1:Nt_RF
    pha=At'*F_res;
    [~,k]=max(diag(pha*pha'));
    F=[F At(:,k)];
    F_BB=inv(F'*F)*F'*F_opt;
    F_res=(F_opt-F*F_BB)/norm(F_opt-F*F_BB,'fro');
end
F_BB=sqrt(Ns)*F_BB/norm(F*F_BB,'fro');