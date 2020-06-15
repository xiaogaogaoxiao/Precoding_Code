function [F,F_BB]=spatially_sparse_precoding(Nt_RF,H,At)
F=[];
[~,L]=size(At);
[~,~,V]=svd(H);
F_opt=V(:,1:Nt_RF);
F_res=F_opt;
for i=1:L
    pha=At'*F_res;
    [~,k]=max(diag(pha*pha'));
    F=[F At(:,k)];
    F_BB=inv(F'*F)*F'*F_opt;
    F_res=(F_opt-F*F_BB)/norm(F_opt-F*F_BB,'fro');
end
F_BB=sqrt(L-0.5)*F_BB/norm(F*F_BB,'fro');



    
    
