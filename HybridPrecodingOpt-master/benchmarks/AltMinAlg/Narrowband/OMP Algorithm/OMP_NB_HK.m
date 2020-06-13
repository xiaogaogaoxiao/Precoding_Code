function [FRF, FBB, stats] = OMP_NB_HK( Fopt, NRF, At )
% [1] Omar El Ayach, Member, IEEE, Sridhar Rajagopal Spatially Sparse Precoding
% in Millimeter Wave MIMO Systems, IEEE TRANSACTIONS ON WIRELESS COMMUNICATIONS, VOL. 13, NO. 3, MARCH 2014
    FRF = [];
    Fres = Fopt;
    
    % set start time
    start_time = tic();
    
    for k = 1:NRF
        PU = At' * Fres;                     % line 4 of Alg.1 of [1]
    %     [aa1,bb1] = max(diag(PU * PU'));
        [aa,bb] = max(sum( abs(PU).^2, 2 )); % line 5 of Alg.1 of [1]
        FRF = [FRF , At(:,bb)];              % line 6 of Alg.1 of [1]
        FBB = pinv(FRF) * Fopt;  %% line 7 of Alg.1 of [1] use pseudoinverse to avoid the inverse of a possible singular matrix
        Fres = (Fopt - FRF * FBB) / norm(Fopt - FRF * FBB,'fro');% line 8 of Alg.1
    end
    
    % measure elapsed time
    elapsed_time = toc(start_time);  
    
    cost_val = norm(Fopt(:,:) - FRF * FBB,'fro')^2;    
    
    cost(1) = cost_val;
    time(1) = 0;      
    cost(2) = cost_val;
    time(2) = elapsed_time;
    
    stats.cost = cost;
    stats.time = time;       

end

