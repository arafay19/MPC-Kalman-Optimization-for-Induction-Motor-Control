function [Qkf, Rkf] = Compute_Qkf_Rkf(kf,Q,S,R)
Rkf=R; % add first element to Rkf
Qkf=Q; % add first element to Qkf
for ii=0:kf-1
    if ii==kf-1
        Qkf=blkdiag(Qkf,S); % add last element to Qkf for the terminal state
    else
        Rkf=blkdiag(Rkf,R); % add element to Rkf 
        Qkf=blkdiag(Qkf,Q); % add element to Qkf 
    end
end
end

