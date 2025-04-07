function [W_cal_x, Omega_x, W_cal_u, Omega_u] = Compute_constraints(kf,W_x,W_xf,W_u,omega_x,omega_xf,omega_u)
W_cal_u=W_u; % add first element to condensed W_u 
Omega_u=omega_u; % add first element to condensed Omega_u

W_cal_x=W_x; % add first element to condensed W_x 
Omega_x=omega_x; % add first element to condensed Omega_x

for ii=0:kf-1
    if ii==kf-1
        W_cal_x=blkdiag(W_cal_x,W_xf); % add last element to condensed W_x for the terminal state
        Omega_x=[Omega_x;omega_xf]; % add last element to condensed Omega_x for the terminal state
    else
        W_cal_u=blkdiag(W_cal_u,W_u); % add element to condensed W_u 
        Omega_u=[Omega_u;omega_u]; % add element to condensed Omega_u
        
        W_cal_x=blkdiag(W_cal_x,W_x); % add element to condensed W_x 
        Omega_x=[Omega_x;omega_x]; % add element to condensed Omega_x
    end
end
end

