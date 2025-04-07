function [D_cal, E_cal, Ckf] = Compute_D_cal_E_cal_Ckf(kf, B_d, C_new)
[~, m]=size(B_d); % determine number of states and inputs

E_cal = eye(m); % first block entry of E_cal
D_cal = zeros(m, kf * m);
Ckf = blkdiag(C_new, C_new);

for ii=1:kf-1
    D_cal_line = [];
    for zz=1:kf
        if zz==ii
            D_cal_line = [D_cal_line, eye(m)];
        else
            D_cal_line = [D_cal_line, zeros(m, m)];
        end
    end
    
    D_cal = [D_cal; D_cal_line;];
    E_cal = [E_cal; zeros(m, m)];
    Ckf=blkdiag(Ckf, C_new);
end

end

