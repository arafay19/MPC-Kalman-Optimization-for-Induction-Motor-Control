function [Akf, Bkf] = Compute_Akf_Bkf(kf,A_d,B_d)
[n, m]=size(B_d); % determine number of states and inputs

Akf=eye(n); % first block entry of Akf
Bkf=zeros(n,m*kf); % first block line of Bkf

for ii=1:kf
    Akf=[Akf;A_d^ii]; % add block entry to Akf
    Bkf_line=A_d^(ii-1)*B_d; % block line of Bkf
    for zz=2:ii
        Bkf_line=[Bkf_line, A_d^(ii-zz)*B_d]; % add block entry to block line
    end
    Bkf_line=[Bkf_line, zeros(n,m*(kf-ii))]; % fill des rest of the block line with zeros
    Bkf=[Bkf;Bkf_line]; % add block line to Bkf
end

end

