function Mt = getMt_vec(Ns,lmax,r,S,Mt_intra)

Mt = zeros(1,lmax+1);

%lmax_2d = lmax;
lmax_2d = length(Mt_intra)-1;

for vind = 1:S
    v = vind-1;
    a = S-v;
    b = 2;
    c = zeros(1,lmax+1);
    
    if (v == 0)
        b = b - 1; % accounts for the delta[v] term
        c(1:lmax_2d+1) = Mt_intra;
    else
        startInd = v*r+1;
        stopInd = lmax_2d + startInd;
        c(startInd:stopInd) = Mt_intra;
    end
    
    c = a*b*c;
    c(v*r+1) = c(v*r+1)/b*(b-1); % accounts for the delta(l-vr) term. +1 because of matlab indexing
    
    Mt = Mt + c;

end