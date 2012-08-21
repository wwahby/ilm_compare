function N_nonstart_top = getN_nonstart_top(Ns, l)

if ((0 <= l) && (l<= sqrt(Ns)/2))
    N_nonstart_top = l;
elseif ( (sqrt(Ns)/2 < l) && (l <= sqrt(Ns) ) )
    N_nonstart_top = l + (l -sqrt(Ns)/2 -1)*(l -sqrt(Ns)/2);
elseif ( (sqrt(Ns) < l) && (l <= 3*sqrt(Ns)/2) )
    N_nonstart_top = l*sqrt(Ns) - 3*Ns/4 + sqrt(Ns)/2;
elseif (3*sqrt(Ns)/2 < l) && (l <= 2*sqrt(Ns) ) % Think J has a typo here. 2sqrt(Ns)?
    N_nonstart_top = Ns - (2*sqrt(Ns)-l)*(2*sqrt(Ns)-l-1);
else
    N_nonstart_top = 0;
end

    
