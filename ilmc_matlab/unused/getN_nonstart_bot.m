function N_nonstart_bot = getN_nonstart_bot(Ns, l, S, r)

A = getAltNg_pyr(r, max(0,l-sqrt(Ns)-r), min(S-1, floor( (l-sqrt(Ns))/r) ) );
B = -2*getAltNg_pyr(r, max(0, l-3/2*sqrt(Ns)-r), min(S-1, floor( (l-3/2*sqrt(Ns))/r) ) );

N_nonstart_bot = A+B;