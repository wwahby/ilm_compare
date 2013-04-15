function N_nonstart_bot = getN_nonstart_bot_vec(Ns, lmax, S, r)
l = 0:lmax;

a1 = max(0,l-sqrt(Ns)-r);
a2 = min(S-1, floor( (l-sqrt(Ns))/r));
a2 = max(0,a2); % need to ensure a2 doesn't go below 0 or else the compact form for Ngpyr doesn't hold

b1 = max(0,l-3/2*sqrt(Ns)-r);
b2 = min(S-1,floor( (l-3/2*sqrt(Ns))/r));
b2 = max(0,b2); % need to ensure b2 doesn't go below 0 or else the compact form for Ngpyr doesn't hold

A = getAltNg_pyr_vec2(r, a1, a2);
B = -2*getAltNg_pyr_vec2(r,b1,b2);

N_nonstart_bot = A + B;

end
