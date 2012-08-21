function Mt = getMt(Ns,l,r,S,Mt_intra_corr,use_corrected)

Mt = 0;
for v = 0:S-1
    a = S-v;
    b = 2-dd(l-v*r)-dd(v);
    c = getMt_int_ra(Ns,l-v*r,Mt_intra_corr,use_corrected);
    Mt = Mt + a*b*c;
    %Mt = Mt + (S-v)*(2-dd(l-v*r)-dd(v) ).* Mt_int_ra(l-v*r);
end
