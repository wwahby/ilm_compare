function [Nc, N_ns_t, N_ns_b] = getNc_vec(Nt, Ns, lmax, Mt, S, r)


[Nstart, N_ns_t, N_ns_b] = getNstart_vec(Nt,Ns,lmax,S,r);
Nc = Mt/Nstart;


end