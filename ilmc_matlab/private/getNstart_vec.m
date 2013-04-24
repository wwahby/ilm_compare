function [Nstart, N_nonstart_top, N_nonstart_bot] = getNstart_vec(Nt, Ns, lmax, S, r)
N_nonstart_top = getN_nonstart_top_vec(Ns, lmax);
N_nonstart_bot = getN_nonstart_bot_vec(Ns, lmax, S, r);

Nstart = Nt - N_nonstart_top - N_nonstart_bot;
end