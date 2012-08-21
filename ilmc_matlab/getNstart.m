function Nstart = getNstart(Nt, Ns, l, S, r)
N_nonstart_top = getN_nonstart_top(Ns, l);
N_nonstart_bot = getN_nonstart_bot(Ns, l, S, r);
Nstart = Nt - N_nonstart_top - N_nonstart_bot;
