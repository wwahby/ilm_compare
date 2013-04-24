function Nc = getNc(Nt, Ns, l, Mt, S, r)

%Mt = getMt(Ns,lmax,r,S);
Nstart = getNstart(Nt, Ns, l, S, r);

Nc = Mt/Nstart;
end