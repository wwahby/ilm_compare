% Calculate the 3D IDF using Joyner's method

Ng = 1e6; % number of gates
N_strata = 4;
v = 100; % interstratal pitch (in gate pitches)

Nt = Ng; % total number of gates
S = N_strata; % number of strata
Ns = Nt/S; % Avg number of gates per stratum

lmax_2d = 2*sqrt(Ng); % max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + N_strata*v; % max 3D wirelength (assuming optimal routing)

k = 4;
p = 0.5;



N_nonstart_top = getN_nonstart_top(Ns,lmax);
N_nonstart_bot = getN_nonstart_bot(Ns,lmax);
Nstart = getNstart(Nt, N_nonstart_top, N_nonstart_bot);

Mt = Mt(r,S);
Nc = getNc(Mt,Nstart);
Nb = getNb(Nc);

iexp = Iexp(alpha,k,p,Nb,Nc);

I = Iidf(iexp,Mt);