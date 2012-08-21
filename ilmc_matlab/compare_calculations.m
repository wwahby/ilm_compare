function [Iidf Mt Iexp l T t Nsc Nxc N_tsv_gates N_tsvs_1d Nc Nb N_nonstart_top N_nonstart_bot t2 t3 t2c t3c] = compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w)

Nt = Ng; % total number of gates
S = N_strata; % number of strata
Ns = Nt/S; % Avg number of gates per stratum
Nx = sqrt(Ns);

lmax_2d = 2*sqrt(Ns); % max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + (N_strata-1)*r; % max 3D wirelength (assuming optimal routing)
lmax = lmax_3d;

l = 0:lmax;
Nb= zeros(1,length(l));
Nc= zeros(1,length(l));
Nstart= zeros(1,length(l));
N_nonstart_top= zeros(1,length(l));
N_nonstart_bot= zeros(1,length(l));
Iidf = zeros(1,length(l));
Mt = zeros(1,length(l));
Iexp = zeros(1,length(l));


%% get corrected mt_intra
N_tsv_gates = Ns* max_area_ratio/(1-max_area_ratio);

N_tsv_gates_1d = floor(sqrt(N_tsv_gates));
N_tsvs_1d = floor(N_tsv_gates_1d/w);

% corrected gate counts (length of chip in gates, as opposed to actual
% number of gates)
Nsc = ceil(Ns + N_tsv_gates);
Nxc = floor(sqrt(Nsc));
Ntc = S*Nsc;

T = floor(Nxc/N_tsvs_1d);
t = 1/2*(T-w);

if(use_corrected == 1)
    [Mt_intra_corr t2 t3 t2c t3c] = calc_Mtcorr_newcalc(Nsc,Nxc,t,w,T,N_tsvs_1d);
else
    Mt_intra_corr = 0;  
    Nxc = Nx;
    Nsc = Ns;
    Ntc = Nt;
end

%% calculate everything!
for i=1:length(l)
    Mt(i) = getMt(Nsc,l(i),r,S,Mt_intra_corr,use_corrected);
    
    Nb(i) = sum(Nc) - Nc(1); % want the sum of every Nc from l=1 to l=i-1
    
    N_nonstart_top(i) = getN_nonstart_top(Nsc, l(i));
    N_nonstart_bot(i) = getN_nonstart_bot(Nsc, l(i), S, r);
    Nstart(i) = Nt - N_nonstart_top(i) - N_nonstart_bot(i);
    
    Nc(i) = getNc(Ntc, Nsc, l(i), Mt(i), S, r);

    Iexp(i) = getIexp(alpha,k,p, Nb(i), Nc(i) );
    Iidf(i) = Iexp(i)*Mt(i);
    
end

%% Cleanup for returning stuff
l = 0:lmax_3d;

