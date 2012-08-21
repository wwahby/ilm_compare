% test stuff
clear all
close all
use_corrected = 1;
Ng = 16e6; % number of gates
N_strata = 4;

Nt = Ng; % total number of gates
S = N_strata; % number of strata
Ns = Nt/S; % Avg number of gates per stratum
r = 1; % strata to gate pitch ratio

lmax_2d = 2*sqrt(Ns); % max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + (N_strata-1)*r; % max 3D wirelength (assuming optimal routing)
lmax = lmax_3d;


%k = 4; %rent constant
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha;

Ns = Nt/S;

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
max_area_ratio = 0.10; % ratio of total TSV area to chip area
N_tsv_gates = Ns* max_area_ratio/(1-max_area_ratio);
w = 100;

max_length_ratio = sqrt(max_area_ratio);

N_tsv_gates_1d = floor(sqrt(N_tsv_gates));
N_tsvs_1d = floor(N_tsv_gates_1d/w);

% corrected gate counts (length of chip in gates, as opposed to actual
% number of gates)
Nsc = Ns + N_tsv_gates;
Nxc = floor(sqrt(Nsc));

T = floor(Nxc/N_tsvs_1d);
t = 1/2*(T-w);

Mt_intra_corr = calc_Mtcorr_newcalc(Ns,2*sqrt(Ns),t,w,T,N_tsvs_1d);

%% calculate everything!
for i=1:length(l)
    Mt(i) = getMt(Ns,l(i),r,S,Mt_intra_corr,use_corrected);
    
    Nb(i) = sum(Nc) - Nc(1); % want the sum of every Nc from l=1 to l=i-1
    
    N_nonstart_top(i) = getN_nonstart_top(Ns, l(i));
    N_nonstart_bot(i) = getN_nonstart_bot(Ns, l(i), S, r);
    Nstart(i) = Nt - N_nonstart_top(i) - N_nonstart_bot(i);
    
    Nc(i) = getNc(Nt, Ns, l(i), Mt(i), S, r);

    Iexp(i) = getIexp(alpha,k,p, Nb(i), Nc(i) );
    Iidf(i) = Iexp(i)*Mt(i);
    
end

l = 0:lmax_3d;

figure(1)
clf
loglog(l,Iidf)
xlabel('length (gates)')
ylabel('Interconnect density')
title('I_{idf}')

figure(2)
clf
loglog(l,Mt)
xlabel('length (gates)')
ylabel('Number of gate pairs')
title('Mt')

figure(3)
clf
loglog(l,Iexp)
xlabel('length (gates)')
ylabel('2d interconnect distribution')
title('I_{exp}')

figure(4)
clf
loglog(l,Nb)
xlabel('length (gates)')
ylabel('Nb')
title('Nb')

figure(5)
clf
loglog(l,Nc)
xlabel('length (gates)')
ylabel('Nc')
title('Nc')

figure(6)
clf
loglog(l,Nstart)
xlabel('length (gates)')
ylabel('Nstart')
title('Nstart')
