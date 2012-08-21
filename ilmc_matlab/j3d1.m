% Calculate the 3D IDF using Joyner's method
%Np = 1e3; % number of points for plotting
clear all
close all

Ng = 16e6; % number of gates
N_strata = 4;

Nt = Ng; % total number of gates
S = N_strata; % number of strata
Ns = Nt/S; % Avg number of gates per stratum
r = 1; % strata to gate pitch ratio

lmax_2d = 2*sqrt(Ng); % max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + N_strata*r; % max 3D wirelength (assuming optimal routing)
lmax = lmax_3d;


%k = 4; %rent constant
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha;



[Iidf Iexp Mt Nb Nc Nstart] = getIidf(Nt,lmax,S,r,alpha,k,p);

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




