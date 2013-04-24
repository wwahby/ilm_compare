function g = calc_g(lx,t,w,T,Ntsvs,lxmax)
% lx - 1d tsv length
% t - TSV offset within cell (in gate lengths)
% w - TSV width (in gate lengths)
% T - cell period (in gate lengths)
% Ntsvs - number of tsvs in 1d

g = calc_h(lx,T-t-w,w,T,Ntsvs,lxmax);
