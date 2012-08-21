function [Iidf Iexp Mt Nb Nc Nstart] = getIidf(Nt,lmax,S,r,alpha,k,p)

Ns = Nt/S;

l = 0:lmax;
Nb= zeros(1,length(l));
Nc= zeros(1,length(l));
Iidf = zeros(1,length(l));
Mt = zeros(1,length(l));
Iexp = zeros(1,length(l));

for i=1:length(l)
    Mt(i) = getMt(Ns,l(i),r,S);
    
    Nb(i) = sum(Nc) - Nc(1); % want the sum of every Nc from l=1 to l=i-1
    
    Nc(i) = getNc(Nt, Ns, l(i), Mt(i), S, r);

    Iexp(i) = getIexp(alpha,k,p, Nb(i), Nc(i) );
    Iidf(i) = Iexp(i)*Mt(i);
    
    %% delete stuff below this
    Nstart(i) = getNstart(Nt, Ns, l(i), S, r);
end

