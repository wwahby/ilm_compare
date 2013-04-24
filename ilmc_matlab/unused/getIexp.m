function iexp = getIexp(alpha,k,p,Nb,Nc)

iexp = alpha*k/Nc*( (1+Nb)^p + (Nb + Nc)^p - (1+Nb+Nc)^p - Nb^p);