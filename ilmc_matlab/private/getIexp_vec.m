function Iexp = getIexp_vec(alpha,k,p,Nb,Nc)

Iexp = alpha*k./Nc.*( (1+Nb).^p + (Nb + Nc).^p - (1+Nb+Nc).^p - Nb.^p);


end