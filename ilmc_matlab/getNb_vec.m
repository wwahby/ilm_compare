function Nb = getNb_vec(Nc)
Nb = zeros(1,length(Nc));

for l=2:length(Nc)-1
    Nb(l) = Nb(l-1) + Nc(l-1);
end



end