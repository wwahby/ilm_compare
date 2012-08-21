function Nb = getNb(Nc)

Nb = zeros(1,length(Nc));
for l=0:length(Nc)
    Nb(l) = sum(Nc(0:l));
end
