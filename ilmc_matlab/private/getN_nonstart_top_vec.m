function Ns_top = getN_nonstart_top_vec(Ns, lmax)

l = 0:lmax;
Ns_top = zeros(1,lmax+1);

l1_ind = floor(sqrt(Ns)/2)+1;
l2_ind = floor(sqrt(Ns))+1;
l3_ind = floor(3*sqrt(Ns)/2)+1;
l4_ind = 2*sqrt(Ns)+1;

l1 = 1:l1_ind;
l2 = (l1_ind+1):l2_ind;
l3 = (l2_ind+1):l3_ind;
l4 = (l3_ind+1):lmax+1;

Ns_top(1:l1_ind) = l1;
Ns_top(l1_ind+1:l2_ind) = l2 + (l2 - sqrt(Ns)/2 -1).*(l2 - sqrt(Ns)/2);
Ns_top(l2_ind+1:l3_ind) = l3*sqrt(Ns) - 3*Ns/4 + sqrt(Ns)/2;
Ns_top(l3_ind+1:end) = Ns - (2*sqrt(Ns)-l4).*(2*sqrt(Ns)-l4-1);

end