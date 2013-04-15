function Mt_intra = getMt_intra_vec(Ns,lmax)

% set up array slices
l1min = 1;
l1max = round(sqrt(Ns)) - 1;
l2min = l1max+1;
l2max = lmax; 

Mt_intra = zeros(1,lmax+1);

%+1 because of Matlab's start-at-1 indexing
l1 = l1min:l1max;
l2 = l2min:l2max;
l1_ind = (l1min:l1max)+1;
l2_ind = (l2min:l2max)+1;

Mt_intra(1) = Ns;
Mt_intra(l1_ind) = ( 2*Ns*l1 - 2*sqrt(Ns)*(l1.^2) + 1/3*l1.^3 );
Mt_intra(l2_ind) = 1/3 * (2*sqrt(Ns) - l2).^3;

end
