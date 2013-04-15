function [ j2d ] = joyner_lin_vec(lmax,Nt)
l = 0:lmax;
j2d = zeros(1,lmax+1);

% define the start and end points of the two slices of l
l1min = 1;
l1max = floor(sqrt(Nt))-1;
l2min = l1max+1;
l2max = lmax;

l1 = l1min:l1max;
l2 = l2min:l2max;

l1_ind = l1+1;
l2_ind = l2+1;

j2d(1) = 1/Nt;

% includes the division by 2 to correct for overcounting
j2d(l1_ind) = 1/2/Nt^2 * (4*Nt*l1 - 4*sqrt(Nt)*l1.^2 + 2/3*l1.^3);
j2d(l2_ind) = 1/3/Nt^2 * (2*sqrt(Nt) - l2).^3;

end

