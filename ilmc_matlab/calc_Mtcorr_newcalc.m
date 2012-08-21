function [Mt term2 term3 term2_conv term3_conv] = calc_Mtcorr_newcalc(Nt,Nx,t,w,T,Ntsvs_1d)

%% setup
lmax = 2*Nx;

%% translate

Ntsvs = Ntsvs_1d;
lxmax = Nx;
lymax = Nx; % symmetric

%% precalculate g and h
% we'll need these to speed up the calculation for all l

g = zeros(1,lxmax+1);
h = zeros(1,lxmax+1);

for lx=0:lxmax
    lx_ind = lx+1;
    g(lx_ind) = calc_g(lx,t,w,T,Ntsvs,lxmax);
    h(lx_ind) = calc_h(lx,t,w,T,Ntsvs,lxmax);
end

%% Calculate!

j2d = zeros(1,lmax+1);
f3d = zeros(1,lmax+1);
term2 = zeros(1,lmax+1);
term3 = zeros(1,lmax+1);

for l = 0:lmax
    l_ind = l+1;
    j2d(l_ind) = joyner_lin(l,Nt);
    
    if(l < Nx)
        lxmax_inner = l;
        lxmin_inner = 0;
    else
        lxmax_inner = lxmax;
        lxmin_inner = l-Nx;
    end
    
    for lx = 0:lxmax_inner
        lx_ind = lx+1;
        ly_ind = l-lx+ 1;
        
        if (lx_ind > lxmax_inner+1) || (lx_ind < lxmin_inner+1)
            gx = 0;
            hx = 0;
        else
            gx = g(lx_ind);
            hx = h(lx_ind);
        end
        
        if (ly_ind > lxmax_inner +1) || (ly_ind < lxmin_inner + 1) % using lxmax and lxmin because we assume symmetric chip
            gy = 0;
            hy = 0;
        else
            gy = g(ly_ind);
            hy = h(ly_ind);
        end
        
        term2(l_ind) = term2(l_ind) + gx*gy;
        term3(l_ind) = term3(l_ind) + hx*hy;
        %term3(l_ind) = term3(l_ind) +
        %calc_h(lx,t,w,T,Ntsvs,lxmax)*calc_h(l-lx,t,w,T,Ntsvs,lxmax);
    end
    
    if l_ind == 1
        j2dnum = j2d(l_ind);
    else
        j2dnum = j2d(l_ind)/2;
    end
    
    f3d(l_ind) = j2dnum - 1/Nt^2*(term2(l_ind) + term3(l_ind));
end
 
term2_conv = conv(g,g);
term3_conv = conv(h,h);

Mt = f3d*Nt^2;
%Mt = f3d*Nt^2/2;
%Mt(1) = Mt(1)*2;

