function Ng_pyr = getAltNg_pyr(r,r_step,N_step)
% using this function because Joyner has a typo somewhere 

Ng_pyr = 0;

for v_step = 0:N_step-1
    a = r_step-v_step*r;
    Ng_step = 2*a*(a+1);
    Ng_pyr = Ng_pyr + Ng_step;
end