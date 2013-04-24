function Ng_pyr = getNg_pyr(r,r_step,N_step)

A = r^2/6*(N_step-1)*N_step*(2*N_step-1);
B = -r/2*(2*r_step+1)*(N_step-1)*N_step;
C = N_step*(r_step^2 + r_step);
Ng_pyr = 2*(A + B + C);