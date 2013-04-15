function Ng_pyr = getAltNg_pyr_vec2(r,r_step,N_step)

Ng_pyr = 2*r_step.*(1+r_step).*N_step + 1/2*r.^2.*(N_step-1).*(2*N_step-1).*N_step - r*(1+2*r_step).*N_step.*(N_step-1);

end
