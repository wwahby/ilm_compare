%function [iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs iidf_rewire] = gen_design(chip,tsv,gate,transistor,wire)
function [chip power wire repeater] = gen_design(chip,tsv,gate,transistor,wire,simulation)

%% Unpack inputs from objects
Ng = chip.num_gates;
alpha = chip.alpha;
k = chip.rent_k;
p = chip.rent_p;
S = chip.num_layers;


Ach_m2 = chip.area_total;
chi = chip.chi;
Tclk = chip.clock_period;
alpha_t = chip.delay_constant;
a = chip.logic_activity_factor;
Vdd = chip.Vdd;

gate_pitch = chip.gate_pitch;

w_trans = transistor.gate_length;
eps_ox = transistor.oxide_rel_permittivity;
tox = transistor.oxide_thickness;
Ioff = transistor.leakage_current_per_micron;

Ro = gate.output_resistance;
N_trans_per_gate = gate.num_transistors;


use_joyner = simulation.use_joyner;
redo_wiring = simulation.redo_wiring_after_repeaters;

AR_tsv = tsv.aspect_ratio;
Atf_max = tsv.max_area_fraction;
h_tsv_m = tsv.height;

rho_m = wire.resistivity;
epsr_d = wire.dielectric_epsr;

%% Update some objects where necessary
wire.layer_area = chip.area_total/chip.num_layers;

%% Presize the chip and TSVs
Ns = Ng/S;
Lx = round(sqrt(Ns));
Ach_tier_gp = Ach_m2/gate_pitch^2/S; % Force chip to use specified area 
Ach_tier_m2 = Ach_m2/S;

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;

% Size the TSVs
h_tsv = ceil(h_tsv_m/gate_pitch);
w_tsv = ceil(h_tsv/AR_tsv);

%% Size the chip so we have a nicely divisible number of unit cells per side
if ((S == 1) || (use_joyner == 1)) % no TSVs for single layer device, or when using original Joyner model
    w_tsv = 0;
    Lxc = Lx;
    Nsc = Ns;
    Ngc = Ng;
    Nuc_1d = 1;
    N_tsvs = 0;
else
    %Nsp = floor( (1+Atf_max)*Ns );
    Nsp = floor(Ns/(1-Atf_max));
    Lxp = floor(sqrt(Nsp));
    Nsp = Lxp^2;
    Tp = ceil(w_tsv/sqrt(Atf_max));

    slack = 0.2;
    [Lxc Tc Nuc_1d gfrac_L gfrac_T] = find_LT_combination(Lxp,Tp,slack);
    N_tsvs = Nuc_1d^2;
end

Nsc = Lxc^2;
Ngc = Nsc*S;

g_tsv = (Nuc_1d*w_tsv)^2; % number of gates displaced by TSVs
Atf_act = g_tsv/Nsc;

Ns_act = Nsc - g_tsv;
Ng_act = Ns_act*S;

repstr1 = sprintf('Ng_nom: %.4g\tNg_cor: %.4g\tNg_act: %.4g\tAtf_act: %.4g',Ng,Ngc,Ng_act,Atf_act);
disp(repstr1)

%% Calculate WLD
iidf = calc_Iidf_corrected(alpha,k,p,Lx,S,h_tsv,Nuc_1d,w_tsv);
%iidf = calc_Iidf(alpha,k,p,round(sqrt(Ng)),1,h_tsv);

%% Cleanup - Get rid of NaNs
iidf(isnan(iidf)) = 0;
lmax = length(iidf) - 1;
l = 0:lmax;

chip.iidf = iidf;
chip.lengths = l;

%% Determine wire pitch and layer assignment

Ach_wla = Ach_tier_gp; % Reduce the chip area by a factor of S when we're folding a design across S layers
%[Ln pn pn_orig Nm] = wire_layer_assignment_alt(iidf,lmax,Ach_wla,chi,rho_m,epsr_d,Tclk,alpha_t);
%[Cxc Ltot Cn] = calc_total_wiring_capacitance2(pn,Ln,Nm,iidf,epsr_d,gate_pitch);

%[Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved(iidf,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area,rho_m_vec,epsr_d,alpha_t,Beta,Tclk,Rc);
wire = wla_improved(chip,wire);
[Cxc Cn] = calc_wiring_capacitance_from_area(wire);
[Cxc2 Cn2] = calc_wiring_capacitance_from_area_old(wire.pn,wire.layers_per_tier,wire.wire_area,wire.via_area,wire.dielectric_epsr);
wire.capacitance_tot = Cxc;
wire.capacitance_per_layer = Cn;

%% Power estimates

eps0 = 8.854e-12; % (F/m)

Ilk = Ioff*(w_trans*1e6);
Cox = eps_ox*eps0*w_trans^2/tox;
Co = Cox; % Need to include parasitics for realistic estimate

Nt = N_trans_per_gate * Ng;
f = 1/Tclk;

Pdyn = 1/2*a*Co*Vdd^2*f*Nt;
Plk = (1-a)*Vdd*Ilk*Nt;
Pw = 1/2*a*Cxc*Vdd^2*f;

%% Repeater insertion

Ach_ri = Ach_tier_m2;
Ainv_min = gate_pitch^2*9; % assume 3:1 W/L for nmos, 3x that for pmos
rho_xcn = rho_m;
Co = N_trans_per_gate*Cox;

pn = wire.pn;
Ln = wire.Ln;

%[iidf_rep h_vec k_vec Arep_used num_vec size_vec] = repeater_insertion(iidf,Ach_ri,Ainv_min,wire.pn,wire.Ln,Cn,rho_xcn,Ro,Co,gate_pitch);
[iidf_rep h_vec k_vec Arep_used num_vec size_vec] = repeater_insertion(iidf,Ach_ri,Ainv_min,pn,Ln,Cn,rho_xcn,Ro,Co,gate_pitch);

Co_rep = Cox*size_vec;
Ilk_rep = Ilk*size_vec;
Plk_rep_vec = (1-a)*Vdd*Ilk_rep.*num_vec;
Plk_rep = sum(Plk_rep_vec);

Pdyn_rep_vec = 1/2*a*Vdd^2*f.*num_vec.*Co_rep;
Pdyn_rep = sum(Pdyn_rep_vec);

Prep = Pdyn_rep + Plk_rep;

Arep_used_mm2 = Arep_used*(1e3)^2;

%% Redo wiring now that we've changed Iidf
if(redo_wiring == 1)
    iidf_rewire = [0 iidf_rep]; % add zero-length value back in
    [Ln pn pn_orig Nm] = wire_layer_assignment_alt(iidf_rewire,lmax,Ach_wla,chi,rho_m,epsr_d,Tclk,alpha_t);
    [Cxc Ltot Cn] = calc_total_wiring_capacitance2(pn,Ln,Nm,iidf_rewire,epsr_d,gate_pitch);

    Pw = 1/2*a*Cxc*Vdd^2*f;
else
    iidf_rewire = 0;
end

%% update output objects
chip.iidf = iidf;
chip.iidf_rewire = iidf_rewire;
chip.lengths = l;

wire.pitch = pn;
wire.longest_per_tier = Ln;
wire.pitch_vec_orig = pn_orig;
wire.capacitance_total = Cxc;
wire.length_total = Ltot;
wire.capacitance_per_tier = Cn;

repeater.size_vs_minv = h_vec;
repeater.num_per_xc = k_vec;
repeater.area_total = Arep_used_mm2;
repeater.number_per_tier = num_vec;
repeater.size_per_tier = size_vec;
repeater.iidf = iidf_rep;

power.dynamic = Pdyn;
power.leakage = Plk;
power.wiring = Pw;
power.repeater = Prep;


