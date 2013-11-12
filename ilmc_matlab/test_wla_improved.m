%clear all
%close all
%% Chip descriptors

% Stack parameters
Ng = 2.5e8;
S = 4;
Ach_mm2 = 100;
Ach_m2 = Ach_mm2*1e-6;

% Tsv parameters
Atf_max = 0.10; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 300e-9;
h_tsv_m = 10e-6;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant


%% Presize the chip and TSVs
Ns = Ng/S;
Lx = round(sqrt(Ns));

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;

% Size the TSVs
h_tsv = ceil(h_tsv_m/gate_pitch);
w_tsv = ceil(h_tsv/AR_tsv);

%% Size the chip so we have a nicely divisible number of unit cells per side
%Nsp = floor( (1+Atf_max)*Ns );
Nsp = floor(Ns/(1-Atf_max));
Lxp = floor(sqrt(Nsp));
Nsp = Lxp^2;
Tp = ceil(w_tsv/sqrt(Atf_max));

slack = 0.2;
[Lxc Tc Nuc_1d gfrac_L gfrac_T] = find_LT_combination(Lxp,Tp,slack);

Nsc = Lxc^2;
Ngc = Nsc*S;

N_tsvs = Nuc_1d^2;
g_tsv = (Nuc_1d*w_tsv)^2; % number of gates displaced by TSVs
Atf_act = g_tsv/Nsc;

repstr1 = sprintf('Ng_nom: %.4g',Ng);
repstr2 = sprintf('Ng_act: %.4g',Ngc);
repstr3 = sprintf('Atf_act: %.4g',Atf_act);
disp(repstr1)
disp(repstr2)
disp(repstr3)

%% Calculate WLD
iidf_2d = calc_Iidf(alpha,k,p,round(sqrt(Ng)),1,h_tsv);
iidf_c = calc_Iidf_corrected(alpha,k,p,Lxc,S,h_tsv,Nuc_1d,w_tsv);
iidf_j = calc_Iidf(alpha,k,p,Lx,S,h_tsv);

%AA = calc_Nstart(Lx,S,r);
% Mt3dj = Mt_3d_joyner(Lx,S,r);
% Mt3dc = Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
% Mt2dj = Mt_2d_joyner(Lx);
% Mt2dc = Mt2d_corrected(Lx, Nuc_1d, w_tsv);
% Nstart = calc_Nstart(Lx,S,r,g_tsv);
% Nnst = calc_Nnst(Lx,S,r,g_tsv);
% Nnsb = calc_Nnsb(Lx,S,r,g_tsv);
% Nc = calc_Nc(Mt3dc,Lx,S,r,g_tsv);
% h = calc_h(Lx, Nuc_1d, w_tsv);
% term4 = zeros(1,length(Mt2dj));
% term4(1:2*w_tsv+1) = N_tsvs*Mt_2d_joyner(w_tsv);
% iexp = calc_Iexp(alpha,k,p,Mt3dc,Lx,S,r,g_tsv);

%% Get rid of NaNs
iidf_2d(isnan(iidf_2d)) = 0;
lmax_2d = length(iidf_2d) - 1;
l2d = 0:lmax_2d;

iidf_c(isnan(iidf_c)) = 0;
lmax_c = length(iidf_c) - 1;
lc = 0:lmax_c;

iidf_j(isnan(iidf_j)) = 0;
lmax_j = length(iidf_j) - 1;
lj = 0:lmax_j;

%% Estimate total wiring capacitance

fmax = 1e9;
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 3.0; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

min_pitch = 60e-9;
%min_pitch = gate_pitch;
Ach = Ach_m2; % Force chip to use specified area 
Ach_gp = Ach_m2/gate_pitch^2;
Ach_mp = Ach_m2/min_pitch^2;
layers_per_tier = 2;
routing_efficiency = 0.4;
layer_area = Ach;


Beta = [0.25 0.9];
A_layer_max = layer_area*routing_efficiency*layers_per_tier;
Rc = 0;
%Rc = [1e2 1e1 1e0 1e-1];
rho_m_vec = [1 1 1] * rho_m;
%rho_m_vec = rho_m;

[Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved(iidf_2d,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area,rho_m_vec,epsr_d,alpha_t,Beta,Tclk,Rc);
[Ln2d pn2d pn2d_orig] = wire_layer_assignment_alt(iidf_2d,lmax_2d,Ach_gp,chi,rho_m,epsr_d,Tclk,alpha_t);
[Ln2dm pn2dm pn2d_origm] = wire_layer_assignment_alt(iidf_2d,lmax_2d,Ach_mp,chi,rho_m,epsr_d,Tclk,alpha_t);

% [Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved(iidf_c,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area/S,rho_m,epsr_d,alpha_t,Beta,Tclk,Rc);
% [Ln2d pn2d pn2d_orig] = wire_layer_assignment_alt(iidf_c,lmax_c,Ach_gp/S,chi,rho_m,epsr_d,Tclk,alpha_t);
% [Ln2dm pn2dm pn2d_origm] = wire_layer_assignment_alt(iidf_c,lmax_c,Ach_mp/S,chi,rho_m,epsr_d,Tclk,alpha_t);
%% Capacitance
[Cxc2d Ltot2d Cn2d] = calc_total_wiring_capacitance(pn2d,Ln2d,iidf_2d,l2d,epsr_d,gate_pitch);
[Cxc2dm Ltot2dm Cn2dm] = calc_total_wiring_capacitance(pn2dm,Ln2dm,iidf_2d,l2d,epsr_d,gate_pitch);
[Cxc Ltot Cno] = calc_total_wiring_capacitance(pn_vec,Ln_vec,iidf_2d,l2d,epsr_d,gate_pitch);
[Ctot Cnn] = calc_wiring_capacitance_from_area(pn_vec,layers_per_tier,A_wires,A_vias,epsr_d);

%% Size gates and see how far off we were
Vt = 0.2;
Vdd = 1;
epsr_ox = 30;
t_ox = 1e-9;
mobility_n = 1300; % (cm^2/Vs);
mobility_p = mobility_n/3;
mobility_av = 1/2*(mobility_n + mobility_p);

L_gate_min = 22e-9;
Wmin = L_gate_min;
logic_depth = 5;
Beta_g = 0.75;


Rnand = calc_Rnand(Vdd,Vt,epsr_ox,t_ox,mobility_av);
Cnand = calc_Cnand(L_gate_min,epsr_ox,t_ox);

l_av = get_average_wirelength(iidf_2d);
C_av = get_capacitance_from_length(l_av,Ln_vec,pn_vec,epsr_d,gate_pitch);
[W gate_pitch_new] = size_gates_nand2(Wmin,Rnand,Cnand,100*C_av,fo,logic_depth,Tclk,Beta_g);


%%
figure(1)
clf
plot(pn2d,'b')
hold on
plot(pn2dm,'g--')
plot(pn_vec/min_pitch,'r')
xlabel('Wiring tier')
ylabel('Wire pitch (MP)')
legend('Orig-GP','Orig-MP','Improved')
fixfigs(1,3,14,12)

figure(2)
clf
semilogy(A_wires/A_layer_max,'b');
hold on
semilogy(A_vias/A_layer_max,'r');
xlabel('Wiring tier')
ylabel('Area fraction')
legend('wires','vias')
fixfigs(2,3,14,12)

figure(3)
clf
plot(Ln2d,'b')
hold on
plot(Ln2dm,'g--')
plot(Ln_vec,'r')
xlabel('Wiring tier')
ylabel('Longest wire routed per tier')
legend('Orig-GP','Orig-MP','Improved')
fixfigs(3,3,14,12)

figure(4)
clf
plot(Cn2d.*Ltot2d*gate_pitch,'b')
hold on
plot(Cn2dm.*Ltot2dm*gate_pitch,'g')
plot(Cno.*Ltot.*gate_pitch,'m')
plot(Cnn,'r')
xlabel('Wiring tier')
ylabel('Capacitance per tier (F)')
legend('Orig-GP','Orig-MP','Improved-oldC','Improved-newC')
fixfigs(4,3,14,12)

%% Old stuff
% [Ln2d pn2d pn2d_orig] = wire_layer_assignment_alt(iidf_2d,lmax_2d,Ach,chi,rho_m,epsr_d,Tclk,alpha_t);
% [Cxc2d Ltot2d Cn2d] = calc_total_wiring_capacitance(pn2d,Ln2d,iidf_2d,l2d,epsr_d,gate_pitch);
% 
% [Lnj pnj pnj_orig] = wire_layer_assignment_alt(iidf_j,lmax_j,Ach,chi,rho_m,epsr_d,Tclk,alpha_t);
% [Cxcj Ltotj Cnj] = calc_total_wiring_capacitance(pnj,Lnj,iidf_j,lj,epsr_d,gate_pitch);
% 
% [Lnc pnc pnc_orig] = wire_layer_assignment_alt(iidf_c,lmax_c,Ach,chi,rho_m,epsr_d,Tclk,alpha_t);
% [Cxcc Ltotc Cnc] = calc_total_wiring_capacitance(pnc,Lnc,iidf_c,lc,epsr_d,gate_pitch);


% %% Calculate error
% ldiff1 = length(iidf_c) - length(iidf_j);
% ldiff2 = length(iidf_j) - length(iidf_c);
% iidf_c = [iidf_c zeros(1,ldiff2)];
% iidf_j = [iidf_j zeros(1,ldiff1)];
% 
% err_raw = abs(iidf_c-iidf_j);
% err_norm = err_raw./iidf_j;
% 
% %% WLD comparison
% figure(11)
% clf
% loglog(iidf_2d,'k');
% hold on
% loglog(iidf_j);
% loglog(iidf_c,'r--')
% title('Iidf')
% 
% %% Pitch comparison
% figure(2)
% clf
% plot(pnj,'b')
% hold on
% plot(pnj_orig,'b--')
% plot(pnc,'r')
% plot(pnc_orig,'r--')
% 
% title('Wire pitch')
% 
% %% Longest wire routed
% figure(3)
% clf
% semilogy(Lnj,'b')
% hold on
% semilogy(Lnc,'r')
% title('longest wire routed per layer')
% 
% %% Capacitance per layer
% figure(4)
% clf
% plot(Cnj,'b')
% hold on
% plot(Cnc,'r--')
% title('Capacitance per layer')
% 
% %% Total capacitance
% figure(5)
% clf
% bar([Cxc2d Cxcj Cxcc])
% title('Total wiring capacitance')
% 
