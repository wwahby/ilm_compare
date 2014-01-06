%clear all
%close all
%% Model flags
use_joyner = 0; % 0=Use corrected distribution, 1=use original joyner distribution
redo_wiring = 0; % Redo wire layer assignment after repeater insertion? 0=no, 1=yes

%% Chip descriptors

% Stack parameters
Ng = 2.5e8;
S = 4;
Ach_mm2 = 100;
Ach_m2 = Ach_mm2*1e-6;

% gate parameters
eps_ox = 25; % HfO2
tox = 1e-9;
w_trans = 30e-9;
N_trans_per_gate = 4;
Ioff = 10e-9; %(A/um)

% Tsv parameters
Atf_max = 0.10; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 300e-9;
h_tsv_m_thin = 10e-6;
h_tsv_m_thick = 300e-6;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant

% Wiring parameters
fmax = 1e9;
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 3.0; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

% Repeater parameters
Ro = 1e3; % (Ohm) Gate output resistance

% Power parameters
a = 0.1; % logic activity factor
Vdd = 1.0; % (V)

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity

%% Pack objects with inputs
chip.num_gates = Ng;
chip.alpha = alpha;
chip.rent_k = k;
chip.rent_p = p;
chip.num_layers = S;
chip.min_pitch = 2*w_trans;
chip.gate_pitch = gate_pitch;

chip.area_total = Ach_m2;
chip.chi = chi;
chip.clock_period = Tclk;
chip.logic_activity_factor = a;
chip.Vdd = Vdd;

transistor.gate_length = w_trans;
transistor.oxide_rel_permittivity = eps_ox;
transistor.oxide_thickness = tox;
transistor.leakage_current_per_micron = Ioff;
transistor.capacitance = eps_ox*eps0*w_trans^2/tox;

gate.output_resistance = Ro;
gate.num_transistors = N_trans_per_gate;
gate.capacitance = N_trans_per_gate*transistor.capacitance;
%gate.pitch = gate_pitch;

tsv.aspect_ratio = AR_tsv;
tsv.max_area_fraction = Atf_max;
tsv.height = h_tsv_m_thin;

wire.delay_constant = alpha_t;
wire.resistivity = rho_m;
wire.dielectric_epsr = epsr_d;
wire.layers_per_tier = 2;
wire.routing_efficiency = 0.4;
wire.Beta = [0.25 0.9];
wire.Rc = 0;

simulation.use_joyner = use_joyner;
simulation.redo_wiring_after_repeaters = redo_wiring;


%% Corrected distribution
h_tsv_m = h_tsv_m_thin;
tsv.height = h_tsv_m_thin;

S = 2;
[ iidf_3d2c l_3d2c Ln_3d2c pn_3d2c pn_orig_3d2c Cxc_3d2c Ltot_3d2c Cn_3d2c Pdyn_3d2c Plk_3d2c Pw_3d2c Prep_3d2c Ng_act_3d2c N_tsvs_3d2c ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

% chip.num_layers = 2;
% [chip_3d2a power_3d2a wire_3d2a repeater_3d2a] = gen_design(chip,tsv,gate,transistor,wire,simulation);
% 
% 
% S=4;
% [ iidf_3d4c l_3d4c Ln_3d4c pn_3d4c pn_orig_3d4c Cxc_3d4c Ltot_3d4c Cn_3d4c Pdyn_3d4c Plk_3d4c Pw_3d4c Prep_3d4c Ng_act_3d4c N_tsvs_3d4c ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
% 
% chip.num_layers = 4;
% [chip_3d4a power_3d4a wire_3d4a repeater_3d4a] = gen_design(chip,tsv,gate,transistor,wire,simulation);


chip.num_layers = 4;
[chip_2da power_2da wire_2da repeater_2da] = gen_design(chip,tsv,gate,transistor,wire,simulation);

S=4;
[ iidf_2dc l_2dc Ln_2dc pn_2dc pn_orig_2dc Cxc_2dc Ltot_2dc Cn_2dc Pdyn_2dc Plk_2dc Pw_2dc Prep_2dc Ng_act_2dc N_tsvs_2dc ] = gen_design_old(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

%% calc error
figure(5)
clf
plot(pn_2dc,'b')
hold on
plot(wire_2da.pn/chip.min_pitch,'r')
fixfigs(5,3,14,12)

figure(6)
clf
plot(Ln_2dc,'b')
hold on
plot(wire_2da.Ln,'r')
fixfigs(6,3,14,12)
