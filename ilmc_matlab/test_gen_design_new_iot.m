%clear all
%close all
%% Model flags
use_joyner = 0; % 0=Use corrected distribution, 1=use original joyner distribution
redo_wiring = 0; % Redo wire layer assignment after repeater insertion? 0=no, 1=yes

%% Chip descriptors

% Stack parameters
Ng = 1.0e9;
Ach_mm2 = 100;
Ach_m2 = Ach_mm2*1e-6;

% gate parameters
eps_ox = 25; % HfO2
tox = 1e-9;
w_trans = 32e-9;
N_trans_per_gate = 4;
Ioff = 10e-9; %(A/um)

% Tsv parameters
Atf_max = 0.011; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 90e-9;
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

%% Pack objects with inputs
chip.num_gates = Ng;
chip.alpha = alpha;
chip.rent_k = k;
chip.rent_p = p;
chip.num_layers = 1;
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

gate.output_resistance = Ro;
gate.num_transistors = N_trans_per_gate;
%gate.pitch = gate_pitch;

tsv.aspect_ratio = AR_tsv;
tsv.max_area_fraction = Atf_max;
tsv.height = h_tsv_m_thin;

wire.delay_constant = alpha_t;
wire.resistivity = rho_m;
wire.dielectric_epsr = epsr_d;
wire.layers_per_tier = 2;
wire.routing_efficiency = 0.4;
wire.Beta = [0.25 0.5 0.9];
wire.Rc = 0;

simulation.use_joyner = use_joyner;
simulation.redo_wiring_after_repeaters = redo_wiring;


%% Corrected distribution
h_tsv_m = h_tsv_m_thin;
tsv.height = h_tsv_m_thin;

S = 2;
[ iidf_3d2c l_3d2c Ln_3d2c pn_3d2c pn_orig_3d2c Cxc_3d2c Ltot_3d2c Cn_3d2c Pdyn_3d2c Plk_3d2c Pw_3d2c Prep_3d2c Ng_act_3d2c N_tsvs_3d2c ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

chip.num_layers = 2;
[chip_3d2a power_3d2a wire_3d2a repeater_3d2a] = gen_design(chip,tsv,gate,transistor,wire,simulation);


S=4;
[ iidf_3d4c l_3d4c Ln_3d4c pn_3d4c pn_orig_3d4c Cxc_3d4c Ltot_3d4c Cn_3d4c Pdyn_3d4c Plk_3d4c Pw_3d4c Prep_3d4c Ng_act_3d4c N_tsvs_3d4c ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

chip.num_layers = 4;
[chip_3d4a power_3d4a wire_3d4a repeater_3d4a] = gen_design(chip,tsv,gate,transistor,wire,simulation);


chip.num_layers = 1;
[chip_2da power_2da wire_2da repeater_2da] = gen_design(chip,tsv,gate,transistor,wire,simulation);

S=1;
[ iidf_2dc l_2dc Ln_2dc pn_2dc pn_orig_2dc Cxc_2dc Ltot_2dc Cn_2dc Pdyn_2dc Plk_2dc Pw_2dc Prep_2dc Ng_act_2dc N_tsvs_2dc ] = gen_design_old(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

% %% Joyner distribution
% h_tsv_m = h_tsv_m_thick;
% tsv.height = h_tsv_m_thick;
% 
% chip.num_layers = 2;
% [chip_3d2b power_3d2b wire_3d2b repeater_3d2b] = gen_design(chip,tsv,gate,transistor,wire);
% 
% chip.num_layers = 4;
% [chip_3d4b power_3d4b wire_3d4b repeater_3d4b] = gen_design(chip,tsv,gate,transistor,wire);
% 
% chip.num_layers = 1;
% [chip_2db power_2db wire_2db repeater_2db] = gen_design(chip,tsv,gate,transistor,wire);
% 
% 
% 
% %% WLD comparison
% lw = 2;
% 
% figure(1)
% clf
% loglog(iidf_2dj,'k','linewidth',lw);
% hold on
% %loglog(iidf_2dc,'k--','linewidth',lw,'linesmoothing','on');
% loglog(iidf_3d2j,'b','linewidth',lw);
% loglog(iidf_3d4j,'g','linewidth',lw);
% loglog(iidf_3d2c,'r','linewidth',lw);
% loglog(iidf_3d4c,'m','linewidth',lw);
% xlim([1e4 1e5])
% ylim([1e-1 1e2])
% 
% legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um')
% ylabel('Number of interconnects')
% xlabel('length (gate pitches)')
% %title('Iidf')
% 
% %% Power comparison
% figure(2)
% clf
% Pw2d = [Pw_2dj Pw_2dj];
% Pw3d2 = [Pw_3d2j Pw_3d2c];
% Pw3d4 = [Pw_3d4j Pw_3d4c];
%     
% h = [300 10];
% plot(h,Pw2d,'k-o')
% hold on
% plot(h,Pw3d2,'k--s')
% plot(h,Pw3d4,'k-.d')
% legend('2D','2 tier','4 tier')
% title('Wiring power')
% %set(gca,'xscale','log')
% %ylim([0 1.3*Pw_2dj])
% xlim([10 300])
% xlabel('TSV height')
% 
% %% Pitch comparison
% 
% figure(2)
% clf
% plot(pn_2dj,'k','linewidth',lw)
% hold on
% plot(pn_3d2j,'b','linewidth',lw)
% plot(pn_3d4j,'g','linewidth',lw)
% plot(pn_3d2c,'r','linewidth',lw)
% plot(pn_3d4c,'m','linewidth',lw)
% 
% plot(pn_orig_2dj,'k--','linewidth',lw)
% plot(pn_orig_3d2j,'b--','linewidth',lw)
% plot(pn_orig_3d4j,'g--','linewidth',lw)
% plot(pn_orig_3d2c,'r--','linewidth',lw)
% plot(pn_orig_3d4c,'m--','linewidth',lw)
% 
% legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um')
% xlabel('wiring tier')
% ylabel('Wire pitch (gate pitches)')
% %title('Wire pitch')
% 
% %% Longest wire routed
% figure(3)
% lw = 2;
% clf
% semilogy(Ln_2dj,'k','linewidth',lw)
% hold on
% semilogy(Ln_3d2j,'b','linewidth',lw)
% semilogy(Ln_3d4j,'g','linewidth',lw)
% semilogy(Ln_3d2c,'r','linewidth',lw)
% semilogy(Ln_3d4c,'m','linewidth',lw)
% 
% legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um','location','se')
% 
% ylabel('longest wire routed (gate pitches)')
% xlabel('Wiring tier')
% title('Longest wire routed by metal layer')
% 
% %% Capacitance per layer
% figure(4)
% clf
% plot(Cn_2dj,'k')
% hold on
% plot(Cn_3d2j,'b')
% plot(Cn_3d4j,'g')
% plot(Cn_3d2c,'r')
% plot(Cn_3d4c,'m')
% title('Capacitance per layer')
% 
% 
% %% Total capacitance
% figure(5)
% clf
% % bar([Cxc_2dj Cxc_3d2j Cxc_3d4j Cxc_3d2c Cxc_3d4c])
% % set(gca,'xticklabel',{'2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um'})
% C2d = [Cxc_2dj Cxc_2dj];
% C3d2 = [Cxc_3d2j Cxc_3d2c];
% C3d4 = [Cxc_3d4j Cxc_3d4c];
% h = [300 10];
% plot(h,C2d,'k-o')
% hold on
% plot(h,C3d2,'k--s')
% plot(h,C3d4,'k-.d')
% legend('2D','2 tier','4 tier')
% title('Total wiring capacitance')
% %set(gca,'xscale','log')
% ylim([0 1.3*Cxc_2dj])
% xlim([10 300])
% xlabel('TSV height')
% 
% %% Another power figure
% 
% 
% P2d = [Pw_2dj Prep_2dj];
% P3d2 = [Pw_3d2c Prep_3d2c];
% P3d4 = [Pw_3d4c Prep_3d4c];
% gap = [0 0];
% 
% Pplot = [P2d ; gap ;  P3d2 ; gap ; P3d4];
% Pplot = [P2d;  P3d2 ; P3d4];
% 
% figure(6)
% clf
% barh = bar(Pplot,1.0,'stack');
% set(barh,{'FaceColor'},{'b';'r'})
% %set(gca,'xtick',[1 3 5])
% set(gca,'xticklabel',{'1','2','4'})
% legend('Wiring','Repeaters','location','ne')
% xlabel('# layers')
% ylabel('Power (W)')
% %%
% fixfigs(1:5,2.5,12,12)

