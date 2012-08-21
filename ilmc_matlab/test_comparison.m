% test stuff
clear all
close all
clc
 
tic % start timing!
Ng = 1e9; % number of gates
N_strata = 2;

r = 1e2; % strata to gate pitch ratio
w = 100;

%k = 4; %rent constant
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha;

max_area_ratio = 0.10; % ratio of total TSV area to chip area


use_corrected = 1;
[Iidf_corr Mt_corr Iexp_corr l_corr T_corr t_corr Nsc_corr Nxc_corr N_tsv_gates_corr N_tsvs_1d_corr Nc_corr Nb_corr N_nstc N_nsbc] = compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w);
lmax_single_corr = 2*Nxc_corr;
lmax_corr = lmax_single_corr + r;


use_corrected = 0;
[Iidf Mt Iexp l T t Ns Nx N_tsv_gates N_tsvs_1d Nc Nb N_nst N_nsb ] = compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w);
lmax = 2*sqrt(Ng);

err_Iidf = abs(Iidf_corr - Iidf);
norm_err_Iidf = err_Iidf./Iidf_corr;

toc % stop timing!

%% Plots
figure(1)
clf
loglog(l,Iidf,'b-')
hold on
loglog(l,Iidf_corr,'r-')
xlabel('length (gates)')
ylabel('I_{idf} (interconnect density)')
title_string = sprintf('Interconnect density, %.1e gates, TSV area = %.2d%%, %d strata',Ng,max_area_ratio*100,N_strata);
title(title_string)
legend('3D Joyner','3D Corrected')

figure(2)
clf
loglog(l,err_Iidf)
ylim([1e-2 10*max(err_Iidf)])
xlabel('length')
ylabel('Difference in I_{idf} (interconnects)')
title_string = sprintf('Raw error in I_{idf}, %.1e gates, TSV area = %.2d%%, %d strata',Ng,max_area_ratio*100,N_strata);
title(title_string)

figure(3)
clf
loglog(l,Mt,'b-')
hold on
loglog(l,Mt_corr,'r-')
title_string = sprintf('M_t, %.1e gates, TSV area = %.2d%%, %d strata',Ng,max_area_ratio*100,N_strata);
title(title_string)
legend('Location','South','3D Joyner','3D Corrected')

figure(4)
clf
loglog(l,Iidf,'b-')
hold on
loglog(l,Iidf_corr,'r-')
axis([1 l(end) 1 10*max(Iidf)])
xlabel('length (gates)')
ylabel('I_{idf} (interconnect density)')
title_string = sprintf('Interconnect density, %.1e gates, TSV area = %.2d%%, %d strata',Ng,max_area_ratio*100,N_strata);
title(title_string)
legend('3D Joyner','3D Corrected')

figure(5)
clf
loglog(l,norm_err_Iidf)
ylim([1e-5 1e0])
xlabel('length')
ylabel('Normalized difference in I_{idf} (interconnects)')
title_string = sprintf('Normalized error in I_{idf}, %.1e gates, TSV area = %.2d%%, %d strata',Ng,max_area_ratio*100,N_strata);
title(title_string)

%% Load python data and compare
pydata2

dl_corr = length(py_Iidf_corr) - length(Iidf_corr);
dl = length(py_Iidf) - length(Iidf);
% extend matlab vectors to make direct comparisons
Iidf_corr2 = [Iidf_corr zeros(1,dl_corr)];
Iidf2 = [Iidf zeros(1,dl)];

diff_Iidf_corr = Iidf_corr2 - py_Iidf_corr;
diff_Iidf = Iidf2 - py_Iidf;

figure(6)
clf
semilogy(l_corr,N_nstc,'b')
hold on
semilogy(py_l_corr,py_N_nstc,'r--')
xlabel('interconnect length (gates)')
ylabel('Nonstarting top gates')
title('Nonstarting top gates - Matlab vs python')
legend('MATLAB','Python','Location','Southeast')

figure(7)
clf
semilogy(l_corr,N_nsbc,'b')
hold on
semilogy(py_l_corr,py_N_nsbc,'r--')
xlabel('interconnect length (gates)')
ylabel('Nonstarting bottom gates')
title('Nonstarting bottom gates - Matlab vs python')
legend('MATLAB','Python','Location','Southeast')


figure(8)
clf
semilogy(py_l_corr,abs(diff_Iidf_corr),'b')
hold on
semilogy(py_l,abs(diff_Iidf),'r')
legend('corr','uncorr')
xlabel('length (gates)')
ylabel('Raw error (IDF)')
title('Error comparison - Matlab vs Python')
