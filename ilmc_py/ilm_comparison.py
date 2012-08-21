import math
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl

import ilm_support_functions as ilms

import time

output_data_filename = "pydata.dat"
time_start_c = time.clock()
time_start_t = time.time()

Ng = 1e9 # number of gates
N_strata = 2

r = 1e2 # strata to gate pitch ratio
w = 100

#k = 4 #rent constant
p = 0.6 # rent exponent
fo = 4 # avg fanout
alpha = fo/(fo+1) # input terminal fraction
k = 3/alpha

max_area_ratio = 0.10 # ratio of total TSV area to chip area


use_corrected = 1
[Iidf_corr, Mt_corr, Mt_intra_corr, Iexp_corr, l_corr, T_corr, t_corr, Nsc_corr, Nxc_corr, N_tsv_gates_corr, N_tsvs_1d_corr, Nc_corr, Nb_corr, N_nstc, N_nsbc] = ilms.compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w)
lmax_single_corr = 2*Nxc_corr
lmax_corr = lmax_single_corr + r


use_corrected = 0
[Iidf, Mt, Mt_intra, Iexp, l, T, t, Ns, Nx, N_tsv_gates, N_tsvs_1d, Nc, Nb, N_nst, N_nsb] = ilms.compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w)
lmax = 2*math.sqrt(Ng)

#stretch Iidf to calculate error
length_diff = len(Iidf_corr)-len(Iidf)
err_Iidf = np.append(Iidf, np.zeros(length_diff) )

err_Iidf = abs(Iidf_corr - err_Iidf)
norm_err_Iidf = err_Iidf/Iidf_corr


## Plots
pl.figure(1)
pl.clf()
pl.loglog(l,Iidf,'b-')
pl.hold(True)
pl.loglog(l_corr,Iidf_corr,'r-')
pl.xlabel('length (gates)')
pl.ylabel('I_{idf} (interconnect density)')
pl.savefig('Iidf_comparison.pdf')

pl.figure(2)
pl.clf()
pl.loglog(l_corr,err_Iidf,'b-')
pl.xlabel('length (gates)')
pl.ylabel('Raw error (gates)')
pl.savefig('Iidf_error.pdf')

pl.figure(3)
pl.clf()
pl.loglog(l_corr,norm_err_Iidf,'b-')
pl.xlabel('length (gates)')
pl.ylabel('Normalized error (gates)')
pl.savefig('Iidf_norm_error.pdf')

pl.figure(4)
pl.clf()
pl.loglog(l,Mt,'b-')
pl.hold(True)
pl.loglog(l_corr,Mt_corr,'r-')
pl.xlabel('length (gates)')
pl.ylabel('Number of gates')
pl.title('Mt comparison')
pl.savefig('Mt.pdf')

time_end_c = time.clock()
time_end_t = time.time()

runtime_c = time_end_c - time_start_c
runtime_t = time_end_t - time_start_t

print('Runtime (clock):',runtime_c)
print('Runtime (time):',runtime_t)

## Save data in comma-delimited format for direct comparison to MATLAB method
"""
uncorr = np.vstack( (l,Iidf,Mt) )
corr = np.vstack( (l_corr,Iidf_corr, Mt_corr) )
err = np.vstack( (err_Iidf, norm_err_Iidf) )

np.savetxt('uncorr_data.txt',uncorr,delimiter=',')
np.savetxt('corr_data.txt',corr,delimiter=',')
np.savetxt('err_data.txt',err,delimiter=',')
"""

# print ALL the arrays!
file_handle = open(output_data_filename, 'w')
dl = ","
af = "];"
ilms.printArrayToFile(file_handle,Iidf,dl,"py_Iidf = [",af)
ilms.printArrayToFile(file_handle,Mt,dl,"py_Mt = [",af)
ilms.printArrayToFile(file_handle,Mt_intra,dl,"py_Mt_intra = [",af)
ilms.printArrayToFile(file_handle,Iexp,dl,"py_Iexp = [",af)
ilms.printArrayToFile(file_handle,l,dl,"py_l = [",af)
"""
ilms.printArrayToFile(file_handle,T,dl,"py_T = [",af)
ilms.printArrayToFile(file_handle,t,dl,"py_t = [",af)
ilms.printArrayToFile(file_handle,Ns,dl,"py_Ns = [",af)
ilms.printArrayToFile(file_handle,Nx,dl,"py_Nx = [",af)
ilms.printArrayToFile(file_handle,N_tsv_gates,dl,"py_N_tsv_gates = [",af)
ilms.printArrayToFile(file_handle,N_tsvs_1d,dl,"py_N_tsvs_1d = [",af)
"""
ilms.printArrayToFile(file_handle,Nc,dl,"py_Nc = [",af)
ilms.printArrayToFile(file_handle,Nb,dl,"py_Nb = [",af)
ilms.printArrayToFile(file_handle,N_nstc,dl,"py_N_nstc = [",af)
ilms.printArrayToFile(file_handle,N_nsbc,dl,"py_N_nsbc = [",af)

ilms.printArrayToFile(file_handle,Iidf_corr,dl,"py_Iidf_corr = [",af)
ilms.printArrayToFile(file_handle,Mt_corr,dl,"py_Mt_corr = [",af)
ilms.printArrayToFile(file_handle,Mt_intra_corr,dl,"py_Mt_intra_corr = [",af)
ilms.printArrayToFile(file_handle,Iexp_corr,dl,"py_Iexp_corr = [",af)
ilms.printArrayToFile(file_handle,l_corr,dl,"py_l_corr = [",af)
"""
ilms.printArrayToFile(file_handle,T_corr,dl,"py_T_corr = [",af)
ilms.printArrayToFile(file_handle,t_corr,dl,"py_t_corr = [",af)
ilms.printArrayToFile(file_handle,Ns_corr,dl,"py_Ns_corr = [",af)
ilms.printArrayToFile(file_handle,Nx_corr,dl,"py_Nx_corr = [",af)
ilms.printArrayToFile(file_handle,N_tsv_gates_corr,dl,"py_N_tsv_gates_corr = [",af)
ilms.printArrayToFile(file_handle,N_tsvs_1d_corr,dl,"py_N_tsvs_1d_corr = [",af)
"""
ilms.printArrayToFile(file_handle,Nc_corr,dl,"py_Nc_corr = [",af)
ilms.printArrayToFile(file_handle,Nb_corr,dl,"py_Nb_corr = [",af)
ilms.printArrayToFile(file_handle,N_nstc,dl,"py_N_nstc = [",af)
ilms.printArrayToFile(file_handle,N_nsbc,dl,"py_N_nsbc = [",af)

file_handle.close()

## Plots (old matlab code)
"""figure(1)
clf
loglog(l,Iidf,'b-')
hold on
loglog(l,Iidf_corr,'r-')
xlabel('length (gates)')
ylabel('I_{idf} (interconnect density)')
title_string = sprintf('Interconnect density, #.1e gates, TSV area = #.2d##, #d strata',Ng,max_area_ratio*100,N_strata)
title(title_string)
legend('3D Joyner','3D Corrected')

figure(2)
clf
loglog(l,err_Iidf)
ylim([1e-2 10*max(err_Iidf)])
xlabel('length')
ylabel('Difference in I_{idf} (interconnects)')
title_string = sprintf('Raw error in I_{idf}, #.1e gates, TSV area = #.2d##, #d strata',Ng,max_area_ratio*100,N_strata)
title(title_string)

figure(3)
clf
loglog(l,Mt,'b-')
hold on
loglog(l,Mt_corr,'r-')
title_string = sprintf('M_t, #.1e gates, TSV area = #.2d##, #d strata',Ng,max_area_ratio*100,N_strata)
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
title_string = sprintf('Interconnect density, #.1e gates, TSV area = #.2d##, #d strata',Ng,max_area_ratio*100,N_strata)
title(title_string)
legend('3D Joyner','3D Corrected')

figure(5)
clf
loglog(l,norm_err_Iidf)
ylim([1e-5 1e0])
xlabel('length')
ylabel('Normalized difference in I_{idf} (interconnects)')
title_string = sprintf('Normalized error in I_{idf}, #.1e gates, TSV area = #.2d##, #d strata',Ng,max_area_ratio*100,N_strata)
title(title_string)
"""
