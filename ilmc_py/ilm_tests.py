import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pylab

import ilm_support_functions as ilms

import time


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

Nt = Ng # total number of gates
S = N_strata # number of strata
Ns = Nt/S # Avg number of gates per stratum
Nx = int(math.floor(math.sqrt(Ns)))

lmax_2d = 2*math.sqrt(Ns) # max 2D wirelength (assuming optimal routing)
lmax_3d = lmax_2d + (N_strata-1)*r # max 3D wirelength (assuming optimal routing)
lmax = int(lmax_3d)

## get corrected mt_intra
N_tsv_gates = Ns* max_area_ratio/(1-max_area_ratio)

N_tsv_gates_1d = math.floor(math.sqrt(N_tsv_gates))
N_tsvs_1d = math.floor(N_tsv_gates_1d/w)

# corrected gate counts (length of chip in gates, as opposed to actual
# number of gates)
Nsc = math.ceil(Ns + N_tsv_gates)
Nxc = math.floor(math.sqrt(Nsc))
Ntc = S*Nsc

T = math.floor(Nxc/N_tsvs_1d)
t = 1/2*(T-w)


# test h
lxmax = int(Nxc)
lengths_x = range(lxmax+1)
h = np.zeros(len(lengths_x))
g = np.zeros(len(lengths_x))

time_hs = time.clock()
for lx in lengths_x:
	#g[lx] = ilms.calc_g(lx,t,w,T,N_tsvs_1d,lxmax)
	h[lx] = ilms.calc_h(lx,t,w,T,N_tsvs_1d,lxmax)

time_he = time.clock()
time_ht = time_he-time_hs

print("calc_h ran ",len(lengths_x)," times in ", time_ht, " seconds")

# test g
time_gs = time.clock()
for lx in lengths_x:
	g[lx] = ilms.calc_g(lx,t,w,T,N_tsvs_1d,lxmax)
	#h[lx] = ilms.calc_h(lx,t,w,T,N_tsvs_1d,lxmax)

time_ge = time.clock()
time_gt = time_he-time_hs
print("calc_g ran ",len(lengths_x)," times in ", time_gt, " seconds")

# now for the longer part, the actual convolution
## Calculate!
time_huge_s = time.clock()

lengths = range(lmax+1)
j2d = np.zeros(len(lengths))
f3d = np.zeros(len(lengths))
term2 = np.zeros(len(lengths))
term3 = np.zeros(len(lengths))

for l in lengths:
	j2d[l] = ilms.joyner_lin(l,Nt)
	
	if(l < Nx):
		lxmax_inner = l
		lxmin_inner = 0
	else:
		lxmax_inner = lxmax
		lxmin_inner = l-Nx

	lengths_x_inner = range(lxmin_inner,lxmax_inner+1)
	
	for lx in lengths_x_inner:
		ly = l-lx

		gx = g[lx]
		hx = h[lx]
		
		if ly not in lengths_x_inner:  # using lengths_x_inner because we assume symmetric chip   #(ly_ind > lxmax_inner +1) or (ly_ind < lxmin_inner + 1):
			gy = 0
			hy = 0
		else:
			gy = g[ly]
			hy = h[ly]
		
		term2[l] = term2[l] + gx*gy
		term3[l] = term3[l] + hx*hy
		#term3(l_ind) = term3(l_ind) +
		#calc_h(lx,t,w,T,Ntsvs,lxmax)*calc_h(l-lx,t,w,T,Ntsvs,lxmax)
	
	if l == 0:
		j2dnum = j2d[l]
	else:
		j2dnum = j2d[l]/2

	
	f3d[l] = j2dnum - 1/Nt**2 *(term2[l] + term3[l])

Mt = f3d*Nt**2

time_huge_e = time.clock()
time_huge_t = time_huge_e - time_huge_s

print("slow convolution ran ",len(lengths)," times in ", time_huge_t, " seconds")
print("Full calculation would take ", lmax/len(lengths)*time_huge_t, " seconds")
print(lmax)


## numpy convolution
time_hconv_s = time.clock()
hconv = np.convolve(h,h,'same')
time_hconv_e = time.clock()
time_hconv_t = time_hconv_e - time_hconv_s

print("numpy convolution ran fully in ", time_hconv_t, " seconds")
