import math
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

## Support functions for corrected/joyner 3D interconnect length model
## William Wahby, May 25 2012

# kroenecker delta selector. Returns 1 if argument is 0, 0 else
def dd(b):
	a = np.zeros(len(b))
	for i in range(len(b)):
		if b[i] == 0:
			a[i] = 1
		else:
			a[i] = 0
	return a

# boolean kroenecker delta selector. Returns True if 0, False else
def ddb(b):
	a = (b == 0 )
	return a

def ddf(b):
	if (b == 0):
		a = 1
	else:
		a = 0
	return a


def calc_g(lx,t,w,T,Ntsvs,lxmax):
# lx - 1d tsv length
# t - TSV offset within cell (in gate lengths)
# w - TSV width (in gate lengths)
# T - cell period (in gate lengths)
# Ntsvs - number of tsvs in 1d

	g = calc_h(lx,T-t-w,w,T,Ntsvs,lxmax)
	return g


def calc_h(lx,t,w,T,Ntsvs,lxmax):
# lx - 1d tsv length
# t - TSV offset within cell (in gate lengths)
# w - TSV width (in gate lengths)
# T - cell period (in gate lengths)
# Ntsvs - number of tsvs in 1d

	lmt = lx % T
	nt = Ntsvs-math.floor(lx/T)

	if (lx < 0) or (lx > lxmax):
		h = 0
	elif lmt < t:
		h = nt*w
	elif lmt > t+w:
		h = (nt-1)*w
	else:
		h = nt*w-(lmt-t)
	
	return h



def calc_Mtcorr_newcalc(Nt,Nx,t,w,T,Ntsvs_1d):

## setup
	lmax = 2*Nx
	lengths = range(lmax+1)

## translate

	Ntsvs = Ntsvs_1d
	lxmax = Nx
	lymax = Nx # symmetric
	lengths_x = range(lxmax+1)
	lengths_y = range(lymax+1)

## precalculate g and h
# we'll need these to speed up the calculation for all l
	g = np.zeros(len(lengths_x))
	h = np.zeros(len(lengths_y))

	for lx in lengths_x:
		g[lx] = calc_g(lx,t,w,T,Ntsvs,lxmax)
		h[lx] = calc_h(lx,t,w,T,Ntsvs,lxmax)

## Calculate!

	j2d = np.zeros(len(lengths))
	f3d = np.zeros(len(lengths))
	term2 = np.zeros(len(lengths))
	term3 = np.zeros(len(lengths))

	# try convolving instead of faking it with messy code
	term2 = np.convolve(g,g)
	term3 = np.convolve(h,h)

	for l in lengths:
		j2d[l] = joyner_lin(l,Nt)
		""" # commenting this out since we're convolving above
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
		"""
		if l == 0:
			j2dnum = j2d[l]
		else:
			j2dnum = j2d[l]/2

		
		f3d[l] = j2dnum - 1/Nt**2 *(term2[l] + term3[l])

	Mt = f3d*Nt**2
#Mt = f3d*Nt**2/2
#Mt(1) = Mt(1)*2

	return Mt


def compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w):
#function [Iidf Mt Iexp l T t Nsc Nxc N_tsv_gates N_tsvs_1d] = compare_calculations(use_corrected, Ng, N_strata,r,p,fo,alpha,k,max_area_ratio,w)

	Nt = Ng # total number of gates
	S = N_strata # number of strata
	Ns = Nt/S # Avg number of gates per stratum
	Nx = math.sqrt(Ns)

	lmax_2d = 2*math.sqrt(Ns) # max 2D wirelength (assuming optimal routing)
	lmax_3d = lmax_2d + (N_strata-1)*r # max 3D wirelength (assuming optimal routing)
	lmax = int(lmax_3d)

	lengths = range(lmax+1)

	# predefine some things we'll need to calculate
	Nb= np.zeros(len(lengths))
	Nc= np.zeros(len(lengths))
	Nstart= np.zeros(len(lengths))
	N_nonstart_top= np.zeros(len(lengths))
	N_nonstart_bot= np.zeros(len(lengths))
	Iidf = np.zeros(len(lengths))
	Mt = np.zeros(len(lengths))
	Iexp = np.zeros(len(lengths))


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

	if(use_corrected == 1):
		Mt_intra_corr = calc_Mtcorr_newcalc(Nsc,Nxc,t,w,T,N_tsvs_1d)
	else:
		Mt_intra_corr = 0 
		Nxc = Nx
		Nsc = Ns
		Ntc = Nt

	## calculate everything!
	for l in lengths:
		Mt[l] = getMt(Nsc,l,r,S,Mt_intra_corr,use_corrected)
		
		Nb[l] = sum(Nc[1:]) # want the sum of every Nc from l=1 to l=i-1
		
		N_nonstart_top[l] = getN_nonstart_top(Ns, l)
		N_nonstart_bot[l] = getN_nonstart_bot(Ns, l, S, r)
		Nstart[l] = Nt - N_nonstart_top[l] - N_nonstart_bot[l]
		
		Nc[l] = getNc(Nt, Ns, l, Mt[l], S, r)

		Iexp[l] = getIexp(alpha,k,p, Nb[l], Nc[l] )
		Iidf[l] = Iexp[l]*Mt[l]

	## Cleanup for returning stuff
	return (Iidf, Mt, Iexp, lengths, T, t, Nsc, Nxc, N_tsv_gates, N_tsvs_1d)



# Get number of gates in the connectible pyramid
def getAltNg_pyr(r,r_step,N_step):
# using this function because Joyner has a typo somewhere 
	Ng_pyr = 0

	for v_step in range(N_step):
		a = r_step-v_step*r
		Ng_step = 2*a*(a+1)
		Ng_pyr = Ng_pyr + Ng_step
	
	return Ng_pyr



# Calculate probability of interconnection
def getIexp(alpha,k,p,Nb,Nc):

	iexp = alpha*k/Nc*( (1+Nb)**p + (Nb + Nc)**p - (1+Nb+Nc)**p - Nb**p)
	return iexp


# Calculate the interconnect distribution function. This makes the calls that get everything rolling
def getIidf(Nt,lmax,S,r,alpha,k,p):

	Ns = Nt/S

	l = range(0,lmax+1)
	Nb= zeros(1,len(l))
	Nc= zeros(1,len(l))
	Iidf = zeros(1,len(l))
	Mt = zeros(1,len(l))
	Iexp = zeros(1,len(l))

	for i in l:
		Mt[i] = getMt(Ns,l[i],r,S)
		
		Nb[i] = sum(Nc[1:]) # want the sum of every Nc from l=1 to l=i-1
		
		Nc[i] = getNc(Nt, Ns, l[i], Mt[i], S, r)

		Iexp[i] = getIexp(alpha,k,p, Nb[i], Nc[i] )
		Iidf[i] = Iexp[i]*Mt[i]
		
		## delete stuff below this
		Nstart[i] = getNstart(Nt, Ns, l[i], S, r)

	return (Iidf, Iexp, Mt, Nb, Nc, Nstart)


def getMt(Ns,l,r,S,Mt_intra_corr,use_corrected):

	Mt = 0

	for v in range(S):
		a = S-v
		b = 2-ddf(l-v*r)-ddf(v)
		c = getMt_int_ra(Ns,l-v*r,Mt_intra_corr,use_corrected)
		Mt = Mt + a*b*c

		return Mt


def getMt_int_ra(Ns,l,Mt_intra_corr,use_corrected):

	if (use_corrected == 0):
		if (l == 0):
			Mt_int_ra = Ns
		elif ( (0 < l) and (l < math.sqrt(Ns) ) ):
			Mt_int_ra = ( 2*Ns*l - 2*math.sqrt(Ns)*(l**2) + 1/3*l**3 )
		elif ( (math.sqrt(Ns) <= l) and (l < 2*math.sqrt(Ns)-1) ):
			Mt_int_ra = 1/3 * (2*math.sqrt(Ns) - l)**3
		else:
			Mt_int_ra = 0
	else:
		if (l < 0):
			Mt_int_ra = 0
		elif (l >= 2*math.sqrt(Ns) ):
			Mt_int_ra = 0
		else:
			Mt_int_ra = Mt_intra_corr[l]

	return Mt_int_ra



def getN_nonstart_bot(Ns, l, S, r):

	A = getAltNg_pyr(r, max(0,l-math.sqrt(Ns)-r), min(S-1, math.floor( (l-math.sqrt(Ns))/r) ) )
	B = -2*getAltNg_pyr(r, max(0, l-3/2*math.sqrt(Ns)-r), min(S-1, math.floor( (l-3/2*math.sqrt(Ns))/r) ) )

	N_nonstart_bot = A+B

	return N_nonstart_bot


def getN_nonstart_top(Ns, l):

	if ((0 <= l) and (l<= math.sqrt(Ns)/2)):
		N_nonstart_top = l
	elif ( (math.sqrt(Ns)/2 < l) and (l <= math.sqrt(Ns) ) ):
		N_nonstart_top = l + (l -math.sqrt(Ns)/2 -1)*(l -math.sqrt(Ns)/2)
	elif ( (math.sqrt(Ns) < l) and (l <= 3*math.sqrt(Ns)/2) ):
		N_nonstart_top = l*math.sqrt(Ns) - 3*Ns/4 + math.sqrt(Ns)/2
	elif (3*math.sqrt(Ns)/2 < l) and (l <= 2*math.sqrt(Ns) ): # Think J has a typo here. 2math.sqrt(Ns)?
		N_nonstart_top = Ns - (2*math.sqrt(Ns)-l)*(2*math.sqrt(Ns)-l-1)
	else:
		N_nonstart_top = 0
	
	return N_nonstart_top



def getNb(Nc):

	Nb = np.zeros(1,len(Nc))

	for l in range(Nc):
		Nb[l] = sum(Nc[:l])

	return Nb


			
def getNc(Nt, Ns, l, Mt, S, r):

	Nstart = getNstart(Nt, Ns, l, S, r)
	Nc = Mt/Nstart

	return Nc


def getNg_pyr(r,r_step,N_step):

	A = r**2/6*(N_step-1)*N_step*(2*N_step-1)
	B = -r/2*(2*r_step+1)*(N_step-1)*N_step
	C = N_step*(r_step**2 + r_step)
	Ng_pyr = 2*(A + B + C)

	return Ng_pyr


def getNstart(Nt, Ns, l, S, r):
	N_nonstart_top = getN_nonstart_top(Ns, l)
	N_nonstart_bot = getN_nonstart_bot(Ns, l, S, r)
	Nstart = Nt - N_nonstart_top - N_nonstart_bot

	return Nstart


def joyner_lin(l,Nt):
	if (l < 0):
		a = 0
	elif (l==0):
		a = 1/Nt
	elif (l < math.sqrt(Nt) ):
		a = 1/Nt**2 * (4*Nt*l - 4*math.sqrt(Nt)*l**2 + 2/3*l**3)
	elif ( l < 2*math.sqrt(Nt) - 1):
		a = 2/3/Nt**2 * (2*math.sqrt(Nt) - l)**3
	else:
		a = 0

	return a



def rect(b):

	a = np.zeros(1,len(b))
	for i in range(len(b)):
		if ( abs(b) <= 0.5):
			a[i] = 1

	return a




