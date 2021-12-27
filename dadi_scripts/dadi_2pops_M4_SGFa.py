#!/usr/bin/env python
# coding: utf-8


import matplotlib
matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration


data = dadi.Spectrum.from_file("pungens-rigida.sfs")

ns=data.sample_sizes


# # Speciation with assymetrical gene flow, M4 (SGFa)

# In[4]:


## Set grid space to explore the parameters
pts_1 = [200,220,240]

def pureSGF(params, ns, pts):
	nuP, nuR, mPR, mRP, T1,T2 = params
	"""
	2-populations, 'pure' ancient gene flow (speciation with gene flow) model.

	
	nuP:  Current size of PUNG population, after split.
	nuR:  Current size of RIG population.
	mPR:  Migration from PUNG to RIG
	mRP:  Symmetric migration parameter, for ancestral migration during split
	T1:   Divergence time for split between rigida and pungens species
	T2:	  Divergence time duration for no migration 
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP,nu2=nuR, m12=mPR, m21=mRP)
	phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuP, nu2=nuR, m12=0, m21=0)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = pureSGF


## Setting parameter bounds and starting values
upper_bound = [ 100, 100, 100, 100, 10, 10] 
lower_bound = [1e-5,1e-5, 0.05, 0.05, 1e-6, 1e-6]
p0 = array([0.5,0.5, 40, 40, 0.015, 0.015])


## Run BFGS optimizer
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=25)



print('Best-fit parameters: {0}'.format(popt))



model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))

theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))


## plotting
plt= dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range = 2, pop_ids = ("pungens","rigida"))
matplotlib.pyplot.savefig("M4_SGFasym.pdf")
