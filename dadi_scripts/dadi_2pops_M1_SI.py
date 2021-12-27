#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration


# In[2]:


data = dadi.Spectrum.from_file("pungens-rigida.sfs")


# In[3]:


ns=data.sample_sizes


# # Strict Isolation model (M1)

# In[4]:


pts_1 = [200,220,240]

def SI(params, ns, pts):
	nuP, nuR, T1 = params
	"""
	3-populations, strict isolation model with no gene flow.

	nuP:  Current size of PUNG population, after split.
	nuR:  Current size of RIG population.
	T1:   Divergence time for split between rigida and pungens species
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP, nu2=nuR)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = SI




#### changed upper and lower bounds to reflect first run being closer to lower bound 
upper_bound = [0.1, 0.1, 0.05] 
lower_bound = [1e-7,1e-7, 1e-7]
p0 = array([0.005,0.005, 0.001])




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



plt=dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=2, pop_ids =("pungens","rigida"))

matplotlib.pyplot.savefig("M1_SI.pdf")

