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


# ## Gene flow, symmetrical migration change at T2, pop size change (PSCMIGCs)


# In[4]:


pts_1 = [200,220,240]

def MigChange(params, ns, pts):
	nuP1, nuR1, T1, mS1, nuP2, nuR2, T2, mS2 = params
	"""
	2-populations, changes in pop size and symmetrical gene flow at T2.
    
	nuP1:  effective size of PUNG population, after split.
	nuR1:  effective size of RIG population, after split.
	T1:   Divergence time for split between rigida and pungens species
	nuP2:  effective size of PUNG population, after pop size change.
	nuR2:  effective size of RIG population, after pop size change.
	T2:   The time scale between secondary contact and present.
	mS1:  migration during T1
	mS2:  migration during T2
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP1, nu2=nuR1, m12=mS1, m21=mS1)
	phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuP2, nu2=nuR2, m12=mS2, m21=mS2)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = MigChange



## Setting parameter bounds and starting values
upper_bound = [100, 100, 10, 200, 100, 100, 10, 200] 
lower_bound = [1e-3, 1e-3, 1e-3, 0.1, 1e-3, 1e-3, 1e-3,  0.1]
p0 = array([1, 1, 0.1, 50, 1, 1, 0.1,  50])



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


plt=dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=2, pop_ids =("pungens","rigida"))

matplotlib.pyplot.savefig("M14_MigPopChangeT2.pdf")

