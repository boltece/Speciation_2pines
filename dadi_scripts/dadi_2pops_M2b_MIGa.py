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


# ## Ongoing asymmetrical migration, M2b (MIGa)



pts_1 = [200,220,240]

def asym_mig(params, ns, pts):
	nuP, nuR, T1, mPR, mRP = params
	"""
	2-populations, ongoing gene flow asymmetrically

	nuP:  Current size of PUNG population, after split.
	nuR:  Current size of RIG population.
	T1:   Divergence time for split between rigida and pungens species
	mPR:  Migration from rigida into pungens
    mRP:  Migration from pungens into rigida
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP, nu2=nuR, m12=mPR, m21=mRP)

	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = asym_mig


## Setting parameter bounds and starting values
upper_bound = [100, 100, 10, 100, 100] 
lower_bound = [1e-3,1e-3,1e-3, 0.1, 0.1]
p0 = array([0.05,0.05,0.015, 40, 40])



## Run BFGS optimizer
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=25)



print('Best-fit parameters: {0}'.format(popt))



## Find the SFS based on the model and the optmized paramters
model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))

theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))



plt=dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=2, pop_ids =("pungens","rigida"))
matplotlib.pyplot.savefig("M2b_asymMIG.pdf")
