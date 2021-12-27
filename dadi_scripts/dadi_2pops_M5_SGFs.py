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


# # Speciation with symmetrical gene flow model, M5 (SGFs)


## Set grid space to explore the parameters
pts_1 = [200,220,240]

def pureSGFsym(params, ns, pts):
	nuP, nuR, T1, T2, mA = params
	"""
	2-populations, 'pure' ancient gene flow (speciation with gene flow) model.

	
	nuP:  Current size of PUNG population, after split.
	nuR:  Current size of RIG population.
	mA:  Symmetrical Migration 
	T1:   Divergence time for split between rigida and pungens species
	T2:	  Divergence time duration for no migration 
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP,nu2=nuR, m12=mA, m21=mA)
	phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuP, nu2=nuR, m12=0, m21=0)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = pureSGFsym


## Setting parameter bounds and starting values
upper_bound = [10, 10, 10, 10, 100] 
lower_bound = [1e-5,1e-5,1e-5,1e-5, 0.25]
p0 = array([0.5,0.5,0.15,0.15, 40])


# Run BFGS optimizer
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt01 = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=25)


print('Best-fit parameters: {0}'.format(popt01))


model01 = func_ex(popt01, ns, pts_1)
ll_model01 = dadi.Inference.ll_multinom(model01, data)
print('Maximum log composite likelihood: {0}'.format(ll_model01))

theta01 = dadi.Inference.optimal_sfs_scaling(model01, data)
print('Optimal value of theta: {0}'.format(theta01))


## plotting
plt= dadi.Plotting.plot_2d_comp_multinom(model01, data, vmin=1, resid_range = 2, pop_ids = ("pungens","rigida"))
matplotlib.pyplot.savefig("M5_SGFsym.pdf")

