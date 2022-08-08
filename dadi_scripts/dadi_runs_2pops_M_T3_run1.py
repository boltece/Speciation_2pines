
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



data = dadi.Spectrum.from_file("pungens-rigida.sfs")


# In[3]:


ns=data.sample_sizes


# ## Gene flow, symmetrical migration and pop size change across three epochs

# ### M_T3, Run 1

# In[4]:


pts_1 = [220,240,260]

def MigChange(params, ns, pts):
	nuP1, nuR1, T1, mS1, nuP2, nuR2, T2, mS2, nuP3, nuR3, T3, mS3 = params
	"""
	2-populations, secondary contact model with asymmetrical gene flow.
    
	nuP1:  effective size of PUNG population, after split.
	nuR1:  effective size of RIG population, after split.
	T1:   First epoch. Divergence between rigida and pungens species
	nuP2:  effective size of PUNG population, after pop size change.
	nuR2:  effective size of RIG population, after pop size change.
	T2:   Second epoch.
    nuP3:  effective size of PUNG population, third epoch.
    nuR3:  effective size of RIG population, third epoch.
    T3:   Third epoch, most recent.
	mS1:  migration during T1
	mS2:  migration during T2
    mS3:  migration during T3
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuP1, nu2=nuR1, m12=mS1, m21=mS1)
	phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuP2, nu2=nuR2, m12=mS2, m21=mS2)
	phi = dadi.Integration.two_pops(phi,xx,T3,nu1=nuP3, nu2=nuR3, m12=mS3, m21=mS3)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = MigChange




## Setting parameter bounds and starting values
upper_bound = [50, 50, 5, 50, 1, 1, 1, 50, 50, 50, 1, 50]
lower_bound = [1, 1, 0.1, 1, 1e-3, 1e-3, 1e-3,  1, 1e-3, 1e-3, 1e-3, 0.1]
p0 = array([20, 20, 1, 20, 0.1, 0.1, 0.1,  20, 5, 5, 0.01, 20])


# In[6]:


p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=25)


# In[7]:


print('Best-fit parameters: {0}'.format(popt))


# In[8]:


model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))

theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))


# In[9]:


plt=dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=2, pop_ids =("pungens","rigida"))

matplotlib.pyplot.savefig("M_T3_MigPopChange_run1.pdf")

