#! /usr/bin/env python

"""
This is an example about how to use the ecgmmPy module.

Created by: Jiangang Hao @ Phys @ U of Michigan 3/9/2009

"""

import numpy as np
import pylab as pl
import ecgmmPy as ec

# generate two gaussian data
ncl=80    # number of each mixture    
nbg=50

mucl=0.5  # mean of each mixture
mubg=0

sdcl=0.04 # width of each mixture
sdbg=0.3


cl=np.random.normal(mucl,sdcl,ncl)
bg=np.random.normal(mubg,sdbg,nbg)

xerr=np.random.uniform(0,0.1,ncl+nbg)

x=np.append(cl,bg)+np.random.uniform(0,1,ncl+nbg)*xerr


ec.wstat(x,xerr)   # calculate the weighted statistics

# assign initial guess as two mixtures
alpha=np.array([0.5,0.5])
mu=np.array([0.1,0.3])
sigma=np.array([0.2,0.05])

#aic=ec.aic_ecgmm(x,xerr,alpha,mu,sigma)
bic=ec.bic_ecgmm(x,xerr,alpha,mu,sigma)

#make plots
pl.figure(figsize=(12,6))

pl.subplot(1,2,1)
y=np.arange(-1,1,0.0001)
ec.ecgmmplot(y,alpha,mu,sigma)
pl.xlabel='x'
pl.text(-0.9,3, r'$\mu=$'+str(mu))
pl.text(-0.9,2.6, r'$\sigma=$'+str(sigma))
pl.title("ECGMM")



# non-error corrected GMM
xerr[:]=0
alpha=np.array([0.5,0.5])
mu=np.array([0.1,0.3])
sigma=np.array([0.2,0.05])
bic=ec.bic_ecgmm(x,xerr,alpha,mu,sigma)


#make plots
pl.subplot(1,2,2)
ec.ecgmmplot(y,alpha,mu,sigma)
pl.xlabel='x'
pl.text(-0.9,3, r'$\mu=$'+str(mu))
pl.text(-0.9,2.6, r'$\sigma=$'+str(sigma))
pl.title("GMM")



#----------fixed component ------

# assign initial guess as two mixtures
alpha=np.array([0.5,0.5])
mu=np.array([0.1,0.3])
sigma=np.array([0.2,0.05])

bic=ec.bic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma,1,-1)
# here, 1,1 means the second component of mu and sigma is fixed. 
# by default, bic=ec.bic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma), 
# it is 0,0, i.e. the first component of mu and sigma is fixed.
# if want only one sigma fixed, use bic=ec.bic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma,-1,0)
# set it to an negative integer will unfix it. 
# note that the relative weight cannot be fixed because they need to sum up to 1. 


#make plots
pl.figure(figsize=(6,6))
y=np.arange(-1,1,0.0001)
ec.ecgmmplot(y,alpha,mu,sigma)
pl.xlabel='x'
pl.text(-0.9,3.4, r'$\alpha=$'+str(alpha))
pl.text(-0.9,3, r'$\mu=$'+str(mu))
pl.text(-0.9,2.6, r'$\sigma=$'+str(sigma))
pl.title("ECGMM")
pl.show()
