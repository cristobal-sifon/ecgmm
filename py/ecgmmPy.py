""" 
This is an python implementation of the error corrected Gaussian Mixture Model.
It is based on a wrapped C++ implementation of ECGMM using SWIG. 



General Input:

    -------------------------------------------------------------------
    x:     1D numpy array for the data
    xerr:  1D numpy array for the measurement errors of the data
    alpha: 1D numpy array containing initial guess of the weights of 
           each Gaussian Mixture
    mu:    1D numpy array containing initial guess of the means of 
           each Gaussian Mixture
    sigma: 1D numpy array containing initial guess of the standard 
           deviations of each Gaussian Mixture
    -------------------------------------------------------------------

Functions:
     
     ------------------------------------------------------------------
     bic_ecgmm:

         Purpose: perform the ECGMM and return Bayesian Information 
                  Criterion (BIC). The number of mixtures is determined 
                  by the input array of alpha of your initial guess. 

         Call  : bic_ecgmm(x,xerr,alpha,mu,sigma)

         Return: BIC. The input alpha, mu and sigma are also updated with
                 the final fitting values.


     ------------------------------------------------------------------
     aic_ecgmm:

         Purpose: perform the ECGMM and return Akaike Information Criterion 
                  (AIC). The number of mixtures is determined by the input 
                  array of alpha of your initial guess. 

         Call  : aic_ecgmm(x,xerr,alpha,mu,sigma)

         Return: AIC. The input alpha, mu and sigma are also updated with 
                 the final fitting values.


     -------------------------------------------------------------------
     wstat:
         
         Purpose: calculate the weighted mean and standard deviation

         Call: wstat(x,x_err)

         Return: (weighted mean, weighted sd)


     -------------------------------------------------------------------
     ecgmmplot: 

         Purpose: make plot of mixture of gaussians based on the fitting
                  results

         Call: ecgmmplot(x,alpha,mu,sigma)

         Return: a pylab plot object. Use pl.show() to see it


Revision History:

     3/4/2009 Created by Jiangang Hao @ Phyiscs @ Univ. of Michigan, Ann Arbor
     8/6/2009 ecgmmplot is added by Jiangang Hao @ Fermilab, Batavia

""" 

from ecgmm import *
import numpy as np
import pylab as pl



def bic_ecgmm(xx=None,xxerr=None,aalpha=None,mmu=None,ssigma=None):
  
    """
    Functions:
     
     ------------------------------------------------------------------
     bic_ecgmm:

         Purpose: perform the ECGMM and return Bayesian Information 
                  Criterion (BIC). The number of mixtures is determined 
                  by the input array of alpha of your initial guess. 

         Call  : bic_ecgmm(x,xerr,alpha,mu,sigma)

         Return: BIC. The input alpha, mu and sigma are also updated with
                 the final fitting values.
    """
    if xxerr == None:
        xxerr = np.zeros(len(xx))
    M=len(xx)
    N=len(aalpha)
    x=DoubleVector(M)
    xerr=DoubleVector(M)
    alpha=DoubleVector(N)
    mu=DoubleVector(N)
    sigma=DoubleVector(N)

    for i in range(0,M):
        x[i]=np.double(xx[i])
        xerr[i]=np.double(xxerr[i])
    
    for i in range(0,N):
        alpha[i]=np.double(aalpha[i])
        mu[i]=np.double(mmu[i])
        sigma[i]=np.double(ssigma[i])

    BIC=BICecgmm(x,xerr,alpha,mu,sigma)
      
    for i in range(0,N):
        aalpha[i]=alpha[i]
        mmu[i]=mu[i]
        ssigma[i]=sigma[i]
    return(BIC)




def aic_ecgmm(xx=None,xxerr=None,aalpha=None,mmu=None,ssigma=None):

    """ 
    aic_ecgmm:

         Purpose: perform the ECGMM and return Akaike Information Criterion 
                  (AIC). The number of mixtures is determined by the input 
                  array of alpha of your initial guess. 

         Call  : aic_ecgmm(x,xerr,alpha,mu,sigma)

         Return: AIC. The input alpha, mu and sigma are also updated with 
                 the final fitting values.
    """
    if xxerr == None:
        xxerr = np.zeros(len(xx))    
    M=len(xx)
    N=len(aalpha)
    x=DoubleVector(M)
    xerr=DoubleVector(M)
    alpha=DoubleVector(N)
    mu=DoubleVector(N)
    sigma=DoubleVector(N)

    for i in range(0,M):
        x[i]=np.double(xx[i])
        xerr[i]=np.double(xxerr[i])
    
    for i in range(0,N):
        alpha[i]=np.double(aalpha[i])
        mu[i]=np.double(mmu[i])
        sigma[i]=np.double(ssigma[i])

    AIC=AICecgmm(x,xerr,alpha,mu,sigma)
    for i in range(0,N):
        aalpha[i]=alpha[i]
        mmu[i]=mu[i]
        ssigma[i]=sigma[i]
        
    return(AIC)



  
def wstat(xx=None,xxerr=None):

    """
    wstat:
         
         Purpose: calculate the weighted mean and standard deviation

         Call: wstat(x,x_err)

         Return: (weighted mean, weighted sd, AIC, BIC)
    """
    if xxerr == None:
        xxerr = np.zeros(len(xx))
    M=len(xx)
    N=1
    x=DoubleVector(M)
    xerr=DoubleVector(M)
    alpha=DoubleVector(N)
    mu=DoubleVector(N)
    sigma=DoubleVector(N)

    for i in range(0,M):
        x[i]=np.double(xx[i])
        xerr[i]=np.double(xxerr[i])
    
    for i in range(0,N):
        alpha[i]=1.
        mu[i]=np.mean(xx)
        sigma[i]=np.std(xx)

    AIC=AICecgmm(x,xerr,alpha,mu,sigma)
    BIC=BICecgmm(x,xerr,alpha,mu,sigma)  
    return mu[0],sigma[0],AIC,BIC




def ecgmmplot(x,alpha,mu,sigma):

    """
    ecgmmplot: 

         Purpose: make plot of mixture of gaussians based on the fitting
                  results

         Call: ecgmmplot(x,alpha,mu,sigma)

         Return: a pylab plot object. Use pl.show() to see it
    """

    color=['r','g','b','c','m','y','k']
    if len(alpha)<=len(color):
        pl.hold(True)
        for i in range(0,len(alpha)):
            pl.plot(x,alpha[i]*np.exp(-(x-mu[i])**2/2./sigma[i]**2)/np.sqrt(2*3.14159265)/sigma[i],'.',color=color[i])
    else:
        print "Number of mixture exceeds 7, all will be in same color"
        pl.hold(True)
        for i in range(0,len(alpha)):
            pl.plot(x,alpha[i]*np.exp(-(x-mu[i])**2/2./sigma[i]**2)/np.sqrt(2*3.14159265)/sigma[i],'b.')

    return(0)


        




    
    
