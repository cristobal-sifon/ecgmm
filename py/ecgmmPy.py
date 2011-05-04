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

    color=['r','b','g','c','m','y','k']
    if len(alpha)<=len(color):
        pl.hold(True)
        for i in range(0,len(alpha)):
            pl.plot(x,alpha[i]*np.exp(-(x-mu[i])**2/2./sigma[i]**2)/np.sqrt(2*3.14159265)/sigma[i],'-',color=color[i])
    else:
        print "Number of mixture exceeds 7, all will be in same color"
        pl.hold(True)
        for i in range(0,len(alpha)):
            pl.plot(x,alpha[i]*np.exp(-(x-mu[i])**2/2./sigma[i]**2)/np.sqrt(2*3.14159265)/sigma[i],'b-')

    return(0)


        
def bsecgmm(xx=None,xxerr=None,aalpha=None,mmu=None,ssigma=None,nboot=50,InfoCriteria= 'AIC'):

    """
    bsecgmm: 

         Purpose: increase the reliability of the ecgmm by bootstrapping

         Call: bsecgmm(x,alpha,mu,sigma,nboot = 50,InfoCriteria= 'AIC')
               nboot: the number of bootstrap you specify.Default is 50
               InfoCriteria: the information criteria used. Can be 'AIC' or 'BIC'. 
                             Default is 'AIC'.
         Return: The calculated information criteria as well as alpha, mu, sigma updated.
    """
    alphaB=[]
    muB=[]
    sigmaB=[]
    if InfoCriteria == 'BIC':
        BIC=np.zeros(nboot)
        for i in range(nboot):
            print i
            ok=np.random.randint(0,nboot-1,nboot)
            x=xx[ok]
            xerr=xxerr[ok]
            alpha=aalpha
            mu=mmu
            sigma=ssigma
            BIC[i]=ec.bic_ecgmm(xx=x,xxerr=xerr,aalpha=alpha,mmu=mu,ssigma=sigma)
            alphaB.append(aalpha)
            muB.append(mmu)
            sigmaB.append(ssigma)
    else:
        AIC=np.zeros(nboot)
        for i in range(nboot):
            ok=np.random.randint(0,nboot-1,nboot)
            x=xx[ok]
            xerr=xxerr[ok]
            alpha=aalpha[ok]
            mu=mmu[ok]
            sigma=ssigma[ok]
            AIC[i]=ec.aic_ecgmm(xx=x,xxerr=xerr,aalpha=alpha,mmu=mu,ssigma=sigma)
            alphaB.append(alpha)
            muB.append(mu)
            sigmaB.append(sigma)
    return(0)


def bic_ecgmmFixedComponent(xx=None,xxerr=None,aalpha=None,mmu=None,ssigma=None,fixMu=0,fixSigma=0):
  
    """
    Functions:
     
     ------------------------------------------------------------------
     bic_ecgmmFixedComponent:

         Purpose: perform the ECGMM and return Bayesian Information 
                  Criterion (BIC) with one component fixed. The index of the
                  fixed component is indicated by fixAlpha, fixeMu, fixSigma,
                  which are three integers to specify the component. The number
                  of mixtures is determined by the input array of alpha of your
                  initial guess. If you want the component unfixed, input the 
                  fixAlpha/fixMu/fixSigma as negative integer

         Call  : bic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma,0,0)
                 fix the first component (default)
                 bic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma,-1,0)
                 fix the first mu and sigma, but allow alpha change.
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

    BIC=BICecgmmFixedComponent(x,xerr,alpha,mu,sigma,fixMu,fixSigma)
      
    for i in range(0,N):
        aalpha[i]=alpha[i]
        mmu[i]=mu[i]
        ssigma[i]=sigma[i]
    return(BIC)


def aic_ecgmmFixedComponent(xx=None,xxerr=None,aalpha=None,mmu=None,ssigma=None,fixMu=0,fixSigma=0):

    """ 
    aic_ecgmmFixedComponent:

         Purpose: perform the ECGMM with one component fixed and return 
                  Akaike Information Criterion 
                  (AIC). The number of mixtures is determined by the input 
                  array of alpha of your initial guess. 

         Call  : aic_ecgmmFixedComponent(x,xerr,alpha,mu,sigma,0,0)
                 see the bic_ecgmmFixedComponent part for more details.
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

    AIC=AICecgmmFixedComponent(x,xerr,alpha,mu,sigma,fixMu=0,fixSigma=0)
    for i in range(0,N):
        aalpha[i]=alpha[i]
        mmu[i]=mu[i]
        ssigma[i]=sigma[i]
        
    return(AIC)




    
    
