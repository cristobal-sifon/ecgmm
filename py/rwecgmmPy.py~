import numpy as np
import ecgmmPy as ec

"""
this package is for the radially weighted ecgmm. 
The useful functions are:
bic1EM() and bic2EM(), each of them will return the bic as well as the estimated parameters when considering the radial weights. 
J. Hao
4/20/2011

"""
def poisson(r):
    const=2.
    res=const*r
    return res

def nfw(r):
    rhos=104.5
    rs=0.15
    x=r/rs
    if type(x) is float:
        if x < 1:
            resf = 1-2./np.sqrt(1-x**2)*np.arctanh(np.sqrt((1-x)/(x+1)))
        if x > 1:
            resf = 1-2./np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1)/(x+1)))
            res = 2*rhos*rs/(x**2-1)*resf*r
        if x == 1:
            resf = 0.
        res = 2*rhos*rs/(x**2-1)*resf*r
    else:
        resf=np.zeros(len(x))
        lt1= x < 1
        resf[lt1] = 1-2./np.sqrt(1-x[lt1]**2)*np.arctanh(np.sqrt((1-x[lt1])/(x[lt1]+1)))
        gt1= x > 1
        resf[gt1] = 1-2./np.sqrt(x[gt1]**2-1)*np.arctan(np.sqrt((x[gt1]-1)/(x[gt1]+1)))
        eq1= x == 1
        resf[eq1] = 0.
    res=2*rhos*rs/(x**2-1)*resf*r
    return res

def nfwTruncted(r):
    """
    the nfw and poisson profile all normalized at r=3Mpc.
    """
    rhos=104.5/3.0287545869954959
    rs=0.15
    x=r/rs
    if type(x) is float:
        if x < 1:
            resf = (1-2./np.sqrt(1.0001**2-1)*np.arctan(np.sqrt((1.0001-1)/(1.0001+1))))
            res = 2*rhos*rs/(1.0001**2-1)*resf
        if x > 1:
            resf = 1-2./np.sqrt(x**2-1)*np.arctan(np.sqrt((x-1)/(x+1)))
            res = 2*rhos*rs/(x**2-1)*resf*r
        if x == 1:
            resf = 0.
            res = 2*rhos*rs/(x**2-1)*resf*r
    else:
        resf=np.zeros(len(x))
        res=np.zeros(len(x))
        xtmp=np.zeros(len(x))
        lt1= x < 1
        xtmp[lt1]=x[lt1]*0+1.0001
        resf[lt1] = (1-2./np.sqrt(xtmp[lt1]**2-1)*np.arctan(np.sqrt((xtmp[lt1]-1)/(xtmp[lt1]+1))))
        res[lt1] = 2*rhos*rs/(xtmp[lt1]**2-1)*resf[lt1]*xtmp[lt1]*rs
        gt1= x > 1
        resf[gt1] = 1-2./np.sqrt(x[gt1]**2-1)*np.arctan(np.sqrt((x[gt1]-1)/(x[gt1]+1)))
        res[gt1] = 2*rhos*rs/(x[gt1]**2-1)*resf[gt1]*r[gt1]
        eq1= x == 1
        resf[eq1] = 0.
        res[eq1]=0.
    return res

def pcz1(c,delta,r,mu1,sigma1):
    res=nfw(r)*np.exp(-(c-mu1)**2/2./(sigma1**2+delta**2))/np.sqrt(2.*3.1415926*(sigma1**2+delta**2))
    return res

def pcz2(c,delta,r,mu2,sigma2):
    res=poisson(r)*np.exp(-(c-mu2)**2/2./(sigma2**2+delta**2))/np.sqrt(2.*3.1415926*(sigma2**2+delta**2))
    return res

def pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    res=pcz1(c,delta,r,mu1,sigma1)*w1/(pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2)
    return res

def pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    res=pcz2(c,delta,r,mu2,sigma2)*w2/(pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2)
    return res

def w1new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    res=sum(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2))/float(len(c))
    return res

def w2new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    res=sum(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2))/float(len(c))
    return res

def mu1new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    upper=sum(c*pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1**2/(sigma1**2+delta**2))
    lower=sum(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1**2/(sigma1**2+delta**2))
    res=upper/lower
    return res

def mu2new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    upper=sum(c*pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2**2/(sigma2**2+delta**2))
    lower=sum(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2**2/(sigma2**2+delta**2))
    res=upper/lower
    return res

def sigma1new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    upper=sum((c-mu1)**2*pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1**4/(sigma1**2+delta**2)**2)
    lower=sum(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1**4/(sigma1**2+delta**2)**2)
    res=np.sqrt(upper/lower)
    return res

def sigma2new(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    upper=sum((c-mu2)**2*pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2**2/(sigma2**2+delta**2)**2)
    lower=sum(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2**2/(sigma2**2+delta**2)**2)
    res=np.sqrt(upper/lower)
    return res

def lkhood2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    res=pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2
    res=sum(np.log(res))
    return res

def bic1EM(c,delta,r):
    """
    Input: color, color errors, radius to the center
    BIC,mu,sigma
    """
    mu,sigma,aic,bic=ec.wstat(c,delta)
    lkhood_b=sum(np.log(pcz2(c,delta,r,mu,sigma)))
    BIC=-2.*lkhood_b + 2.*np.log(len(c))
    return BIC,mu,sigma
            

def bic2EM(c,delta,r,alpha,mu,sigma):
    """
    Input: color, color errors, radius to the center, initial guess of alpha, mu and sigma
    Output: BIC,alpha,mu,sigma
    
    """    
    w1=alpha[0]
    w2=alpha[1]
    mu1=mu[0]
    mu2=mu[1]
    sigma1=sigma[0]
    sigma2=sigma[1]
    NIter=100
    acc=0.00001
    for k in range(NIter):
        w1_new=w1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        w2_new=w2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma1_new=sigma1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma2_new=sigma2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu1_new=mu1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu2_new=mu2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_a=lkhood2(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_b=lkhood2(c,delta,r,w1_new,w2_new,mu1_new,mu2_new,sigma1_new,sigma2_new)
        alpha[0]=w1_new
        alpha[1]=w2_new
        mu[0]=mu1_new
        mu[1]=mu2_new
        sigma[0]=sigma1_new
        sigma[1]=sigma2_new
        if abs(lkhood_a - lkhood_b) <= acc:
            break
    BIC=-2.*lkhood_b + 5.*np.log(len(c))
    return BIC,alpha,mu,sigma
def nfwRandom(num=None):
    """
    generate random points following the nfw(r)
    """
    r=0.5
    scale=5.
    accept=0
    reject=0
    rr=np.zeros(num)
    i=0
    while accept < num:
        rt=r+np.random.randn()*scale
        if rt < 0 or rt > 1.:
            continue
        Ratio=min([1,nfw(rt)/nfw(r)])
        U=np.random.rand()
        if U <= Ratio:
            rr[i] = rt
            accept = accept + 1
            i = i+1
        else:
            reject = reject +1
    accRatio=accept/float(accept+reject)
    return rr,accRatio

def poissonRandom(num=None):
    """
    generate random points following the poisson(r)
    """
    u=np.random.rand(num)
    r=np.sqrt(u)
    return r
