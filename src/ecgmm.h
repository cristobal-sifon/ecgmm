/*----------------------------------------------------------------------------
In this .h file, I defined two CONVENIENT functions based on the ecGMMClass.h. 
To use these functions, you must include ecGMMClass.h and ecgmm.h in your main progrom. The two functions, BICecgmm and AICecgmm implement the ECGMM and return
BIC and AIC. 

9/20/2008, Jiangang Hao, Dept. of Physics, University of Michigan, Ann Arbor
-----------------------------------------------------------------------------*/

#include "ecGMMClass.h"
using namespace std;

//----- function to calculate the ECGMM parameters and BIC-----------


double BICecgmm(vector<double> &x, vector<double> &xErr,vector<double> &alpha,vector<double> &mu,vector<double> &sigma)
{

  double ep = 10e-20;    //set precision goal
  long maxIter = 1000;   //set the maximum iterations

  ecGMM f((int)alpha.size());
  f.x = x;
  f.xErr = xErr;
  f.W = alpha;
  f.Mu = mu;
  f.Sigma = sigma;
  f.M=(int)x.size();    // # of data points
  f.N=(int)alpha.size();// # of mixutures

  for(int k=0;k<maxIter;k++)
    {
      double lhood = f.lhood();
      f.W_update();
      f.Mu_update();
      f.Sigma_update();
      f.update();

      if(log(f.lhood())-log(lhood)<= ep)
	{
	  alpha = f.W;
          mu = f.Mu;
          sigma = f.Sigma;
	  break;
          
	}

    }
   return(f.BIC());
}


//----- function to calculate the ECGMM parameters and AIC-----------

double AICecgmm(vector<double> &x, vector<double> &xErr,vector<double> &alpha,vector<double> &mu,vector<double> &sigma)
{

  double ep = 10e-20;    //set precision goal
  long maxIter = 1000;   //set the maximum iterations

  ecGMM f((int)alpha.size());
  f.x = x;
  f.xErr = xErr;
  f.W = alpha;
  f.Mu = mu;
  f.Sigma = sigma;
  f.M=(int)x.size();    // # of data points
  f.N=(int)alpha.size();// # of mixutures

  for(int k=0;k<maxIter;k++)
    {
      double lhood = f.lhood();
      f.W_update();
      f.Mu_update();
      f.Sigma_update();
      f.update();

      if(log(f.lhood())-log(lhood)<= ep)
	{
	  alpha = f.W;
          mu = f.Mu;
          sigma = f.Sigma;
	  break;
          
	}

    }
   return(f.AIC());
}


