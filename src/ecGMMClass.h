/*---------------------------------------------------------------------------

  NAME:
    ecGMMClass.h
  
  
  PURPOSE:

    This C++ class implement the EM algorithm for the one dimensional 
    Error Corrected Gaussian Mixture Model as specified in 
    Jiangang Hao et al, 2008

  PLATFORM:

    This code is tested on Redhat Linux and Cygwin with gcc version 3.4.4
    in 2008. You can update your system to that time to make sure it works. 
  
  
  NOTE: since different people will use the mixture model very differnetly, I did
        not provide a simple function, instead, i put the most basic brick: the 
        ecGMMClass that you can customize very freely by following the example.
  

  DEPENDENCIES:  

    The class itself does not require other libraries than the standard C++
    librariers. If you want to use random number to test it, I recommend you 
    use the GNU gsl random number generators. In the ecGMM.cpp file, the random
    generator is used. So, make sure you have gsl properly installed. 

  USAGE: 
    This class produce an object, in which you can specify the data as 
    well as the EM iterations. The input part including:

    W : your initial guess of the weights of Gaussian components
    Mu: your initial guess of the locations of Gaussian componets
    Sigma: your initial guess of the width of Gaussian componets

    x: the data array
    xErr: the measurement errors corresponding to x

    N: the number of mixtures.

    Please refer to ecGMMexamle.cpp for how to use it. 
   


  REVISION HISTORY:
    Created September-2008: Jiangang Hao, University of Michigan
    9/15/2009: Bug in the sigma estimation is fixed. JG HAO
   

  Copyright (C) 2008 Jiangang Hao, Dept.of Phys., University of Michigan
                     jghao at umich dot edu  or jianganghao at gmail dot com

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


  ---------------------------------------------------------------------------*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

using namespace std;

class ecGMM

{
 public: 

  vector<double> x, xErr, W, Mu, Sigma;
  double GaussErrFun(double y, double yErr, double mean, double sd);
  double pSum(int j);
  double pIJ(int i, int j);
  void W_update();
  void Mu_update();
  void Sigma_update();
  void update();//replace W,Mu and Sigma with the new one
  void PartialUpdate(int fixMu, int fixSigma);//replace some of the W,Mu and Sigma
  double lhood();// Calculate the likelihood
  double BIC();// Calculate the BIC of the maximum likelihood
  double AIC();// Calculate the AIC of the maximum likelihood
  int N;  // number of mixture
  int M;  // number of data points
 
 
  double ep; // precision goal

  ecGMM(int N){for(int k=0; k<N; k++){Wt.push_back(0.); Mut.push_back(0.);Sigmat.push_back(0.);}}//constructor to initialize Wt, Mut, Sigmat

  private: 
  vector<double> Wt, Mut, Sigmat; // store the update
   
};


/*************definition of member functions******************/

double ecGMM::GaussErrFun(double y, double yErr, double mean, double sd)
{
  double PI=3.14159265;
  double res=0.;
  res = exp(- pow((y-mean),2.0)/(yErr*yErr+sd*sd)/2.0)/sqrt(2.0*PI*(yErr*yErr+sd*sd));
  return res; 
}


//--------define \sum_k=0^N p(y_j|z_j=k,\theta^{(t)})-------------------//
double ecGMM::pSum(int j)
{
  double res = 0.;
  for(int k = 0; k < N; k++)
    {
      res = res + GaussErrFun(x[j],xErr[j],Mu[k],Sigma[k])*W[k];
    }
  return res;
}


//-------define p(z_j=i|y_j,\theta^{(t)})--------------------//

double ecGMM::pIJ(int i, int j)

{
  double res;
  res = GaussErrFun(x[j],xErr[j],Mu[i],Sigma[i])*W[i]/pSum(j);
  return res;

}



//-------Update Weight and store it in Wt--------------------//


void ecGMM::W_update()
{
  double res;
  for(int i = 0; i < N; i++)
    {
      res = 0.;
      for(int j = 0; j < M; j++)
	{
          res = res + pIJ(i,j);
	}
      Wt[i] = res/double(M);
    }
}


//-------Update Mu and store it in Mut--------------------//

void ecGMM::Mu_update()
{
  double resu, resd;
  
  for(int i = 0; i < N; i++)
    {
      resu = 0.;
      resd = 0.;
 
      for(int j = 0; j < M; j++)
	{
          resu = resu + x[j]*pIJ(i,j)*pow(Sigma[i],2)/(pow(Sigma[i],2)+pow(xErr[j],2));
          resd = resd + pIJ(i,j)*pow(Sigma[i],2)/(pow(Sigma[i],2)+pow(xErr[j],2));
	}
      Mut[i] = resu/resd;
     }
 }

//-------Update Sigma and store it in Sigmat--------------------//

void ecGMM::Sigma_update()
{
   for(int i = 0; i < N; i++)
    {
     double resu = 0.;
     double resd = 0.;
      for(int j = 0; j < M; j++)
	{
          resu = resu + pow((x[j] - Mu[i]),2)*pIJ(i,j)*pow(pow(Sigma[i],2)/(pow(Sigma[i],2)+pow(xErr[j],2)),2);
          resd = resd + pIJ(i,j)*pow(Sigma[i],2)/(pow(Sigma[i],2)+pow(xErr[j],2));
	}
      Sigmat[i] = sqrt(resu/resd);
    }
}


//-------update all parameters --------------------
void ecGMM::update()
{
  for(int i=0; i<N; i++)
    {
      W[i] = Wt[i];
      Mu[i] = Mut[i];
      Sigma[i] = Sigmat[i];
    }
}

//-------Partial updates --------------------
void ecGMM::PartialUpdate(int fixMu, int fixSigma)
{
  for(int i=0; i<N; i++)
    {
      W[i] = Wt[i];
      if(i != fixMu){Mu[i] = Mut[i];}
      if(i != fixSigma){Sigma[i] = Sigmat[i];}
    }
}



//-----calculate the likelihood------------
double ecGMM::lhood()
{
  double res1, res2;

  res1=1.;
  for(int j=0;j<M;j++)
    {
      res2=0.;
      for(int i=0; i<N; i++)
	{
	  res2 = res2 + W[i]*GaussErrFun(x[j],xErr[j],Mu[i],Sigma[i]);
	}
      res1=res1*res2;
    }
  return(res1);
      
}


//-----return the BIC-----------------------

double ecGMM::BIC()
{
  double res;
  res=-2.0*log(lhood())+double(3.*N-1.)*log(double(M));
  return(res);
}


//-----return the AIC-----------------------

double ecGMM::AIC()
{
  double res;
  res=-2.0*log(lhood())+double(3.*N-1.)*2.;
  return(res);
}






