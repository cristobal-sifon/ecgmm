/*----------------------------------------------------------------
  NAME: ecGMMexample.cpp

  PURPOSE: demonstrate how to use ecGMMClass to performe 1 dim fitting. Firstly, 
           2 Gaussian data were generated and then fit it with 1,2 and 3 mixutres 
           respectively. Comparing their BIC, one can see the 2 Gaussian 
	   is favored. 
  NOTE: since different people will use the mixture model very differnetly, I did
        not provide a simple function, instead, i put the most basic brick: the 
        ecGMMClass that you can customize very freely by following the example.
        
 

  COMPILING: g++ -Wall ecGMMexample.cpp -lgsl -lgslcblas -lm -o ecGMMexample

  


  REVISION HISTORY:
    Created September-2008: Jiangang Hao, University of Michigan
   

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

 --------------------------------------------------------------------*/

#include "ecGMMClass.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int main()
{ 
  vector<double> x,xErr;
  double ep = 10e-20;    //set precision goal
  long maxIter = 1000;   //set the maximum iterations
  double mean,var,lh,BIC; 
  double MuCL = 0.5;
  double SdCL = 0.04;
  double MuBG = 0.5;
  double SdBG = 0.04;

  double SdErr = 0.;
  int NBG = 60;
  int NCL = 20;


  //--------------- random data generation-----------------------//

     const gsl_rng_type * T;
     gsl_rng * r;
     gsl_rng_env_setup();
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);

     
     for(int i = 0; i< NBG; i++)
       {
	 x.push_back(gsl_ran_gaussian(r,1.)*SdBG+MuBG);// background component

       }

     for(int i = 0; i< NCL; i++)
       {
	 x.push_back(gsl_ran_gaussian(r,1.)*SdCL+MuCL);//cluster component
	
       }

     for(int i=0; i< NBG+NCL; i++)
       {
       	 xErr.push_back(gsl_rng_uniform (r)*SdErr);
         x[i]=x[i]+gsl_ran_gaussian(r,1.)*xErr[i];
       }

     gsl_rng_free (r);
   
 //-----------------------foreword-------------------------
     cout << "------------------------------------------------" << endl;
     cout << "-     Welcome to the ECGMM example program     -" << endl;
     cout << "- Email Jiangang Hao: jghao at umich dot edu   -" << endl;
     cout << "- for help! (Last updated 9/21/2008)           -" << endl;
     cout << "------------------------------------------------" << endl;
  
     cout <<"    " << endl;
     
     cout <<"----------------True Values-----------------" << endl;
     cout << "True Background Number=" << NBG << endl;
     cout << "True Background Location=" << MuBG << endl;
     cout << "True Background Width=" << SdBG << endl;
     cout <<"    " << endl;
  
     cout << "True Cluster Number=" << NCL << endl;
     cout << "True Cluster Location=" << MuCL << endl;
     cout << "True Cluster Width=" << SdCL << endl;
     cout <<"    " << endl;
  
     cout << "Measurement Errors: N(0,1)*Unif(0,1)*" <<SdErr<< endl;
     cout <<"    " << endl;
  


  //----------------fitting with 1 componet Gaussian --------------//

      ecGMM one(1);              // create object
      
      for(int i=0;i<(int)x.size();i++)
        {
	  one.x.push_back(x[i]);
	  one.xErr.push_back(xErr[i]);  // load data
	 // one.xErr.push_back(0.);       // no error correction
	}

      one.W.push_back(1.);       // initialize parameters
      one.Mu.push_back(0.8);
      one.Sigma.push_back(0.4);
      one.M = one.x.size();
      one.N = 1;

     
      for(int k=0; k < maxIter; k++)  //EM iterations
	{
          double lhood = one.lhood();
	  one.W_update();
	  one.Mu_update();
	  one.Sigma_update();
	  one.update();
          if(log(one.lhood())-log(lhood)<= ep) 
	    {
	      cout << "------ Fit with one Gaussian Mixtures----" << endl;
	      cout << "W = " << one.W[0] << endl;
	      cout << "Mu = " << one.Mu[0] <<endl;
	      cout << "Sigma = " << one.Sigma[0] <<endl;
	      cout << "Maximum Likelihood=" <<one.lhood()<<endl;
	      cout << "BIC="<<one.BIC() <<endl;
	      cout <<"interations="<<k<< endl; 
	      break;
	    }
	}     
  
      
      

 //----------------fitting with 2 componet Gaussian --------------//

      ecGMM two(2);              // create object
      
      for(int i=0;i<(int)x.size();i++)
        {
	  two.x.push_back(x[i]);
	  two.xErr.push_back(xErr[i]);  // load data
	  //two.xErr.push_back(0.);       // no error correction
	}

      two.W.push_back(0.5);       // initialize parameters
      two.W.push_back(0.5);
      two.Mu.push_back(0.03);
      two.Mu.push_back(1.3);

      two.Sigma.push_back(0.4);
      two.Sigma.push_back(0.06);
      two.M = two.x.size();
      two.N = 2;


      

      for(int k=0; k < maxIter; k++)  //EM iterations
	{
      double lhood = two.lhood();
	  two.W_update();
	  two.Mu_update();
	  two.Sigma_update();
	  two.update();
          if(log(two.lhood())-log(lhood)<= ep) 
	    {
	      cout << "------ Fit with two Gaussian Mixtures----" << endl;
	      cout << "W = " << two.W[0] << " "<<two.W[1]<<endl;
	      cout << "Mu = " << two.Mu[0] <<" " << two.Mu[1] <<endl;
	      cout << "Sigma = " << two.Sigma[0] <<" " << two.Sigma[1]<<endl;
	      cout << "Maximum Likelihood=" <<two.lhood()<<endl;
	      cout << "BIC="<<two.BIC() <<endl;
	      cout <<"interations="<<k<< endl; 
	      break;
	    }
	}

   
  //----------------fitting with 3 componet Gaussian --------------//


      ecGMM three(3);              // create object
      
      for(int i=0;i<(int)x.size();i++)
        {
	  three.x.push_back(x[i]);
	  three.xErr.push_back(xErr[i]);  // load data
	  //three.xErr.push_back(0.);    //no error correction
	}

      three.W.push_back(0.3);       // initialize parameters
      three.W.push_back(0.3);
      three.W.push_back(0.4);

      three.Mu.push_back(0.1);
      three.Mu.push_back(1.1);
      three.Mu.push_back(2.1);

      three.Sigma.push_back(0.1);
      three.Sigma.push_back(0.6);
      three.Sigma.push_back(0.6);

      three.M = three.x.size();
      three.N = 3;

      for(int k=0; k < maxIter; k++)  //EM iterations
	{
      double lhood = three.lhood();
	  three.W_update();
	  three.Mu_update();
	  three.Sigma_update();
	  three.update();

	  if(log(three.lhood())-log(lhood)<= ep) 
	    {
	      cout << "------ Fit with three Gaussian Mixtures----" << endl;
	      cout << "W = " << three.W[0] << " "<<three.W[1]<< " " << three.W[2]<<endl;
	      cout << "Mu = " << three.Mu[0] <<" " << three.Mu[1] <<" " <<three.Mu[2]<<endl;
	      cout << "Sigma = " << three.Sigma[0] <<" " << three.Sigma[1]<< " " << three.Sigma[2] <<endl;
	      cout << "Maximum Likelihood=" <<three.lhood()<<endl;
	      cout << "BIC="<<three.BIC() <<endl;
          cout <<"interations="<<k<< endl; 
          break;
 
	    }

	}

  return(0);

}















  

 
