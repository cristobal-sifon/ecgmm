Error Corrected Gaussian Mixture Model is an extension of Gaussian Mixture Model when the data have Gaussian measurement errors.

# Introduction #

Gaussian Mixture Model (GMM) has wide applications in statistical analysis. However, traditional GMM consider only the case that each data point is measurement error free. In many applications, such as in astrophysics, the data points of interests normally have non-negligible measurement errors. In these situations, we need to model the measurement errors into the likelihood function of the GMM. This leads to a generalized GMM with measurement error corrections. We call it Error Corrected Gaussian Mixture Model (ECGMM). Its details and application to galaxy cluster analysis can be found in Hao et al, ApJ, 2009.

The ECGMM using EM algorithm is implemented in C++ and wrapped into a python package using SWIG. The method has been used to measure the properties of galaxy cluster red sequence and optical galaxy cluster detection (part of my PhD thesis). In the following, I will show how to get the code work on your computer. Please note that I am currently not able to provide any support for Win and Mac OS. I can only support linux. The codes have been successfully tested on ubuntu linux (32 and 64) and redhat linux (Scientific linux).


For details, please see this webpage:

https://sites.google.com/site/jiangangecgmm/