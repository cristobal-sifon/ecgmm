# Error Corrected Gaussian Mixture Model #

Traditional Gaussian Mixture Model does not handle the measurement errors of each data point. In many applications, the data point themselves are uncertain to certain level and then a error corrected (or weighted if you would like) Gaussian Mixture Model is desirable.

In **_Hao et al, Astrophys.J.702:745 2009 (arXiv:0907.4383)_**, we introduced an error corrected Gaussian Mixture Model and an fast EM algorithm for its implementation. The technique has been applied to precision measurements of the evolution of the galaxy cluster red sequence.

In this project, I developed a C++ class to implement the error corrected Gaussian Mixture Model. A python version of the code was developed by wrapping the C++ function using SWIG.

Currently, it is only for 1 dimensional case. Higher dimensional extension is algorithmically trivial, though technically still need some special cares.

For more details, please see the wiki tab.