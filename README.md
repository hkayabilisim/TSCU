TSCU
====
Time Series Classification Utility (TSCU) is a MATLAB script
that is used to classify time series by using various distance
metrics and classification methods.

Author: Huseyin Kaya, hkayabilisim@gmail.com, 2013/03/02

Web Page: http://web.itu.edu.tr/~huseyinkaya/tscu
Source Page:  https://github.com/hkayabilisim/TSCU.git

SOURCE CODE
===========
Source files are available at https://github.com/hkayabilisim/TSCU.git

INSTALLATION
============
Please see http://web.itu.edu.tr/huseyinkaya/tscu.

EXAMPLES
========
I included some examples on web page. You can find more on the user manual.

DOCUMENTATION
=============
Please see http://web.itu.edu.tr/huseyinkaya/tscu and User Manual.

TODO
====
* Optimization subroutines should be added as an option.

FEATURES
========
* Dynamic Time Warping, Constrained Dynamic Time Warping and Signal
  Alignment via Genetic Algorithm (SAGA) are available for candidate
  alignment methods.
* Paralel computing option is available. You can speed-up the computation
  if you enable it and use K-NN classifier. 
* Classification performance measures such as user, producer accuracies, 
  estimated labels, Kappa statistics are provided in the function output.
* Displays input training and testing data with class labels.
KNOWN ISSUES
============
* TSCU can not use multi-channel time series. 
* Althoug Jcost1 is MEX counterpart of Jcost0, it is still slower than Jcost0.
