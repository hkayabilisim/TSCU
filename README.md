TSCU
====
Time Series Classification Utility (TSCU) is a MATLAB script
that is used to classify time series by using various distance
metrics and classification methods.

* Author: Huseyin Kaya, hkayabilisim@gmail.com, 2013/03/02
* Web Page: http://tscu.blogspot.com
* Source Page:  https://github.com/hkayabilisim/TSCU.git

SOURCE CODE
===========
Source files are available at https://github.com/hkayabilisim/TSCU.git

INSTALLATION
============
Please see http://tscu.blogspot.com

EXAMPLES
========
You can execute the test scripts "tscu-testXX.m" to see TSCU in action
or examine their published outputs on the web page.

DOCUMENTATION
=============
Please see http://tscu.blogspot.com

TODO
====
* Linear SVM will be added as an alternative classification scheme.

FEATURES
========
* Dynamic Time Warping, Constrained Dynamic Time Warping and Signal
  Alignment via Genetic Algorithm (SAGA) are available for candidate
  alignment methods.
* Paralel computing option is available. You can speed-up the computation
  if you enable the corresponding option and use K-NN classifier. 
* Classification performance measures such as user, producer accuracies, 
  estimated labels, Kappa statistics are provided in the function output.
* Displays input training and testing data with class labels.
* Nelder-Mead Simplex method can now be used as an optimization alternative
  if alignment is carried out with SAGA.

KNOWN ISSUES
============
* TSCU can not use multi-channel time series. 
* Althoug Jcost1 is MEX counterpart of Jcost0, it is still slower than Jcost0
  if 'SAGAOptimizationMethod' is 'Simplex'.

CHANGE LOG
==========
* Jcost0 now accepts the time variable (t). Making time discretization 
  outside of Jcost0 makes it approximately 5% faster. 2013-07-03
