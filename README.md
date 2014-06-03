TSCU
====
Time Series Classification Utility (TSCU) is a MATLAB script
that is used to classify time series by using various distance
metrics and classification methods.

* Author: Huseyin Kaya, hkayabilisim@gmail.com, 2013/03/02
* Web Page: http://www.timewarping.org
* Source Page:  https://github.com/hkayabilisim/TSCU.git

SOURCE CODE
===========
Source files are available at https://github.com/hkayabilisim/TSCU.git

INSTALLATION
============
Please see http://www.timewarping.org

EXAMPLES
========
You can execute the test scripts "tscu-testXX.m" to see TSCU in action
or examine their published outputs on the web page.

DOCUMENTATION
=============
Please see http://www.timewarping.org

TODO
====
* export fig is not compatible with Octave

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
* Nelder-Mead Simplex, Genetic Algorithm Toolbox and a simplified MEX version of 
  Genetic Algorithm can be used as an optimization routine if you want to use
  SAGA alignment method.

KNOWN ISSUES
============
* TSCU can not use multi-channel time series. 
* export fig is not compatible with Octave

CHANGE LOG
==========
2014-06-04
* A simplified and fast MEX version of Genetic Algorithm is added.
* SAGACostfunction is removed.
* dtw.c is renamed as tscudtw.c
* export fig package and the curve registration toolbox of 
  Ramsay & Silverman is added

2013-07-03
* Jcost0 now accepts the time variable (t). Making time discretization 
  outside of Jcost0 makes it approximately 5% faster. 2013-07-03
