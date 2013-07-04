%% Apple Mac OSX  
%mex  Jcost1.c  CLIBS="\$CLIBS -framework Accelerate" -largeArrayDims 

%% Locally compiled LAPACK with gfortran (Mac ports)
mex Jcost1.c CLIBS="\$CLIBS -L lapack/MACOSXgfortran  -llapacke -llapack -lrefblas  -L /opt/local/lib/gcc44 -lgfortran "

%% UHEM
%First compile LAPACK with ifort,icc and -fPIC
%National High Performance Computing Center of Turkey
%mex Jcost1.c CLIBS="\$CLIBS -L lapack/Linuxifort -llapacke -llapack -lrefblas -L /RS/progs/intel/lib/intel64 -lifcore -lirc "

