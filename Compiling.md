Polymer libraries is developed using g++ (version >=4.6).

Some routine used in polyerm\_lib.h use LAPACK and FFTW3.  Moreover polyerm\_lib.h use routine that were independently developed by Julien Dorier.


After having included polymer\_lib.h in your code you can compile with:

**g++ name\_of\_your\_code.cpp Tstatistic.cpp -O3 -fopenmp -llapack -lfftw3 -o nme\_of\_the\_output\_executable**

It is strongly suggested to put the various libraries and option after the files that need to be compile, in fact some system does not recognize the option if these are put before the files.