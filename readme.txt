ilm_compare
Version 1.3
William Wahby
2012

ilmc predicts the density of interconnects on a 3DIC via the Joyner method, then compares that result to a modified Joyner method which accounts for the placement of TSVs the array.
This version includes stable copies of both the matlab and python versions of the project. Both versions should produce identical results, but the python version runs faster.
The matlab version currently runs serially, while the python version is vectorized and uses a simple convolution operation to calculate Mt_intra_corr,
rather than the original (arcane) process (still used in the matlab code).

ilmc_python requires python 3.2 or better, numpy, and pylab/matplotlib. ilmc_matlab requires matlab 2007a or better.
Future development will take place primarily on the python version, with features being ported over to the MATLAB version on an as-needed basis.

To run, simply run ilm_comparison.py or test_comparison.m with the appropriate input parameters.

Currently the program is set up to run a test case used by Joyner in his PhD thesis.
