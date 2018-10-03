This program is intended to find the optimal scaling of a pair of <a href="https://github.com/lh3/psmc">PSMC</a> demographic history curves relative to one another such that contemporaneous slope differences are minimized, presumably reflecting the true parallel histories of species expected a priori to have been ecologically interdependent or under the same environmental pressures. See <a href="http://rspb.royalsocietypublishing.org/content/285/1888/20181032">Hecht, Thompson and Rosenthal (2018)</a> for a fuller explanation. 

C-PSMC.py takes as input the PSMC output files (.psmc) for two samples at a time. First, inflection points in both PSMC curves are identified. Then, one of the curves is held constant (“independent_sample”) while the other ("dependent_sample") is rescaled such that the timing (horizontal axis value) of one of its inflection points precisely aligns with that of an inflection point in the independent_sample curve. This rescaling is done for every pairwise combination of inflection points. For each scalar, the average difference in slope direction is calculated (for the overlapping portion of the curves). Scalar values represent relative evolutionary rate, the product of the between-sample ratios in mutation rate and generation time. For example, a scalar of 0.33 would imply 3 dependent_sample generations for every independent_sample generation, given equal mutation rates. The output result of C-PSMC.py is a list of the best relative scalars and their associated scores (lower is better).

**USAGE**

The only dependencies are Python (v2.7+) and the Numpy library for Python. The required variables are the file names (without extension) of the .PSMC output files for <independent_sample> and <dependent_sample>, and the desired number of top results <top_N> (default 5) sorted by slope differential. These are specified as command-line arguments:

> python C-PSMC.py <independent_sample> <dependent_sample> <top_N>


NOTE: The folder 'C-PSMC_original' contains the code used in the 2018 paper; however, we recommend using the version of C-PSMC.py published in the main directory for major improvements in run speed and ease of setup. A separate README can be found in C-PSMC_original.






