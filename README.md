This program is intended to find the optimal way of scaling a pair of PSMC demographic history curves relative to one another such that contemporaneous slope differences are minimized, presumably reflecting the true parallel histories of species expected a priori to have been ecologically interdependent or under the same environmental pressures. See Hecht, Thompson and Rosenthal (2018) for a fuller explanation.

C-PSMC.py takes as input the PSMC output files for two samples at a time. In short, one of the PSMC curves is held constant while the fit of the other is scored by aligning it to the first at a range of points for a range of scalar values. Sample files can be found in this directory, and are already plugged in with appropriate variables in C-PSMC.py. 

There are eight variables for the user to specify for each analysis, which are explained in detail below. These are set by editing C-PSMC.py directly. This can be done using any plain text editor, but one designed for script editing, like Sublime Text, is recommended. Once the variables are set, run the program like any other Python script. The following commands should work on any Unix system:

> cd */C-PSMC-master

python C-PSMC.py

______________________________________


independent_samples =
This variable contains a list of the names of the samples to be held constant and compared against. Generally, these will be for a host species. For sample input file HomoSapiens.psmc, this variable should read independent_samples = ['HomoSapiens'].


dependent_samples = 
As above, but a list of the file names for samples that will be scaled relative to the independent_samples; usually symbiont/parasite/pathogen species.


timepoints = 
This variable contains a list of three integers specifying the range of time points (in years) at which to compare slope agreement between the independent_samples and dependent_samples. For example, to compare the rate of change in effective population size at every thousand-year interval between the present and 400,000 years ago, the variable should read: timepoints = ['0','401000','1000']. Note that the quotations are necessary, as the values are read as strings. The oldest time point must not be older than the oldest one included in the reconstructed history of all independent_samples. The oldest time point can be easily read from the bottom of the <sample_name>.psmc.csv file that will be produced during the course of the run, if you do not already have an idea based on a raw PSMC plot. 


grange = 
Range of generation times (relative evolutionary rates) at which to score each dependent_sample relative to each independent_sample. For example, to consider generation times between one one-thousandth and twice as long as the independent_sample in intervals of one one-thousandth, set grange = ['0.001','2.001','0.001'].


min_timepoints = 
This variable specifies the minimum number of the dependent_sample PSMC file's atomic time units (RS lines) need to be involved in the calculation of the average slope difference score for a given value of g. This is basically an arbitrary threshold for statistical significance, but effectively sets a lower bound on the relative generation time. We recommend experimenting with different values if the species you are comparing are expected to have very different evolutionary rates (e.g. human vs Drosophila). By default, this variable is set to 10.


mutation_rates = 
Assumed per-generation mutation rates, in the format '2.5e-8' per site per generation, for each independent and dependent sample, respectively. For example, for independent_samples = ['HomoSapiens'] and dependent_samples = ['Agamb','Pfalc'], we set mutation_rates = ['2.5e-8','1e-9','2.5e-8']. 


independent_sample_scalars = 
The generation time to assume for each independent sample, if known (if unknown, set to 1). Each time point for the independent_sample will be multiplied (stretched) by this factor before comparison with the dependent_sample. For example, for humans, this is set to independent_sample_scalars = [25].


plausible_ranges = 
This variable specifies for each dependent_sample the plausible range of generation times, in years. For example, a mouse is unlikely to experience 100 generations in a year, although it may experience a few dozen. This variable simply guides the analysis, helping to rule out anomalous good fits at extreme and implausible relative scales. Correct usage 













