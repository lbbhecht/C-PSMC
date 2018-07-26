This program is intended to find the optimal scaling of a pair of PSMC demographic history curves relative to one another such that contemporaneous slope differences are minimized, presumably reflecting the true parallel histories of species expected a priori to have been ecologically interdependent or under the same environmental pressures. See Hecht, Thompson and Rosenthal (2018) for a fuller explanation.

C-PSMC.py takes as input the PSMC output files (.psmc) for two samples at a time. In short, one of the PSMC curves is held constant (“independent_sample”) and the average difference in slope direction of the other curve (“dependent_sample”) is calculated by aligning it to the first at a range of points along the curves. The dependent curve is then re-scaled by increasing the relative generation time and the average slope difference is recalculated. This process is repeated across a wide range of possible generation times resulting in a graph of average slope differences. The smallest average slope difference is the optimal generation time for the dependent species.

Sample files are included in this directory, and have already been plugged in with appropriate variables in C-PSMC.py. There are eight variables for the user to specify for each analysis, which are explained in detail below. These are set by editing C-PSMC.py directly using any plain text editor. The only dependencies are Python (v2.7+), R, and the Numpy library for Python. Once the variables are set, run the program like any other Python script. The following commands should work on any Unix system:

> cd */C-PSMC-master
>
> python C-PSMC.py

______________________________________


> independent_samples = [‘\<HostSampleName>’]  

This variable contains a list of the names of the samples to be held constant and compared against. Generally, these will be for a host species. For sample input file HomoSapiens.psmc, this variable should read independent_samples = ['HomoSapiens'].


> dependent_samples = [‘\<SampleName1>’,’\<SampleName2>’,…]  

As above, but a list of the file names for samples that will be scaled relative to the independent_samples; usually symbiont/parasite/pathogen species.


> timepoints = [‘\<Beginning_of_range>’,’\<End_of_range>’,’\<interval>’]  

This variable contains a list of three integers specifying the range of time points (in years) at which to compare slope agreement between the independent_samples and dependent_samples. For example, to compare the rate of change in effective population size at every thousand-year interval between the present and 400,000 years ago, the variable should read: timepoints = ['0','401000','1000']. Note that the quotations are necessary, as the values are read as strings. The oldest time point must not be older than the oldest one included in the reconstructed history of all independent_samples. The oldest time point can be easily read from the bottom of the <sample_name>.psmc.csv file that will be produced during the course of the run, if you do not already have an idea based on a raw PSMC plot. 



> grange = [‘\<Beginning_of_range>’,’\<End_of_range>’,’\<interval>’]  

Range of generation times (relative evolutionary rates) at which to score each dependent_sample relative to each independent_sample. For example, to consider generation times between one one-thousandth and twice as long as the independent_sample in intervals of one one-thousandth, set grange = ['0.001','2.001','0.001'].



> min_timepoints = \<value>  

This variable specifies the minimum number of the dependent_sample PSMC file's atomic time units (RS lines) required before calculating the average slope difference score for a given value of g. In effect, this sets a lower bound on the relative generation time. We recommend experimenting with different values if a priori estimates of relative generation time are lacking and if large differences are expected in the evolutionary rates of each comparator (e.g. human vs Drosophila). By default, this variable is set to 10.



> mutation_rates = [‘\<MutationRateIndependent>’,’\<MutationRateDependent1>’,’\<MutationRateDependent2>’,…]  

Assumed per-generation mutation rates, in the format '2.5e-8' per site per generation, for each independent and dependent sample, respectively. For example, for independent_samples = ['HomoSapiens'] and dependent_samples = ['Agamb','Pfalc'], we set mutation_rates = ['2.5e-8','1e-9','2.5e-8']. 



> independent_sample_scalars = [\<GenerationTimeInYears>]  

The generation time to assume for each independent sample, if known (if unknown, set to 1). Each time point for the independent_sample will be multiplied (stretched) by this factor before comparison with the dependent_sample. For example, for humans, this is set to independent_sample_scalars = [25], based on a priori estimates of 1 human generation per 25 years.



> plausible_ranges = [[\<min_GenerationTimeInYears>,\<max_GenerationTimeInYears>]]  

This variable specifies for each dependent_sample the plausible range of generation times, in years. For example, a mouse is unlikely to experience 100 generations in a year, although it may experience a few dozen. The purpose of the variable is simply to guide the analysis, helping to rule out anomalous good fits at extreme and implausible relative scales. As an example, correct usage when paired with dependent_samples = ['Agamb','Pfalc1'] (mosquito and malaria) would be: plausible_ranges = [[0.01,0.12],[0.3,0.8]]. This states that only generation times between 0.01 and 0.12 years (0.5-6 weeks) for mosquito, or 0.3-0.8 years for malaria, should be considered for 'best fit'.













