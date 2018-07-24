import csv
import os
import itertools
import math
import numpy as np
import subprocess
import time



### Specify sample PSMC file names, without file extensions (should be .psmc)
### 'Dependent samples' are ones whose demographic histories are to be scaled against an 'independent sampele'. Normally, these will be symbionts and hosts, respectively.
independent_samples = ['HomoSapiens']
dependent_samples = ['Agamb','Pfalc1']



### Curve-fitting range settings
timepoints = ['0','401000','1000'] ## ['min','max','increment'] ## Range of time points (in years) at which to compare slope agreement over PSMC-reconstructed history. The oldest time point must not be older than the oldest one included in the reconstructed history of all independent_samples.
grange = ['0.001','2.001','0.001'] ## ['min','max','increment'] ## Range of generation times at which to score each dependent_sample relative to each independent_sample.



### Optimal scalar (g) constraints
min_timepoints = 10 ## What minimum number of the PSMC file's atomic time units (RS lines) does the average slope difference score for a particular g need to be based on?
mutation_rates = ['2.5e-8','1e-9','2.5e-8'] ## Assumed per-generation mutation rates, in the format '2.5e-8', for each independent and dependent sample, respectively, in the above lists.
independent_sample_scalars = [25] ## The generation time to assume for each independent sample, if known (if unknown, set == 1). Each time point will be multiplied (stretched) by this factor.
plausible_ranges = [[0.01,0.12],[0.3,0.8]] ## The plausible range of generation times, in years, for each dependent sample (min,max).




##################################################################
##################################################################
####################     Comparative PSMC     ####################
##################################################################
#####     Luke Hecht, Peter Thompson, Benjamin Rosenthal     #####
##################################################################
##################################################################




### Generating base scaled .psmc.csv
def extract_par(sample_name):
	outfile = open('%s.psmc.csv' % (sample_name), 'w')
	infile = open('%s_plot.0.txt' % (sample_name), 'r')
	infile = csv.reader(infile, delimiter='\t')
	row_list = []
	for row in infile:
		row_list.append(row[:2])
	for row in row_list:
		outfile.write('%s,%s\n' % (row[0], row[1]))
	outfile.close()


### Curve-fitting
def curvefitting(independent_samples,dependent_samples,grange,timepoints):
	timepoints_start = int(timepoints[0])
	timepoints_end = int(timepoints[1])
	timepoints_increment = int(timepoints[2])
	timepoints = range(timepoints_start,timepoints_end,timepoints_increment)

	grange_start = float(grange[0])
	grange_end = float(grange[1])
	grange_increment = float(grange[2])
	garray = np.arange(grange_start,grange_end,grange_increment)

	global curvefit_output_files
	curvefit_output_files = []
	for independent_sample in independent_samples:
		targets = [independent_sample]
		target_slopes = []
		target_timepoints_atomic_datapoints = []
		for target in targets:
			timepoints_atomic_datapoints = []
			slopes = []
			coordinates = open('%s.psmc.csv' % (target), 'rb')
			coordinates = csv.reader(coordinates, delimiter=',')
			coordinates = [row for row in coordinates]
			for index in range(len(coordinates)-1):
				time0 = float(coordinates[index][0])
				time1 = float(coordinates[index+1][0])
				slope = (float(coordinates[index][1]) - float(coordinates[index+1][1]))/(time0-time1) ## slope per time-point (thousand years)
				count = 0
				for timepoint in timepoints:
					if time0 <= int(timepoint) < time1:
						slopes.append([timepoint,slope])
						count = 1
				timepoints_atomic_datapoints.append(count)
			for slope in slopes:
				if slope[1] == 0:
					slope[1] = 0
				elif slope[1] < 0:
					slope[1] = 1
				elif slope[1] > 0:
					slope[1] = -1
			target_slopes.append(slopes)
			target_timepoints_atomic_datapoints.append(timepoints_atomic_datapoints)
					
		for dependent_sample in dependent_samples:
			samples = [dependent_sample]
			gfactor_sample_timepoints_atomic_datapoints = []	
			gfactor_sample_slopes = []
			for gfactor in garray:
				sample_slopes = []
				sample_timepoints_atomic_datapoints = []
				for sample in samples:
					timepoints_atomic_datapoints = []
					slopes = []
					coordinates = open('%s.psmc.csv' % (sample), 'rb')
					coordinates = csv.reader(coordinates, delimiter=',')
					coordinates = [row for row in coordinates]
					for index in range(len(coordinates)-1):
						time0 = float(coordinates[index][0])*gfactor
						time1 = float(coordinates[index+1][0])*gfactor
						slope = (float(coordinates[index][1]) - float(coordinates[index+1][1]))/(time0-time1)
						count = 0
						for timepoint in timepoints:
							if time0 <= int(timepoint) < time1:
								slopes.append([timepoint,slope])
								count = 1
						timepoints_atomic_datapoints.append(count)
					for slope in slopes:
						if slope[1] == 0:
							slope[1] = 0
						elif slope[1] < 0:
							slope[1] = 1
						elif slope[1] > 0:
							slope[1] = -1
					sample_slopes.append(slopes)
					sample_timepoints_atomic_datapoints.append(timepoints_atomic_datapoints)
				gfactor_sample_slopes.append(sample_slopes)
				gfactor_sample_timepoints_atomic_datapoints.append(sample_timepoints_atomic_datapoints)


			gfactor_scores = []
			count = grange_start
			for sample_slopes in gfactor_sample_slopes:
				gfactor_avglen = []
				gfactor_avgdiff = []
				gfactor_avgdata = []
				samples_atomic_datapoints = gfactor_sample_timepoints_atomic_datapoints[gfactor_sample_slopes.index(sample_slopes)]
				for target_series in target_slopes:
					target_atomic_datapoints = target_timepoints_atomic_datapoints[target_slopes.index(target_series)]
					for sample_series in sample_slopes:
						sample_atomic_datapoints = samples_atomic_datapoints[sample_slopes.index(sample_series)]
						target_x_sample_diff = []
						for index in range(len(sample_series)):
							slope_differential = abs(target_series[index][1] - sample_series[index][1])
							target_x_sample_diff.append(slope_differential)
						gfactor_avglen.append(len(sample_series))
						gfactor_avgdiff.append(np.mean(target_x_sample_diff))
						gfactor_avgdata.append(np.mean((sum(sample_atomic_datapoints))))
				gfactor_scores.append('%.3f,%.4f,%i,%i' % (count, np.mean(gfactor_avgdiff), np.mean(gfactor_avgdata), np.mean(gfactor_avglen)))
				count += grange_increment

			### curvefit_output
			curvefit_output_file = '%s_x_%s.curvefit.csv' % (independent_sample, dependent_sample)
			curvefit_output_files.append(curvefit_output_file)
			curvefit_output_file = open(curvefit_output_file, 'w')
			for item in gfactor_scores:
				curvefit_output_file.write(item+'\n')
			curvefit_output_file.close()
	return curvefit_output_files


### Finding top n generation times
def best_fits(independent_samples,dependent_samples,top_n,input_file_names):
	all_targets = independent_samples
	all_samples = dependent_samples

	target_index = -1
	for target in all_targets:
		targets = [target]
		target_index += 1
		target_g_multiplier = independent_sample_scalars[target_index]
		sample_index = -1
		for sample in all_samples:
			samples = [sample]
			sample_index += 1
			grange = plausible_ranges[sample_index]

			infile = input_file_names[all_samples.index(sample)]
			infile = open(infile, 'rb')
			lines = csv.reader(infile, delimiter=',')
			lines = [row for row in lines]
			lines = [row for row in lines if float(row[2]) >= min_timepoints]
			lines = [row for row in lines if grange[0] <= float(row[0]) <= grange[1]]
			gs = [float(row[0]) for row in lines]
			scores = [float(row[1]) for row in lines]
			best_scores = []
			optimal_gs = []
			for iteration in range(top_n):
				if len(scores) > 0:
					best_score = min(scores)
					optimal_g = gs[scores.index(best_score)]
					best_scores.append(best_score)
					optimal_gs.append(optimal_g)
					lines = [row for row in lines if float(row[0]) <= optimal_g*0.9 or float(row[0]) >= optimal_g*1.1]
					gs = [float(row[0]) for row in lines]
					scores = [float(row[1]) for row in lines]
				else:
					best_score = 2
					optimal_g = 0
					best_scores.append(best_score)
					optimal_gs.append(optimal_g)
				
		
			### optimal_output
			print_list = [target,sample]
			print_list += [[optimal_gs[index],best_scores[index]] for index in range(top_n)]
			print('\t'.join(map(str,print_list)))


##################################################################
##################################################################


### extracting data from .psmc
index = 0
for independent_sample in independent_samples:
	os.system('./psmc_plot.pl -R -u %s -g %i -x 100 -Y 10 %s_plot %s.psmc' % (mutation_rates[index], independent_sample_scalars[independent_samples.index(independent_sample)], independent_sample, independent_sample))
	extract_par(independent_sample)
	index += 1
	for dependent_sample in dependent_samples:
		os.system('./psmc_plot.pl -R -u %s -g 1 -x 100 -Y 10 %s_plot %s.psmc' % (mutation_rates[index], dependent_sample, dependent_sample))
		extract_par(dependent_sample)
		index += 1

### curve-fitting and finding best fits
best_fits(independent_samples,dependent_samples,5,curvefitting(independent_samples,dependent_samples,grange,timepoints))
print

### producing R plots of fit x g
R_args = [output_file for output_file in curvefit_output_files]
os.system('rscript curvefit_plot.r %i ' % (float(grange[1])-float(grange[2])) + ' '.join(R_args))







