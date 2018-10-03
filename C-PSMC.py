import csv
import os
import itertools
import math
import numpy as np
import subprocess
import time
import sys



### Specify sample PSMC file names, without file extensions (should be .psmc)
### 'Dependent samples' are ones whose demographic histories are to be scaled against an 'independent sampele'. Normally, these will be symbionts and hosts, respectively.
independent_samples = [sys.argv[1]]
dependent_samples = [sys.argv[2]]






##################################################################
##################################################################
####################     Comparative PSMC     ####################
##################################################################
#####     Luke Hecht, Peter Thompson, Benjamin Rosenthal     #####
##################################################################
##################################################################





sample_pairs = []
for target in independent_samples:
	for sample in dependent_samples:
		sample_pair = [target,sample]
		sample_pairs.append(sample_pair)



def parse_psmc(sample_name):
	infile = open('%s.psmc' % (sample_name), 'r')
	infile = csv.reader(infile, delimiter='\t')
	row_list = [row for row in infile]
	max_RD = [row for row in row_list if row[:1]==['MM']]
	max_RD = [row[1] for row in max_RD if row[1][:12]=='n_iterations']
	max_RD = str(max_RD)
	max_RD = int(max_RD[max_RD.index(':')+1:max_RD.index(',')])
	for row in row_list:
		if row[:2] == ['RD',str(max_RD)]:
			relevant_start = row_list.index(row)
	relevant_row_list = row_list[relevant_start:]
	RStimeseries = []
	for row in relevant_row_list: 
		if row[0] == 'RS':
			RStimeseries.append(row)
	return RStimeseries




### Curve-fitting
def find_scalars(independent_samples, dependent_samples):
	targets_inflection_points = []
	for target in independent_samples:
		psmc_output = parse_psmc(target)

		times = [float(row[2]) for row in psmc_output] ## can change back to coordinates from smoothed_coordinates
		Nes = [float(row[3]) for row in psmc_output]

		inflection_points = []
		for start in range(0,len(Nes),3):
			slic = Nes[start:start+10]
			if 2 < slic.index(min(slic)) < 7:
				inflection_points.append(start+slic.index(min(slic)))
			elif 2 < slic.index(max(slic)) < 7:
				inflection_points.append(start+slic.index(max(slic)))
		inflection_points = sorted(list(set(inflection_points)))
		inflection_points = [times[point] for point in inflection_points] ## associating 'true' time with inflection point
		#inflection_points = [pos_inflection_points,neg_inflection_points]
		targets_inflection_points.append(inflection_points)

	#####			
	#####
	samples_inflection_points = []
	for sample in dependent_samples:		
		psmc_output = parse_psmc(sample)

		times = [float(row[2]) for row in psmc_output] ## can change back to coordinates from smoothed_coordinates
		Nes = [float(row[3]) for row in psmc_output]

		inflection_points = []
		for start in range(0,len(Nes),3):
			slic = Nes[start:start+10]
			if 2 < slic.index(min(slic)) < 7:
				inflection_points.append(start+slic.index(min(slic)))
			elif 2 < slic.index(max(slic)) < 7:
				inflection_points.append(start+slic.index(max(slic)))
		inflection_points = sorted(list(set(inflection_points)))
		inflection_points = [times[point] for point in inflection_points] ## associating 'true' time with inflection point
		samples_inflection_points.append(inflection_points)


	sample_scalars = []
	for target in independent_samples:
		target_number = independent_samples.index(target)
		target_inflections = targets_inflection_points[target_number]
		#sample_scalars = []
		for sample in dependent_samples:
			sample_number = dependent_samples.index(sample)
			sample_inflections = samples_inflection_points[sample_number]
			scalars = []
			for target_inflection in target_inflections:
				for sample_inflection in sample_inflections:
					scalar = sample_inflection/target_inflection
					scalars.append(scalar)
			sample_scalars.append(scalars)
	return sample_scalars



def binary_slopes(times, Nes, timepoints):
	slopes = []
	for index in range(len(times)-1): ## can change back to coordinates from smoothed_coordinates
		for timepoint in timepoints:
			time0 = float(times[index])
			time1 = float(times[index+1])
			if time0 < timepoint <= time1:
				slope = (float(Nes[index]) - float(Nes[index+1]))/(time0-time1) ## slope per time-point (thousand years)
				if slope < 0:
					slope = 1
				elif slope > 0:
					slope = -1
				else:
					slope = 0
				slopes.append(slope)
	return slopes


				

def curvefitting(independent_samples, dependent_samples, scalars):
	for sample_pair in sample_pairs:
		print " ".join(sample_pair)
		target = sample_pair[0]
		sample = sample_pair[1]
		target_psmc_output = parse_psmc(target)
		sample_psmc_output = parse_psmc(sample)
		target_times = [float(row[2]) for row in target_psmc_output]
		sample_times = [float(row[2]) for row in sample_psmc_output]
		target_Nes = [float(row[3]) for row in target_psmc_output]
		sample_Nes = [float(row[3]) for row in sample_psmc_output]


		top5scalars = []
		scalar_set = scalars[sample_pairs.index(sample_pair)]
		for scalar in scalar_set:
			scaled_sample_times = [time*scalar for time in sample_times]
			max_timepoint = min([max(target_times), max(scaled_sample_times)])
			if max_timepoint == max(scaled_sample_times):
				timepoints = scaled_sample_times
			else: 
				timepoints = target_times
			#target_slopes = weight_slopes(target, binary_slopes(target_times, target_Nes, timepoints))
			#scaled_sample_slopes = weight_slopes(sample, binary_slopes(scaled_sample_times, sample_Nes, timepoints))
			target_slopes = binary_slopes(target_times, target_Nes, timepoints)
			scaled_sample_slopes = binary_slopes(scaled_sample_times, sample_Nes, timepoints)


			slope_differentials = []
			for index in range(len(timepoints)-1):
				#slope_differential = abs(target_slopes[index] - scaled_sample_slopes[index])
				if target_slopes[index] == 1 and scaled_sample_slopes[index] == 1:
					slope_differential = 0
				elif target_slopes[index] == -1 and scaled_sample_slopes[index] == -1:
					slope_differential = 0
				elif target_slopes[index] == 0 and scaled_sample_slopes[index] == 0:
					slope_differential = 0
				else: slope_differential = 2
				#slope_differential = abs(target_slopes[index] - scaled_sample_slopes[index])
				slope_differentials.append(slope_differential)
			top5scalars.append([sample_pair, float("%.3f" % (scalar)), int(np.mean(slope_differentials)*100)])
		top5scalars.sort(key=lambda x: x[2])
		for item in top5scalars[:10]: print item







curvefitting(independent_samples,dependent_samples,find_scalars(independent_samples,dependent_samples))










