#Calculate read ratios at each position from treated and untreated sample files using a sliding window approach. Get highest peak in the region of lowest pvalue from fishers test by the highest Rho readtrough score
import datetime
import scipy.stats as stats
import pathlib
import os
import os.path
import json
print('WARNING: may take 3-4h to score and select peaks with E. coli genome ~4.5 Mbp')
#get parameters from json config file in the smae folder
cwd = str(pathlib.Path(__file__).parent.absolute())
with open(os.path.join(cwd, 'Rho-termination_peak-caller_Input-parameters.json')) as f:
  parameters = json.load(f)

#assign parameters
window = int(parameters["window"])  # window size for summing ups/ds read coverage
peak_window = int(parameters["single_peak_window"])   # min window size between 2 clusters of significant transcriptional radthrough (for calling a signle peak in the region)
pval_threshold = float(parameters["p_value_threshhold"]) #fisher's test max p value threshold
ratio_threshold = int(parameters["ratio_threshhold"])  # 1/(untreated/treated) ratio direction, ratio > 1 means treated sample has more readthrough
genome_end_coord = int(parameters["genome_size"])  # aka genome size
coord_col = int(parameters["coordinate_column"]) - 1 #coordinate column in the input file
read_count_col = int(parameters["read_count_column"]) - 1 #read count column in the input file
try:
	strand_col = int(parameters["optional_strand_column"]) - 1 #gff files will have negative values = no strand column, others might have a seprate col for strand
except:
	strand_col = ""
#modify input/output file paths/names depending on what was entered into the parameters file
input_dir = parameters["optional_input_file_directory"]
if input_dir == "":
	input_dir = cwd
untreated_file = os.path.join(input_dir, parameters["untreated_file_name"])
treated_file = os.path.join(input_dir, parameters["treated_file_name"])
out_dir = parameters["optional_output_file_directory"]
if out_dir == "":
	out_dir = cwd
peak_file = parameters["optional_out_peak_file_name"]
if peak_file == "":
	peak_file = 'Rho_peaks_within_' + str(window) + 'nt_below_pval' + str(pval_threshold) + '.txt'
peak_file = os.path.join(out_dir, peak_file)
window_sum_file = parameters["optional_out_sum-score_file_name"]
if window_sum_file == "":
	window_sum_file='Summary-Rho_score&win_sums_in_' + str(window) + '_below_pval' + str(pval_threshold) + 'nt.txt'
window_sum_file = os.path.join(out_dir, window_sum_file)

#check that file paths are corect for inputs
for f in [untreated_file, treated_file]:
	if os.path.exists(os.path.join(input_dir, f)) == False:
		raise Exception('Check input file directory (leave blank if in the same folder as this program) and/or file names in Input_parameters.json file; need double-backslash for windows')
	if os.path.exists(out_dir) == False:
		raise Exception('Check output file directory (leave blank if in the same folder as this program) in Input_parameters.json file')

print('Start:', datetime.datetime.now())

#read input file and save coordinates and associated read count
def read_gff(file_name):
	with open(file_name, 'r') as f: #creating dictionaries for coord:count pair for each strand separately
		plus_strand = {}
		minus_strand = {}
		if strand_col == "": #no strand col = minus strand are represented as negative values of read count
			for line in f:
				splitline = line.rstrip().split('\t')
				coord = int(splitline[coord_col])
				count = float(splitline[read_count_col])
				if count < 0:
					minus_strand.update({coord:count})
				elif count > 0:
					plus_strand.update({coord: count})
		elif strand_col != "":  #when strand colum is present and count values are opsitive regardless of the strand
			for line in f:
				splitline = line.rstrip().split('\t')
				coord = int(splitline[coord_col])
				if splitline[strand_col] == '-':
					count = float(splitline[read_count_col]) * -1
					minus_strand.update({coord: count})
				elif splitline[strand_col] == '+':
					count = float(splitline[read_count_col])
					plus_strand.update({coord: count})
	return plus_strand, minus_strand

#Not every single coordinate from 1 to x is represented because no value was associated with it, check what's missing
#fix by adding Key for missing coord and value of 0 for each strand dict
def restore_zeros(coord_val_dict, single_nested):
	if single_nested == 'single':
		for i in range(1, genome_end_coord + 1):
			if i in coord_val_dict:
				continue
			else:
				coord_val_dict.update({i: 0})
	elif single_nested == 'nested':
		for i in range(1, genome_end_coord + 1):
			if i in coord_val_dict:
				continue
			else:
				coord_val_dict.update({i: [0,0,0,0,0,0]})
	return coord_val_dict

#make a dictionary that contains only 1 peak within xnt window, do one strand at the time
#get was usefult when don't want to repoppulate the list with 0, but wouldn't work when dict val is a list
def find_local_peaks(cluster_dict):
	local_peaks = {}
	position = 1
	while position <= len(cluster_dict):
		if cluster_dict.get(position)[4] == 0: #pvalue is 4th in the list
			position = position + 1
		elif cluster_dict.get(position)[4] != 0: #looking for the begining of ratios cluster
			peaks = {}
			while cluster_dict.get(position)[4] != 0: #merge cluster into one dict
				peaks.update({position:cluster_dict.get(position)[4]})
				position = position + 1
				if cluster_dict.get(position)[4] == 0:
					for x in range(1, peak_window + 1): #merge until the last peak in the cluster is more than 500nt away
						try: #added this to fix getting to the position close to the end of the chromosome
							if cluster_dict.get(position + x)[4] == 0:
								continue
							else:
								position = position + x
								break
						except TypeError:
							break
				if position > len(cluster_dict):
					break
			max_peak = max(peaks.items(), key=lambda y:y[1])
			local_peaks.update({max_peak[0]: cluster_dict[max_peak[0]]})
	return local_peaks

untreated = read_gff(untreated_file)
plus_strand_none = restore_zeros(untreated[0], 'single')
minus_strand_none = restore_zeros(untreated[1], 'single')
treated = read_gff(treated_file)
plus_strand_BCM = restore_zeros(treated[0], 'single')
minus_strand_BCM = restore_zeros(treated[1], 'single')

print('updated 3end list:', datetime.datetime.now())

#sorts the dict by key, dic[key] returns value for that key in a list (and turns minus strand values to + values)
plus_sort_none = list([plus_strand_none[key]for key in sorted(plus_strand_none)])
minus_sort_none = list([minus_strand_none[key]*-1 for key in sorted(minus_strand_none)])
plus_sort_BCM = list([plus_strand_BCM[key] for key in sorted(plus_strand_BCM)])
minus_sort_BCM = list([minus_strand_BCM[key]*-1 for key in sorted(minus_strand_BCM)])

#save memory by clearing dictionaries that won't be used
del(plus_strand_none)
del(minus_strand_none)
del(plus_strand_BCM)
del(minus_strand_BCM)
#calculate: ratios of none/BCM which are ratios of sums in a specific window up and down from i

with open(window_sum_file, 'w') as w:
	w.write('Coordinate\tStrand\tuntreated_us\tuntreated_ds\tTreated_us\tTreated_ds\tReadthrough Score\tSignificance score(p_value)\n')
	up = genome_end_coord - window
	for coord in range(0, genome_end_coord):
		if window - 1 <= coord < up:
			ds_plus_none = float(sum(plus_sort_none[coord + 1: coord + 1 + window]))
			us_plus_none = float(sum(plus_sort_none[coord + 1 - window: coord + 1]))
			ds_plus_BCM = float(sum(plus_sort_BCM[coord + 1: coord + 1 + window]))
			us_plus_BCM = float(sum(plus_sort_BCM[coord + 1 - window:coord + 1]))
		elif coord < window - 1:  # below this range sum of the window positions before i will be out of range,so in this case add missing postions from the end of the list (list is 'circular')
			ds_plus_none = float(sum(plus_sort_none[coord + 1: coord + 1 + window]))
			us_plus_none = float(sum(plus_sort_none[0: coord + 1]) + sum(plus_sort_none[genome_end_coord-window+coord+1:]))
			ds_plus_BCM = float(sum(plus_sort_BCM[coord + 1: coord + 1 + window]))
			us_plus_BCM = float(sum(plus_sort_BCM[0: coord + 1])+sum(plus_sort_BCM[genome_end_coord-window+coord+1:]))
		elif coord >= up:  # above this range sum of the 'window' positions after i will be out of range,so in this case add missing postions from the beginning of the list (list is 'circular')
			ds_plus_none = float(sum(plus_sort_none[coord + 1:])+sum(plus_sort_none[0:window  + coord - genome_end_coord + 1]))
			us_plus_none = float(sum(plus_sort_none[coord - window + 1:coord + 1]))
			ds_plus_BCM = float(sum(plus_sort_BCM[coord + 1:])+sum(plus_sort_BCM[0:window + coord - genome_end_coord + 1]))
			us_plus_BCM = float(sum(plus_sort_BCM[coord - window + 1:coord + 1]))
		try:
			ratio_plus = 1 / ((ds_plus_none / us_plus_none) / (ds_plus_BCM / us_plus_BCM))
			zero_error_plus = 'No'
		except ZeroDivisionError:  #deal with 0 values
			zero_error_plus = 'Yes'
			temp_sums_plus = {'ds_plus_none':ds_plus_none, 'us_plus_none':us_plus_none, 'ds_plus_BCM':ds_plus_BCM, 'us_plus_BCM':us_plus_BCM} #does't replace origina values if creating a dictionary for them
			for name, sums in temp_sums_plus.items():
				if sums == 0:
					temp_sums_plus[name] = 0.0001
			ratio_plus = 1 / ((temp_sums_plus['ds_plus_none'] / temp_sums_plus['us_plus_none']) / (temp_sums_plus['ds_plus_BCM'] / temp_sums_plus['us_plus_BCM'])) #using dict values instead of originals
		if ratio_plus > ratio_threshold: #this will take care of change in expresion in the wrong direction
			oddsratio, pvalue_plus = stats.fisher_exact([[ds_plus_none, us_plus_none],[ds_plus_BCM, us_plus_BCM]])
			if pvalue_plus < pval_threshold:
				if zero_error_plus == 'No':
					w.write(str(coord + 1) + '\t+\t' + str(us_plus_none) + '\t' + str(ds_plus_none) + '\t' + str(us_plus_BCM) + '\t' + str(ds_plus_BCM) + '\t' + str(ratio_plus) + '\t' + str(pvalue_plus) + '\n')
				elif zero_error_plus == 'Yes': #'adjusted' readtrough score goes into peak calling function but will write N/A in the file
					w.write(str(coord + 1) + '\t+\t' + str(us_plus_none) + '\t' + str(ds_plus_none) + '\t' + str(us_plus_BCM) + '\t' + str(ds_plus_BCM) + '\tN/A\t' + str(pvalue_plus) + '\n')
		if window <= coord <= up:
			ds_minus_none = float(sum(minus_sort_none[coord - window: coord]))
			us_minus_none = float(sum(minus_sort_none[coord: coord + window]))
			ds_minus_BCM = float(sum(minus_sort_BCM[coord - window: coord]))
			us_minus_BCM = float(sum(minus_sort_BCM[coord: coord + window]))
		elif coord < window:
			ds_minus_none = float(sum(minus_sort_none[0:coord])+sum(minus_sort_none[genome_end_coord-window+coord:]))
			us_minus_none = float(sum(minus_sort_none[coord:coord+window]))
			ds_minus_BCM = float(sum(minus_sort_BCM[0:coord])+sum(minus_sort_BCM[genome_end_coord-window+coord:]))
			us_minus_BCM = float(sum(minus_sort_BCM[coord:coord + window]))
		elif coord > up:
			ds_minus_none = float(sum(minus_sort_none[coord-window:coord]))
			us_minus_none = float(sum(minus_sort_none[coord:])+sum(minus_sort_none[0:window + coord - genome_end_coord]))
			ds_minus_BCM = float(sum(minus_sort_BCM[coord-window:coord]))
			us_minus_BCM = float(sum(minus_sort_BCM[coord:])+sum(minus_sort_BCM[0:window + coord - genome_end_coord]))
		try:
			ratio_minus = 1 / ((ds_minus_none / us_minus_none) / (ds_minus_BCM / us_minus_BCM))
			zero_error_minus = 'No'
		except ZeroDivisionError:
			zero_error_minus = 'Yes'
			temp_sums_minus = {'ds_minus_none':ds_minus_none, 'us_minus_none':us_minus_none, 'ds_minus_BCM':ds_minus_BCM, 'us_minus_BCM':us_minus_BCM} #does't replace origina values if creating a dictionary for them
			for name, sums in temp_sums_minus.items():
				if sums == 0:
					temp_sums_minus[name] = 0.0001
			ratio_minus = 1 / ((temp_sums_minus['ds_minus_none'] / temp_sums_minus['us_minus_none']) / (temp_sums_minus['ds_minus_BCM'] / temp_sums_minus['us_minus_BCM']))  # using dict values instead of originals
		if ratio_minus > ratio_threshold:
			oddsratio, pvalue_minus = stats.fisher_exact([[ds_minus_none, us_minus_none],[ds_minus_BCM, us_minus_BCM]])
			if pvalue_minus < pval_threshold:
				if zero_error_minus == 'No':
					w.write(str(coord + 1) + '\t-\t' + str(us_minus_none) + '\t' + str(ds_minus_none) + '\t' + str(us_minus_BCM) + '\t' + str(ds_minus_BCM) + '\t' + str(ratio_minus) + '\t' + str(pvalue_minus) + '\n')
				elif zero_error_minus == 'Yes':
					w.write(str(coord + 1) + '\t-\t' + str(us_minus_none) + '\t' + str(ds_minus_none) + '\t' + str(us_minus_BCM) + '\t' + str(ds_minus_BCM) + '\tN/A\t' + str(pvalue_minus) + '\n')

del(plus_sort_none)
del(minus_sort_none)
del(plus_sort_BCM)
del(minus_sort_BCM)

print('Finished calculating ratios:', datetime.datetime.now())
peak_plus, peak_minus = {}, {}
div_zero_error_plus, div_zero_error_minus = [], []
with open(window_sum_file, 'r') as r:
	next(r)
	for line in r:
		splitline = line.rstrip().split('\t')
		strand = splitline[1]
		ratio = splitline[6]
		coord = int(splitline[0])
		if strand == '+':
			if ratio != 'N/A':
				peak_plus.update({coord: [float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float(ratio), float(splitline[7])]})
			elif ratio == 'N/A':
				peak_plus.update({coord: [float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float('inf'), float(splitline[7])]})
				div_zero_error_plus.append(coord)
		elif strand == '-':
			if ratio != 'N/A':
				peak_minus.update({coord: [float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float(ratio), float(splitline[7])]})
			elif ratio == 'N/A':
				peak_minus.update({coord: [float(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float('inf'), float(splitline[7])]})
				div_zero_error_minus.append(coord)

#Extract highest peak from ratio clusters
peak_plus = restore_zeros(peak_plus, 'nested')
peak_minus = restore_zeros(peak_minus, 'nested')
local_peak_plus = find_local_peaks(peak_plus)
local_peak_minus = find_local_peaks(peak_minus)
del(peak_plus)
del(peak_minus)
print('All highest peaks extracted:', datetime.datetime.now())

#write gff for peaks and a file with sums in each window
with open(peak_file, 'w') as p:
	p.write('Coordinate\tStrand\tUntreated_us\tUntreated_ds\tTreated_us\tTreated_ds\tReadthrough score (R)\tSignificance score (p_value)\n')
	for coord, val in local_peak_plus.items():
		if coord not in div_zero_error_plus:
			pass
		else:
			val[4] = 'N/A'
		p.write(str(coord) + '\t+\t' + '\t'.join(str(n) for n in val) + '\n')
	for coord, val in local_peak_minus.items():
		if coord not in div_zero_error_minus:
			pass
		else:
			val[4] = 'N/A'
		p.write(str(coord) + '\t-\t' + '\t'.join(str(n) for n in val) + '\n')

print('Finished:', datetime.datetime.now())
print('If N/A is reported as Readthrough score(R) - one of the regions used in the calculation is 0, indicating a really high score;\n if "Significance Score" is reported as 0, then p value from fishers exact test was too low to print.')
