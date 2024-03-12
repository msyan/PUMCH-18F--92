
#=======================
#       import
#=======================

import numpy as np
import common as c

#=======================
#
#=======================

#-----------------------
# read in everything

#----------
# BPNDs

subjects = ['subjects-1', 'subjects-2', 'subjects-3', 'subjects-n']

path_file_BPND = '/Path/to/BPnd/results.txt'

paths_BPND = {
	'subjects':subjects,
	'main_path':path_file_BPND
}

path_positivity = '/Path/to/positivity/results.txt'

BPNDs = c.read_in_BPNDs(paths_BPND)
positivity = c.read_in_BPNDs_from_old_result(path_positivity, BPNDs=BPNDs)


#----------
# TAC, SUVr

# all subjects
subjects_all = ['subjects-1', 'subjects-2', 'subjects-3', 'subjects-n']
path_file_TAC = '/Path/to/TAC/results/'

paths_TACs = {
	'subjects':subjects_all,
	'main_path':path_file_TAC
}

TACs = c.read_TAC_4_all(paths_TACs)

# reference roi for all subjects
path_file_ref = '/Path/to/reference/results/'

paths_ref = {
	'subjects':subjects_all,
	'main_path':path_file_ref
}

TACs_ref = c.read_TAC_4_all_reference_roi(paths_ref)

results = c.compute_SUVrs_4_all(TACs, TACs_ref)
c.compute_suvr_global(results)

#----------
# all, AD, HC

SUVr_vs_BPND_all, SUVr_vs_BPND_AD, SUVr_vs_BPND_HC = c.prepare_BPND_SUVr_AD_HC(results, BPNDs, positivity)

#-----------------------
# compute linear regression

results_lg = {
	'times':[]
}

for t in SUVr_vs_BPND_all['time']:
	
	print('\n\nTime: %4s min' %(t))	

	res_linregress = c.linear_regression(SUVr_vs_BPND_all[t]['all_rois']['DVRs'], SUVr_vs_BPND_all[t]['all_rois']['SUVr'])

	print('  %s' %(res_linregress['equation']))
	print('  r = %7.5f p = %7.5f' %(res_linregress['rval'], res_linregress['pval']))

	results_lg['times'].append(t)
	results_lg[t] = res_linregress

c.plot_DVR_VS_SUVr_together(results_lg, SUVr_vs_BPND_AD, SUVr_vs_BPND_HC)
