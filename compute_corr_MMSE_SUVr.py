
#=======================
#       import
#=======================

import numpy as np
import common as c

import matplotlib.pyplot as plt

#=======================
#
#=======================

#-----------------------
# read in everything

#----------
# MMSE

MMSE = {
	'subjects-1':20,
	'subjects-2':20,
	'subjects-n':20
}


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
# sort

MMSE_vs_SUVR = {
	'time':[],
	'ROIs':[],
	'MMSE':[ MMSE[sub] for sub in subjects_all ]
}

selected_ROIs = ['Global', '10011', '10014', '10012', '10009', '10010']

# loop over all subjects:
for sub in subjects_all:

	# loop over all ROIs
	for roi in selected_ROIs:

		if roi not in MMSE_vs_SUVR['ROIs']:

			MMSE_vs_SUVR['ROIs'].append(roi)
			MMSE_vs_SUVR[roi] = {
				'SUVr': [ [] for i in range(6) ],
				'MMSE': [ [] for i in range(6) ]
			}

		# loop over all time frames of this sub
		for (i, t) in enumerate(results[sub]['SUVr']['time']):

			if t not in MMSE_vs_SUVR['time']:
				MMSE_vs_SUVR['time'].append(t)

			MMSE_vs_SUVR[roi]['SUVr'][i].append(results[sub]['SUVr'][roi][i])
			MMSE_vs_SUVR[roi]['MMSE'][i].append( MMSE[sub] )

for roi in selected_ROIs:

	print()

	for (i, t) in enumerate(MMSE_vs_SUVR['time']):

		res_linregress = c.linear_regression(MMSE_vs_SUVR[roi]['MMSE'][i], MMSE_vs_SUVR[roi]['SUVr'][i])

		str_tmp = '%20s %4s %i %s, r = %7.5f p = %7.5f' %((c.roi_definition(roi)), t, len(MMSE_vs_SUVR[roi]['MMSE'][i]), res_linregress['equation'],res_linregress['rval'], res_linregress['pval'])
		print(str_tmp)

