
#=======================
#       import
#=======================

import numpy as np
import common as c

#=======================
#
#=======================

def compute_global_tac_suvr(results):

	global_rois = ['10009', '10010', '10011', '10012', '10014']

	# loop over subjects
	for sub in results['subjects']:

		tac_reference = results[sub]['TACs']['CBGM']

		tac_global = [ [] for i in range(len(results[sub]['TACs']['time']))]

		# loop over global rois
		# sort and store
		for roi in global_rois:

			# loop over all time points
			for (tac_glb, tac_roi) in zip(tac_global, results[sub]['TACs'][roi]):
				tac_glb.append(tac_roi)

		results[sub]['TACs']['global'] = [ np.mean(tac_glb) for tac_glb in tac_global ]
		results[sub]['SUVr']['global'] = [ tac_glb/tac_ref for (tac_glb, tac_ref) in zip(results[sub]['TACs']['global'], tac_reference) ]


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
compute_global_tac_suvr(results)

#-----------------------



suvrs_AD = [ [] for i in range(6) ]
suvrs_HC = [ [] for i in range(6) ]

rois = ['global', '10009', '10010', '10011', '10012', '10014']

for sub in results['subjects']:

	for roi in rois:

		suvr_tmp = results[sub]['SUVr'][roi]

		# if AD
		if positivity[sub] == 0:

			for (i, suvr) in enumerate(suvr_tmp):

				suvrs_AD[i].append(suvr)

		elif positivity[sub] == 1:

			for (i, suvr) in enumerate(suvr_tmp):

				suvrs_HC[i].append(suvr)

	print(sub, positivity[sub])


max_AD = [ max(suvrs) for suvrs in suvrs_AD ]
min_AD = [ min(suvrs) for suvrs in suvrs_AD ]
max_HC = [ max(suvrs) for suvrs in suvrs_HC ]
min_HC = [ min(suvrs) for suvrs in suvrs_HC ]

print('\n\n ROIs:', rois)
print('max  min')
for (maxx, minn) in zip(max_HC, min_AD):
	print('%4.2f %4.2f' %(maxx, minn))



