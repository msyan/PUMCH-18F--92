
#=======================
#       import
#=======================

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
# AD, HC

BPND_AD, BPND_HC = c.prepare_BPND_AD_HC_global_mean(BPNDs, positivity)


#----------
# plot

rois = ['Global', '10011', '10014', '10012', '10009', '10010']

quants_AD = [BPND_AD[roi] for roi in rois]
quants_HC = [BPND_HC[roi] for roi in rois]

c.plot_box_BPND_AD_HC_regional(quants_AD, quants_HC, rois)

