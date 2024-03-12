
import numpy as np
import matplotlib.pyplot as plt

import common as c

# =======================

roi_dict = {
'10001':	'Whole Cerebellum',
'10002':	'Cerebral White Matter',
'10003':	'Thalamus',
'10004':	'Caudate',
'10005':	'Putamen',
'10006':	'Brain Stem',
'10007':	'Hippocampus',
'10008':	'Para Hippocampal',
'10009':	'Anterior Cingulate',
'10010':	'Posterior Cingulate',
'10011':	'Prefrontal',
'10012':	'Precuneus',
'10013':	'Parietal',
'10014':	'Temporal' ,
'10015':	'Occipital',
'10016':	'Cuneus',
# '10000':    'White Matter',
'10000':    'Cerebral White Matter',
'10030':    'Cerebellar Gray Matter',

'global':	'Global'
}


rois = ['10009', '10010', '10011', '10012', '10014']
rois_glb = ['global', '10011', '10014', '10012', '10009', '10010']



num_rois = len(rois)
roi_ref = '10030'

new_times  = [ 100.0*i +    5.0 for i in range(17) ]
new_times += [ 500.0*i + 1705.0 for i in range(1)]
new_times += [ 700.0*i + 2205.0 for i in range(2)]
new_times += [ 300.0*i + 3605.0 for i in range(2)]
new_times += [  50.0*i + 3905.0 for i in range(8)]
new_times += [ 200.0*i + 4405.0 for i in range(3)]

# =======================

def compute_suvr_mean_extreme(paths):

	rtn_stat = []
	rtn_suvr = []

	for (ip, p) in enumerate(paths):

		stats = {
			'ROIs':[]
		}

		TACs = c.read_all_tacs(p)
		TACs['ROIs'] = []

		SUVr = {
			'subjects':TACs['subjects'],
			'ROIs': rois
		}

		# compute SUVr
		# loop over subjects
		for sub in SUVr['subjects']:

			SUVr[sub] = {
				'time':TACs[sub]['time'][3:-1],
			}

			# loop over rois
			for roi in SUVr['ROIs']:

				# compute suvr for all subs and cortical rois
				suvr_tmp = []
				for (tac_roi, tac_ref) in zip(TACs[sub][roi], TACs[sub][roi_ref]):
					if tac_ref == 0.0:
						suvr_tmp.append(0.0)
					else:
						suvr_tmp.append(tac_roi / tac_ref)

				suvr_tmp = suvr_tmp[3:-1]

				# interpolate for a uniform array of time
				suvr_func = c.interpolate(SUVr[sub]['time'], suvr_tmp)
				suvr_intp = [ suvr_func(t) for t in new_times ]

				# make sure no extrapolation happens
				max_time = max(SUVr[sub]['time'])
				for (i, t) in enumerate(new_times):
					if t > max_time:
						suvr_intp[i] = 0

				SUVr[sub][roi] = suvr_intp

			# compute averaged suvr for global cortical area
			suvr_glb_tmp = []
			for (i, t) in enumerate(new_times):

				glb_tmp = [ SUVr[sub][roi][i] for roi in SUVr['ROIs'] ]
				suvr_glb_tmp.append( np.mean(glb_tmp) )

			SUVr[sub]['global'] = suvr_glb_tmp

		# loop over rois to find their min, max, mean
		for roi in rois_glb:

			# initiate roi for stats
			if roi not in stats['ROIs']:

				stats['ROIs'].append(roi)
				stats[roi] = {
					'mean':[],
					'max':[],
					'min':[]
				}

			# loop over time
			for (i, t) in enumerate(new_times):

				# sort suvr statistics by loop over all subjects
				suvr_stats = []
				for sub in SUVr['subjects']:

					if SUVr[sub][roi][i] != 0:
						suvr_stats.append( SUVr[sub][roi][i] )

				stats[roi]['mean'].append( np.mean(suvr_stats) )
				stats[roi]['min'].append( min(suvr_stats) )
				stats[roi]['max'].append( max(suvr_stats) )

		rtn_stat.append(stats)
		rtn_suvr.append(SUVr)

	return rtn_stat, rtn_suvr

# =======================

paths = ['/path/to/TAC/non-AD.txt', '/path/to/TAC/AD.txt']
ADs = ['HC', 'AD']


stats, suvrs = compute_suvr_mean_extreme(paths)

# =======================

colors = ['steelblue', 'darkorange']

i_1st_pt = 1

new_times_min = [ t / 60.0 for t in new_times ]

alpbet = ['A', 'B', 'C', 'D', 'E', 'F']

width_cm_1 = 19.0
width_inch_1 = width_cm_1 / 2.54

width_cm_2 = 12.0
width_inch_2 = width_cm_2 / 2.54

fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(width_inch_1, width_inch_2), layout='constrained', dpi=300)

for (i, roi) in enumerate(rois_glb):

	j = int( ( i - (i%3) ) / 2 )
	k = int(i%3)

	# loop over HC and AD
	# labels = ADs
	handles = []
	for a in range(2):

		hd_plot, = axss[j][k].plot(new_times_min[i_1st_pt:-1], stats[a][roi]['mean'][i_1st_pt:-1], color=colors[a], linestyle='-')
		hd_fill = axss[j][k].fill_between(
			new_times_min[i_1st_pt:-1],
			stats[a][roi]['min'][i_1st_pt:-1],
			stats[a][roi]['max'][i_1st_pt:-1],
			facecolor = 'gray',
			alpha=0.3
			)

		hd_up  = axss[j][k].plot(new_times_min[i_1st_pt:-1], stats[a][roi]['max'][i_1st_pt:-1], color=colors[a], linestyle=':')
		hd_low = axss[j][k].plot(new_times_min[i_1st_pt:-1], stats[a][roi]['min'][i_1st_pt:-1], color=colors[a], linestyle=':')

		handles.append((hd_plot, hd_fill))

	axss[j][k].tick_params(labelsize=8)
	axss[j][k].set_xlabel('Time [min]', fontsize=10, font='arial')
	axss[j][k].set_ylabel('SUVR', fontsize=10, font='arial')
	axss[j][k].set_ylim(0.8, 3.0)
	axss[j][k].set_title(roi_dict[roi], fontsize=12, font='arial')

	axss[j][k].legend(labels=ADs, handles=handles, fontsize=8)

plt.savefig('plot_suvr_dynamic_glb_rgn.png')
plt.show()



