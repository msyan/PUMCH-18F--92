
import numpy as np
import matplotlib.pyplot as plt

import common as c

# =======================

path_weight_dosage = '/path/to/body/weight/and/dosage.txt'
weight_dosage = c.read_weight_dosage(path_weight_dosage)

fac_suvbw = {
    'positive':[],
    'negative':[]
}

for sub in weight_dosage['subjects']:

    print('\nsub: %s , weights: %4.1f kg, dosage: %6.3f mCi' %(sub, weight_dosage[sub]['weight'], weight_dosage[sub]['dosage']))
    print('     %s   1/(dosage / weight) = %6.3f' %(weight_dosage[sub]['positive'], 1.0/(weight_dosage[sub]['dosage']/weight_dosage[sub]['weight'])))

    if weight_dosage[sub]['positive'] == '+':
        fac_suvbw['positive'].append(1.0/(weight_dosage[sub]['dosage']/weight_dosage[sub]['weight']))
    elif weight_dosage[sub]['positive'] == '-':
        fac_suvbw['negative'].append(1.0/(weight_dosage[sub]['dosage']/weight_dosage[sub]['weight']))

print('\nGroup average of 1/(dosage/weight)')
print('  AD: %f +- %f , HC %f +- %f' %(np.mean(fac_suvbw['positive']), np.std(fac_suvbw['positive']),
                                       np.mean(fac_suvbw['negative']), np.std(fac_suvbw['negative'])
                                       ))

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
# '10030':    'Cerebellar Gray Matter'
'10030':    'Cerebellar Cortex'
}

rois = ['10011', '10009', '10000', '10030', '10010', '10014', '10012']

colors = ['steelblue', 'lightsteelblue', 'k', 'grey', 'darkorange', 'yellowgreen', 'forestgreen']

new_times  = [ 10.0*i +   10.0 for i in range(9) ]
new_times += [ 15.0*i +  100.0 for i in range(2) ]
new_times += [ 20.0*i +  130.0 for i in range(3) ]
new_times += [ 30.0*i +  190.0 for i in range(4) ]
new_times += [ 60.0*i +  310.0 for i in range(8) ]
new_times += [180.0*i +  790.0 for i in range(4) ]
new_times += [320.0*i + 1510.0 for i in range(12)]


paths = ['/path/to/TAC/non-AD.txt', '/path/to/TAC/AD.txt']

means = []
for p in paths:

    TACs = c.read_all_tacs(p)
    TACs['ROIs'] = []

    # =======================

    # loop over all rois
    for roi in rois:

        TACs['ROIs'].append(roi)
        TACs[roi] = {
            'time':new_times
        }

        tmp_quants = []

        # loop over subjects
        sum_activities = np.zeros(len(new_times))
        num_subjects = np.zeros(len(new_times))
        for sub in TACs['subjects']:

            # evaluate on new time grid
            tac_func_tmp = c.interpolate(TACs[sub]['time'], TACs[sub][roi])

            tac_tmp = tac_func_tmp(new_times)
            suv_tmp = [ c.compute_SUVbw(t, tac, weight_dosage[sub]['weight'], weight_dosage[sub]['dosage']) for (t, tac) in zip(new_times, tac_tmp)]

            # actually store SUVbw
            TACs[roi][sub] = suv_tmp 

            # make sure no extrapolation happens
            max_time_pt = max(TACs[sub]['time'])
            for (i, t) in enumerate(new_times):
                if t > max_time_pt:
                    TACs[roi][sub][i] = 0
                else:
                    num_subjects[i] += 1

            # accumulate activities
            sum_activities += TACs[roi][sub]

        TACs[roi]['mean'] = [ sum_activities[i] / num_subjects[i] for i in range(len(new_times)) ]

    means.append(TACs)

# =======================

new_times_min = [ t / 60.0 for t in new_times ]

alpbet = ['A', 'B', 'C', 'D', 'D']

width_cm_1 = 14.0
width_inch_1 = width_cm_1 / 2.54

width_cm_2 = 8.0
width_inch_2 = width_cm_2 / 2.54

linewidth = 1.0

markersize = 1.0

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(width_inch_1, width_inch_2), layout='constrained', dpi=300)

handles = []
labels = []
for (i, roi) in enumerate(rois):

    hd, = axs[0].plot(new_times_min, means[0][roi]['mean'], c=colors[i], marker='o', markersize=markersize, linewidth=linewidth)
    axs[1].plot(new_times_min, means[1][roi]['mean'], c=colors[i], marker='o', markersize=markersize, linewidth=linewidth)

    handles.append(hd)
    labels.append(roi_dict[roi])

for (i, ax) in enumerate(axs):

    ax.tick_params(labelsize=8)

    ax.set_xlabel('Time [min]', fontsize=10, font='arial')
    ax.set_ylabel('SUV', fontsize=10, font='arial')
    ax.set_ylim(-0.1, 8.0)
    ax.text(0, 7.5, alpbet[i], fontsize=12, font='helvetica')

plt.savefig('plot_SUVbw_group_mean.png')
plt.show()


# =======================
