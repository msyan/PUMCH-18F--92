
import numpy as np
import matplotlib.pyplot as plt

import common as c

# ===================

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
# '10011':	'Frontal',
'10011':    'Prefrontal',
'10012':	'Precuneus',
'10013':	'Parietal',
'10014':	'Temporal' ,
'10015':	'Occipital',
'10016':	'Cuneus',
# '10000':    'White Matter',
'10000':    'Cerebral White Matter',
'10030':    'Cerebellar Gray Matter',

'global':   'Global'
}

alpbet = ['A', 'B', 'C', 'D', 'E', 'E']

def plot_corr(roi, corr, MMSE, SUVr):

    MMSE_arr = np.linspace(1, 29, 6)
    SUVR_arr = corr['slope'] * MMSE_arr + corr['intercept']

    str_tmp = 'SUVR = %5.2f MMSE + %4.2f , r = %4.2f' %(corr['slope'], corr['intercept'], abs(corr['r_val']))

    plt.scatter(MMSE, SUVr)
    plt.plot(MMSE_arr, SUVR_arr)
    plt.text(5, 1.25, str_tmp)

    plt.title(roi)

    plt.ylim(1.0, 2.75)
    plt.xlim(0, 30)

    plt.show()

def plot_corr_all(MMSE_SUVR):

    width_cm_1 = 19.0	# double column for NeuroImage    
    width_inch_1 = width_cm_1 / 2.54

    width_cm_2 = 14.0
    width_inch_2 = width_cm_2 / 2.54


    fontsize_1 = 8
    fontsize_2 = fontsize_1 + 2
    fontsize_3 = fontsize_2 + 2


    rois = ['global', '10011', '10014', '10012', '10009', '10010']

    fig, axs = plt.subplots(layout='constrained', nrows=2, ncols=3, figsize=(width_inch_1, width_inch_2), dpi=300)

    MMSE_arr = np.linspace(1, 29, 6)

    for (i, roi) in enumerate(rois):

        k = int(i%3)
        j = int( (i - i%3)/3 )
        ax = axs[j][k]

        correlations = c.linear_regress(MMSE_SUVR['MMSE'], MMSE_SUVR[roi])

        SUVR_arr = correlations['slope'] * MMSE_arr + correlations['intercept']

        str_tmp = 'SUVR = %5.2f MMSE + %4.2f\nr = %5.2f' %(
            correlations['slope'], 
            correlations['intercept'], 
            correlations['r_val']
            )

        ax.scatter(MMSE_SUVR['MMSE'], MMSE_SUVR[roi], color='steelblue', s=20)
        ax.plot(MMSE_arr, SUVR_arr, color='k', linestyle=':')
        ax.text(4, 1.16, str_tmp, fontsize=fontsize_1)
        ax.set_title(roi_dict[roi], fontsize=fontsize_3, font='arial')

        ax.set_ylabel('SUVR', fontsize=fontsize_2, font='arial')
        ax.set_ylim(1.0, 2.75)
        ax.set_xlabel('MMSE', fontsize=fontsize_2, font='arial')
        ax.set_xlim(0, 30)
        ax.tick_params(labelsize=fontsize_1)

        ax.set_yticks(np.linspace(1.2, 2.8, 5))

        if roi == '10010':
            for (MMSE_tmp, SUVR_tmp) in zip(MMSE_SUVR['MMSE'], MMSE_SUVR[roi]):
                print(MMSE_tmp, SUVR_tmp)


    plt.savefig('plot_MMSE_vs_SUVR')
    plt.show()

# ===================

path = '/path/to/MMSE/results.txt'

MMSE_SUVR = c.read_file(path)
c.compute_suvr_global(MMSE_SUVR)

plot_corr_all(MMSE_SUVR)
