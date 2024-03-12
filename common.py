
#=======================
#       import
#=======================

import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#=======================
#       define
#=======================

#-----------------------
# linear regression

def linear_regress(x, y):

    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)

    rtn = {
        'slope':slope,
        'intercept':intercept,
        'r_val':r_value,
        'p_val':p_value
    }

    return rtn

#-----------------------
# time definition

time_dict = {
	'7.5' :  '0 - 15 min',
	'22.5': '15 - 30 min',
	'37.5': '30 - 45 min',
	'52.5': '45 - 60 min',
	'67.5': '60 - 75 min',
	'82.5': '75 - 90 min'
}

#-----------------------
# roi definition

roi_dict = {
	'ROIs': [str(10001+i) for i in range(16)],
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
	'10011':	'Prefrontal',
	'10012':	'Precuneus',
	'10013':	'Parietal',
	'10014':	'Temporal' ,
	'10015':	'Occipital',
	'10016':	'Cuneus',
}

def roi_definition(roi_index):
	if roi_index not in roi_dict['ROIs']:
		return roi_index
	else:
		return roi_dict[roi_index]

#-----------------------
# manage TAC, SUVr

def read_in(path, off_sets:list, skiprows=8):
	rtn = np.loadtxt(path, skiprows=skiprows, usecols=(15+off_sets[0], 16+off_sets[1], 17+off_sets[2]), dtype=np.float64)
	return rtn

def sort_TAC(TAC_res, TAC_factor):
	
	TACs = {
		'ROIs':[],
		'time':[]
	}

	for (i, line) in enumerate(TAC_res):

		roi_tmp  = str(int(line[0]))
		time_tmp = line[1]
		tac_tmp  = line[2]

		# intiate an ROI
		if roi_tmp not in TACs['ROIs']:
			TACs['ROIs'].append(roi_tmp)
			TACs[roi_tmp] = []

		TACs[roi_tmp].append( line[2]*TAC_factor )
		TACs['time'].append(line[1])

	npt_time = int( len(TAC_res)/len(TACs['ROIs']) )
	TACs['time'] = TACs['time'][0:npt_time]

	for roi in TACs['ROIs']:

		if True in np.isnan(TACs[roi]):
			print('  Warning: Nan exist in roi %s' %(roi))

		TACs[roi] = np.nan_to_num(TACs[roi])

	return TACs

def compute_SUVr(TAC_res, roi_reference):

	SUVr = {
		'ROIs': TAC_res['ROIs'],
		'time': TAC_res['time']
	}

	npt_time = TAC_res[roi_reference]

	for roi in TAC_res['ROIs']:
		SUVr[roi] = [ TAC_res[roi][i]/TAC_res[roi_reference][i] for i in range(len(npt_time)) ]

	return SUVr

def compute_SUVr_4_all(paths, roi_reference):

	results = {
		'subjects':paths['subjects']
	}

	for sub in paths['subjects']:

		print('computing for %s' %(sub))

		path_tmp = os.path.join(paths['main_path'], sub+'.txt')

		if sub == 'TSAD-09':
			off_sets = [4, 4, 4]
		elif (sub == 'TSAD-53') or (sub == 'TSAD-61'):
			off_sets = [5, 5, 5]
		elif sub == 'TSAD-26':
			off_sets = [0, 1, 1]
		else:
			off_sets = [0, 0, 0]

		if (sub == 'TSAD-53') or (sub == 'TSAD-61') or (sub == 'TSAD-63'):
			TAC_factor = 1.0
		else:
			TAC_factor = 1000.0

		TACs = sort_TAC( read_in(path_tmp, off_sets), TAC_factor )


		results[sub] = {
			'TACs': TACs,
			'SUVr': compute_SUVr(TACs, roi_reference)
		}

	return results

def read_TAC_4_all(paths):

	results = {
		'subjects':paths['subjects']
	}

	for sub in paths['subjects']:

		print('computing for %s' %(sub))

		path_tmp = os.path.join(paths['main_path'], sub+'.txt')

		if sub == 'TSAD-09':
			off_sets = [4, 4, 4]
		elif (sub == 'TSAD-53') or (sub == 'TSAD-61'):
			off_sets = [5, 5, 5]
		elif sub == 'TSAD-26':
			off_sets = [0, 1, 1]
		else:
			off_sets = [0, 0, 0]

		if (sub == 'TSAD-53') or (sub == 'TSAD-61') or (sub == 'TSAD-63'):
			TAC_factor = 1.0
		else:
			TAC_factor = 1000.0

		TACs = sort_TAC( read_in(path_tmp, off_sets), TAC_factor )


		results[sub] = {
			'TACs': TACs,
		}

	return results

def read_TAC_4_all_reference_roi(paths):

	def _read_in_ref(path, skiprows, usecols:list):

		f = open(path, 'r')

		lines = f.readlines()

		tmp_results = []

		for (i, line) in enumerate(lines):

			str_tmp = line.split()

			if i < skiprows:
				pass
			elif len(str_tmp)<1:
				pass
			else:
				tmp_results.append( [str_tmp[c] for c in usecols] )

		f.close()

		return tmp_results

	def _sort_TAC(TAC_result, TAC_factor):

		TACs = {
			'ROIs':[],
			'time':[]
			}

		for res in TAC_result:

			if res[0] not in TACs['ROIs']:
				TACs['ROIs'].append(res[0])
				TACs[res[0]] = []

			if res[1] not in TACs['time']:
				TACs['time'].append(res[1])

			TACs[res[0]].append( float(res[2])*TAC_factor )

		return TACs

	results = {
		'subjects':paths['subjects']
	}

	for sub in paths['subjects']:

		print('computing for %s' %(sub))

		path_tmp = os.path.join(paths['main_path'], sub+'.txt')

		TAC_res = _read_in_ref(path_tmp, skiprows=8, usecols=[-4, -3, -2])

		if (sub == 'TSAD-63') or (sub == 'TSAD-11') or (sub == 'TSAD-09'):
			TAC_factor = 1.0
		else:
			TAC_factor = 1000.0

		TACs = _sort_TAC(TAC_res, TAC_factor)

		results[sub] = {
			'TACs': TACs,
		}

	return results

def compute_SUVrs_4_all(TACs, TACs_ref):

	results = {
		'subjects':TACs['subjects']
	}

	for sub in results['subjects']:

		results[sub] = {
			'TACs':TACs[sub]['TACs']
		}

		# append reference TAC
		results[sub]['TACs']['ROIs'].append( TACs_ref[sub]['TACs']['ROIs'][0] )
		results[sub]['TACs'][ TACs_ref[sub]['TACs']['ROIs'][0] ] = TACs_ref[sub]['TACs'][ TACs_ref[sub]['TACs']['ROIs'][0] ]

		# compute SUVr
		SUVr = {
			'ROIs': results[sub]['TACs']['ROIs'],
			'time': results[sub]['TACs']['time']
		}
		for roi in SUVr['ROIs']:

			SUVr[roi] = []

			for (tac_val, tac_ref) in zip( results[sub]['TACs'][roi], TACs_ref[sub]['TACs'][TACs_ref[sub]['TACs']['ROIs'][0]] ):
				SUVr[roi].append( tac_val/tac_ref )

		results[sub]['SUVr'] = SUVr

	return results

def stat_single_quant(results, quant):

	SUBs = results['subjects']
	ROIs = results[SUBs[0]][quant]['ROIs']
	nSUB = len(SUBs)

	statistics = {
		'ROIs':ROIs,
		'time':[]
	}

	# loop over all rois
	for r in ROIs:

		tmp_dict = {
			'time':[]
		}

		# loop over all subjects
		for s in SUBs:

			# loop over all time
			for (i, t) in enumerate(results[s][quant]['time']):

				# check if t is already included
				if str(t) not in statistics['time']:
					
					statistics['time'].append(str(t))
					statistics[str(t)] = {'time':str(t)}

					for rr in ROIs:
						statistics[str(t)][rr] = []

				# check if t is already included
				if str(t) not in tmp_dict['time']:

					tmp_dict['time'].append(str(t))
					tmp_dict[str(t)] = []

				tmp_dict[str(t)].append( results[s][quant][r][i] )

		for t in tmp_dict['time']:

			statistics[t][r].append( np.mean(tmp_dict[t]) )
			statistics[t][r].append( np.std(tmp_dict[t]) )
			statistics[t][r].append( len(tmp_dict[t]) )

	return statistics

def average_single_quant_old(results, quant):

	SUBs = results['subjects']
	ROIs = results[SUBs[0]]['TACs']['ROIs']
	nSUB = len(SUBs)

	mean_values = {
		'ROIs':ROIs
	}

	for roi in ROIs:

		tmp_dict = {
			'n_time':[],
			'quants':[]
		}

		for sub in SUBs:
			tmp_dict['quants'].append(results[sub][quant][roi])
			tmp_dict['n_time'].append(len(results[sub][quant][roi]))

		min_t_pt = min(tmp_dict['n_time'])
		acc_val = []
		for t in range(min_t_pt):

			tmp_val = 0.0
			for qs in tmp_dict['quants']:
				tmp_val += qs[t]
			acc_val.append(tmp_val)

		mean_values[roi] = {
			'time':results[SUBs[0]]['TACs']['time'][0:min_t_pt],
			'vals':[acvl/nSUB for acvl in acc_val]
		}

	return mean_values

def read_file(path):

    f = open(path, 'r')
    lines = f.readlines()

    MMSE_SUVR = {
        'MMSE':[]
    }

    for (i, line) in enumerate(lines):

        str_tmp = line.split()

        # on 1st line, read roi
        if i == 0:

            # MMSE_SUVR['ROIs'] = str_tmp[1:7]
            MMSE_SUVR['ROIs'] = str_tmp[1:8]

            for roi in MMSE_SUVR['ROIs']:
                MMSE_SUVR[roi] = []

        else:

            # read in MMSE score
            MMSE_SUVR['MMSE'].append( float(str_tmp[0]) )

            # read in SUVR in each roi
            for (j, roi) in enumerate(MMSE_SUVR['ROIs']):

                MMSE_SUVR[roi].append( float(str_tmp[j+1]) )

    return MMSE_SUVR

def compute_suvr_global(MMSE_SUVR):

    num_AD = len(MMSE_SUVR['MMSE'])

    MMSE_SUVR['global'] = []

    # loop over subjects
    for i in range(num_AD):

        # loop over roi
        tmp_suvr_glb = []
        for roi in ['10009', '10010', '10011', '10012', '10014']:

            tmp_suvr_glb.append(MMSE_SUVR[roi][i]) 
            print(roi)

        MMSE_SUVR['global'].append( np.mean(tmp_suvr_glb) )

def read_all_tacs(path):

    # lines = np.loadtxt(path, skiprows=1)
    f = open(path, 'r')
    lines = f.readlines()

    TACs = {
        'subjects':[],
    }

    for (i, line) in enumerate(lines):

        if i > 0:

            str_tmp = line.split()

            roi = str_tmp[1]
            time = float(str_tmp[2])
            tac = float(str_tmp[3])
            sub = str_tmp[5]

            if sub not in TACs['subjects']:
               TACs['subjects'].append(sub)
               TACs[sub] = {
                    'time':[],
                    'ROIs':[]
              }

            if roi not in TACs[sub]['ROIs']:
                TACs[sub]['ROIs'].append(roi)
                TACs[sub][roi] = []

            if time not in TACs[sub]['time']:
                TACs[sub]['time'].append( time )

            TACs[sub][roi].append(tac)

    f.close()

    return TACs

def interpolate(time, activity):
    tac_func = sc.interpolate.CubicSpline(time, activity)
    # tac_func = sc.interpolate.make_interp_spline(time, activity)
    return tac_func

def read_weight_dosage(path):

    f = open(path, 'r')

    lines = f.readlines()

    rtn = {
        'subjects':[]
    }

    for line in lines:

        str_tmp = line.split()

        sub_tmp = str_tmp[0]
        weight_tmp = float(str_tmp[5])  # unit in kg
        dosage_tmp = float(str_tmp[-1]) # unit in mCi
        positive = str_tmp[6]

        # if sub_tmp not in rtn['subjects']:
        rtn['subjects'].append(sub_tmp)
        rtn[sub_tmp] = {
            'weight': weight_tmp,
            'dosage': dosage_tmp,
            'positive': positive
        }

    return rtn

def compute_SUVbw(time, tac, weight, dosage):

    # time      in unit of second
    # tac       in unit of kBq/ml
    # weight    in unit of kg
    # dosage    in unit of mCi

    # use 1 g/ml as brain's density to convert weight into volumn
    density_brain = 1.0 
    weight_norm = (weight * 1000.0 ) / density_brain # now unit is in ml

    # convert dosage from mCi into kBq
    # 1 Ci = 3.7 * 10^10 Bq
    fac_Ci_2_Bq = 3.7E10
    dosage_kBq = ((dosage / 1000.0) * fac_Ci_2_Bq) / 1000.0 

    SUVbw = tac / ( dosage_kBq / weight_norm )
    
    return SUVbw

#-----------------------
# manage BPND

def read_in_BPNDs(paths):

	results = {
		'subjects':paths['subjects']		
	}

	for sub in paths['subjects']:

		results[sub] = {
			'ROIs':[]
		}

		path_tmp = os.path.join(paths['main_path'], sub+'_BPND.txt')

		if sub == 'TSAD-63':
			off_sets = [-15, -15, -16]
		else:
			off_sets = [-15, -15, -15]

		BPND_tmp = read_in(path_tmp, off_sets, skiprows=1)

		for line in BPND_tmp:

			results[sub]['ROIs'].append( str(int(line[0])) )
			results[sub][ str(int(line[0])) ] = line[2]

	return results

def read_in_BPNDs_from_old_result(path, BPNDs=[]):

	def _read_in(path, usecols:list, skiprows=0):

		f = open(path, 'r')
		lines = f.readlines()

		results = []

		for (i, line) in enumerate(lines):

			if i < skiprows:
				pass
			else:
				str_tmp = line.split()
				results.append([str_tmp[c] for c in usecols])

		return results

	#res_BPND = np.loadtxt(path, skiprows=1, usecols=(0, 1, 3, 4), dtype=np.float64)
	res_BPND = _read_in(path, skiprows=1, usecols=[0, 1, 3, 4])


	if BPNDs == []:
		BPNDs = {
			'subjects':[]
		}

	# store the positiveness of subjects
	positiveness = {
		'subjects':[]
	}

	for line in res_BPND:

		sub_tmp = line[0]
		sub = '%s-%s' %(sub_tmp[0:4], sub_tmp[4:6])

		# check if this subject is included
		if sub not in BPNDs['subjects']:
			BPNDs['subjects'].append(sub)
			BPNDs[sub] = {
				'ROIs':[]
			}

		roi = line[1]

		# check if this roi is included for this subject
		if roi not in BPNDs[sub]['ROIs']:
			BPNDs[sub]['ROIs'].append(roi)
			BPNDs[sub][roi] = float(line[3])

		# read in positiveness
		if sub not in positiveness['subjects']:
			positiveness['subjects'].append(sub)
			positiveness[sub] = int(line[2])

		# 63: positive
		# 53, 51: negative
		subs = [53, 61, 63]
		for s in subs:
			positiveness['subjects'].append('TSAD-%2i' %(s))

			# 0 is AD
			if s == 63:
				positiveness['TSAD-%2i' %(s)] = 0
			else:
				positiveness['TSAD-%2i' %(s)] = 1

	return positiveness

def prepare_BPND_SUVr(results_SUVr, results_BPND):

	pairs_SUVr_BPND = {
		'time':[]
	}

	# loop over all subjects
	for sub in results_SUVr['subjects']:

		# loop over all times
		for t in results_SUVr[sub]['SUVr']['time']:

			# check if this time is already included
			if str(t) not in pairs_SUVr_BPND['time']:
				pairs_SUVr_BPND['time'].append(str(t))
				pairs_SUVr_BPND[str(t)] = {
					'all_rois':{
						'SUVr':[],
						'DVRs':[],
						'BPND':[]
						},
					'ROIs':[]
				}

			# pick the index of time piont t
			index_t = results_SUVr[sub]['SUVr']['time'].index(t)

			# loop over all rois
			for roi in results_SUVr[sub]['SUVr']['ROIs']:

				# check if this roi is already included
				if roi not in pairs_SUVr_BPND[str(t)]['ROIs']:
					pairs_SUVr_BPND[str(t)]['ROIs'].append(roi)
					pairs_SUVr_BPND[str(t)][roi] = {
						'SUVr':[],
						'DVRs':[],
						'BPND':[]
					}

				# append a new pair
				pairs_SUVr_BPND[str(t)]['all_rois']['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
				pairs_SUVr_BPND[str(t)]['all_rois']['BPND'].append( results_BPND[sub][roi])
				pairs_SUVr_BPND[str(t)]['all_rois']['DVRs'].append( results_BPND[sub][roi] + 1.0 )

				pairs_SUVr_BPND[str(t)][roi]['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
				pairs_SUVr_BPND[str(t)][roi]['BPND'].append( results_BPND[sub][roi] )
				pairs_SUVr_BPND[str(t)][roi]['DVRs'].append( results_BPND[sub][roi] + 1.0 )

	return pairs_SUVr_BPND

def prepare_BPND_SUVr_AD_HC(results_SUVr, results_BPND, positivity):


	selected_ROIs = ['10009', '10010', '10011', '10014', '10013', '10012']
	# selected_ROIs = ['10002']


	pairs_SUVr_BPND_all = {
		'time':[]
	}

	pairs_SUVr_BPND_AD = {
		'time':[]
	}

	pairs_SUVr_BPND_HC = {
		'time':[]
	}

	# loop over all subjects
	for sub in results_SUVr['subjects']:

		# print('  %s , %i' %(sub, positivity[sub]))

		# loop over all times
		for t in results_SUVr[sub]['SUVr']['time']:

			# pick the index of time piont t
			index_t = results_SUVr[sub]['SUVr']['time'].index(t)

			# positive
			if positivity[sub] == 0:

				# check if this time is already included
				if str(t) not in pairs_SUVr_BPND_AD['time']:
					pairs_SUVr_BPND_AD['time'].append(str(t))
					pairs_SUVr_BPND_AD[str(t)] = {
						'all_rois':{
							'SUVr':[],
							'DVRs':[],
							'BPND':[]
							},
						'ROIs':[]
					}


				# loop over rois
				# for roi in results_BPND[sub]['ROIs']:
				for roi in selected_ROIs:

					# print(sub, roi)


					# check if this roi is already included
					if roi not in pairs_SUVr_BPND_AD[str(t)]['ROIs']:
						pairs_SUVr_BPND_AD[str(t)]['ROIs'].append(roi)
						pairs_SUVr_BPND_AD[str(t)][roi] = {
							'SUVr':[],
							'DVRs':[],
							'BPND':[]
						}

					# append a new pair
					pairs_SUVr_BPND_AD[str(t)]['all_rois']['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
					pairs_SUVr_BPND_AD[str(t)]['all_rois']['BPND'].append( results_BPND[sub][roi])
					pairs_SUVr_BPND_AD[str(t)]['all_rois']['DVRs'].append( results_BPND[sub][roi] + 1.0 )

					pairs_SUVr_BPND_AD[str(t)][roi]['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
					pairs_SUVr_BPND_AD[str(t)][roi]['BPND'].append( results_BPND[sub][roi] )
					pairs_SUVr_BPND_AD[str(t)][roi]['DVRs'].append( results_BPND[sub][roi] + 1.0 )

			# negative
			elif positivity[sub] == 1:

				# check if this time is already included
				if str(t) not in pairs_SUVr_BPND_HC['time']:
					pairs_SUVr_BPND_HC['time'].append(str(t))
					pairs_SUVr_BPND_HC[str(t)] = {
						'all_rois':{
							'SUVr':[],
							'DVRs':[],
							'BPND':[]
							},
						'ROIs':[]
					}

				# loop over all rois
				# for roi in results_BPND[sub]['ROIs']:
				for roi in selected_ROIs:

					# check if this roi is already included
					if roi not in pairs_SUVr_BPND_HC[str(t)]['ROIs']:
						pairs_SUVr_BPND_HC[str(t)]['ROIs'].append(roi)
						pairs_SUVr_BPND_HC[str(t)][roi] = {
							'SUVr':[],
							'DVRs':[],
							'BPND':[]
						}

					# append a new pair
					pairs_SUVr_BPND_HC[str(t)]['all_rois']['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
					pairs_SUVr_BPND_HC[str(t)]['all_rois']['BPND'].append( results_BPND[sub][roi])
					pairs_SUVr_BPND_HC[str(t)]['all_rois']['DVRs'].append( results_BPND[sub][roi] + 1.0 )

					pairs_SUVr_BPND_HC[str(t)][roi]['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
					pairs_SUVr_BPND_HC[str(t)][roi]['BPND'].append( results_BPND[sub][roi] )
					pairs_SUVr_BPND_HC[str(t)][roi]['DVRs'].append( results_BPND[sub][roi] + 1.0 )

			# always include to all
			# check if this time is already included
			if str(t) not in pairs_SUVr_BPND_all['time']:
				pairs_SUVr_BPND_all['time'].append(str(t))
				pairs_SUVr_BPND_all[str(t)] = {
					'all_rois':{
						'SUVr':[],
						'DVRs':[],
						'BPND':[]
						},
					'ROIs':[]
				}


			# loop over all rois
			# for roi in results_BPND[sub]['ROIs']:
			for roi in selected_ROIs:

				# check if this roi is already included
				if roi not in pairs_SUVr_BPND_all[str(t)]['ROIs']:
					pairs_SUVr_BPND_all[str(t)]['ROIs'].append(roi)
					pairs_SUVr_BPND_all[str(t)][roi] = {
						'SUVr':[],
						'DVRs':[],
						'BPND':[]
					}

				# append a new pair
				pairs_SUVr_BPND_all[str(t)]['all_rois']['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
				pairs_SUVr_BPND_all[str(t)]['all_rois']['BPND'].append( results_BPND[sub][roi])
				pairs_SUVr_BPND_all[str(t)]['all_rois']['DVRs'].append( results_BPND[sub][roi] + 1.0 )

				pairs_SUVr_BPND_all[str(t)][roi]['SUVr'].append( results_SUVr[sub]['SUVr'][roi][index_t] )
				pairs_SUVr_BPND_all[str(t)][roi]['BPND'].append( results_BPND[sub][roi] )
				pairs_SUVr_BPND_all[str(t)][roi]['DVRs'].append( results_BPND[sub][roi] + 1.0 )

	return pairs_SUVr_BPND_all, pairs_SUVr_BPND_AD, pairs_SUVr_BPND_HC

def prepare_BPND_AD_HC(results_BPND, positivity):

	# cortex_global_rois = ['10011', '10014', '10009', '10010', '10012', '10013']
	cortex_global_rois = ['10011', '10014', '10009', '10010', '10012']

	BPND_AD = {
		'subjects': [],
		'ROIs': ['Global'],
		'Global': []
	}

	BPND_HC = {
		'subjects':[],
		'ROIs': ['Global'],
		'Global': []
	}

	for roi in results_BPND[results_BPND['subjects'][0]]['ROIs']:

		BPND_AD['ROIs'].append( roi )
		BPND_AD[roi] = []

		BPND_HC['ROIs'].append( roi )
		BPND_HC[roi] = []


	for sub in results_BPND:

		# 0 is AD
		if positivity[sub] == 0:

			BPND_AD['subjects'].append(sub)
			BPND_AD[sub] = results_BPND[sub]

			for roi in BPND_AD[sub]['ROIs']:

				if roi in cortex_global_rois:
					BPND_AD['Global'].append(BPND_AD[sub][roi])
	
				BPND_AD[roi].append(BPND_AD[sub][roi])


		elif positivity[sub] == 1:

			BPND_HC['subjects'].append(sub)
			BPND_HC[sub] = results_BPND[sub]

			for roi in BPND_HC[sub]['ROIs']:

				if roi in cortex_global_rois:
					BPND_HC['Global'].append(BPND_HC[sub][roi])
	
				BPND_HC[roi].append(BPND_HC[sub][roi])

	return BPND_AD, BPND_HC

def prepare_BPND_AD_HC_global_mean(results_BPND, positivity):

	# cortex_global_rois = ['10011', '10014', '10009', '10010', '10012', '10013']
	cortex_global_rois = ['10011', '10014', '10009', '10010', '10012']

	BPND_AD = {
		'subjects': [],
		'ROIs': ['Global'],
		'Global': []
	}

	BPND_HC = {
		'subjects':[],
		'ROIs': ['Global'],
		'Global': []
	}

	for roi in results_BPND[results_BPND['subjects'][0]]['ROIs']:

		BPND_AD['ROIs'].append( roi )
		BPND_AD[roi] = []

		BPND_HC['ROIs'].append( roi )
		BPND_HC[roi] = []


	for sub in results_BPND:

		# 0 is AD
		if positivity[sub] == 0:

			BPND_AD['subjects'].append(sub)
			BPND_AD[sub] = results_BPND[sub]

			BPnd_glb = []
			for roi in BPND_AD[sub]['ROIs']:

				if roi in cortex_global_rois:
					# BPND_AD['Global'].append(BPND_AD[sub][roi])
					BPnd_glb.append(BPND_AD[sub][roi])

				BPND_AD[roi].append(BPND_AD[sub][roi])

			BPND_AD['Global'].append(np.mean(BPnd_glb))


		elif positivity[sub] == 1:

			BPND_HC['subjects'].append(sub)
			BPND_HC[sub] = results_BPND[sub]

			BPnd_glb = []
			for roi in BPND_HC[sub]['ROIs']:

				if roi in cortex_global_rois:
					# BPND_HC['Global'].append(BPND_HC[sub][roi])
					BPnd_glb.append(BPND_HC[sub][roi])
	
				BPND_HC[roi].append(BPND_HC[sub][roi])
			
			BPND_HC['Global'].append(np.mean(BPnd_glb))

	return BPND_AD, BPND_HC

def compute_suvr_global(results):

	# rois_cortical_global = ['10009', '10010', '10011', '10013', '10012', '10014']
	rois_cortical_global = ['10009', '10010', '10011', '10012', '10014']

	# loop over all subjects
	for sub in results['subjects']:

		results[sub]['SUVr']['ROIs'].append('global')
		results[sub]['SUVr']['global'] = []

		# loop over all time frames of this sub
		for (i, t) in enumerate(results[sub]['SUVr']['time']):

			tmp_suvr_glb = [ results[sub]['SUVr'][roi][i] for roi in rois_cortical_global ]

			results[sub]['SUVr']['global'].append( np.mean(tmp_suvr_glb) )

#-----------------------
# statistics

def t_test_calc(quants_1, quants_2):
	t_stat, p_stat = st.ttest_ind(quants_1, quants_2, nan_policy='omit')
	return t_stat, p_stat

#-----------------------
# plot

colors = ['black', 'gray', 'red', 'orange', 'yellow', 'green', 'lime', 'cyan', 'blue', 'purple']
colors_nb = ['darkorange', 'sienna', 'steelblue', 'forestgreen', 'lightsteelblue', 'yellowgreen']


line_styles = ['-', '-.', ':']

markers_1 = ['o', 's', '^', 'v', '<', '>']
markers_2 = ['^', 'v', '<', '>', 'o', 's']

nc = len(colors)
ncb = len(colors_nb)

units = {
	'TACs': 'TAC [kBq/cc]',
	'SUVr': 'SUVr'
}

def plot_subject(sub_id, subject):

	for quant in ['TACs', 'SUVr']:

		handles = []
		labels = []
		title_save = 'plt_%s_%s.pdf' %(sub_id, quant)

		for (i, roi) in enumerate(subject[quant]['ROIs']):

			c_tmp = colors[int(i%nc)]
			ls_tmp = line_styles[ int((i-i%nc)/nc) ]

			labels.append(roi_definition(roi))

			hd, = plt.plot(subject[quant]['time'], subject[quant][roi], c=c_tmp, linestyle=ls_tmp)
			handles.append(hd)

		#if quant == 'SUVr':
		#	plt.ylim(0.6, 3.0)
		#	ncol = 3
		#else:
		#	ncol = 2
		ncol = 3

		plt.legend(handles=handles, labels=labels, ncol=ncol)
		plt.grid(True)
		plt.xlabel('time [min.]')
		plt.ylabel(units[quant])
		plt.savefig(title_save)
		plt.show()

def plot_means(average, quant):

	handles = []
	labels = []
	title_save = 'plt_mean_%s.pdf' %(quant)

	for (i, roi) in enumerate(average['ROIs']):

		c_tmp = colors[int(i%nc)]
		ls_tmp = line_styles[ int((i-i%nc)/nc) ]

		labels.append(roi_definition(roi))

		hd, = plt.plot(average[roi]['time'], average[roi]['vals'], c=c_tmp, linestyle=ls_tmp)
		handles.append(hd)

	if quant == 'SUVr':
		plt.ylim(0.5, 1.6)
		ncol=3
	else:
			ncol = 2

	plt.legend(handles=handles, labels=labels, ncol=ncol)
	plt.grid(True)
	plt.xlabel('time [min.]')
	plt.ylabel(units[quant])
	plt.savefig(title_save)
	plt.show()

def plot_DVR_vs_SUVr_temporal_no_color_diff(res_linregress, pairs, t, title=True):

	xs = np.linspace(1.0, 2.4, 20)
	ys = res_linregress['slope'] * xs + res_linregress['intercept']
	plt.plot(xs, ys, c='k', linestyle='-')

	plt.scatter(pairs['DVRs'], pairs['SUVr'], c='skyblue', marker='o')

	plt.grid(True)
	plt.xlabel('DVR')
	plt.ylabel('SUVR')
	plt.xlim(1.0, 2.4)
	if title:
		plt.title(res_linregress['equation'])
		plt.savefig('plt_DVR_SUVr_t_%s_titled.pdf' %(str(t)))
	else:
		plt.savefig('plt_DVR_SUVr_t_%s.pdf' %(str(t)))
	plt.show()

def plot_DVR_vs_SUVr_temporal(res_linregress, pairs_AD, pairs_HC, t, title=True):

	handles = []
	labels = []

	xs = np.linspace(1.0, 2.4, 20)
	ys = res_linregress['slope'] * xs + res_linregress['intercept']
	plt.plot(xs, ys, c='k', linestyle='-')

	i = 0
	for roi in pairs_AD['ROIs']:

		if (roi=='10002') or (roi=='10003') or (roi=='10004') or (roi=='10005') or (roi=='10006') or (roi=='10007'):

			hd = plt.scatter(pairs_AD[roi]['DVRs'], pairs_AD[roi]['SUVr'], c=colors_nb[i], marker='v')
			plt.scatter(pairs_HC[roi]['DVRs'], pairs_HC[roi]['SUVr'], c=colors_nb[i], marker='^')

			labels.append(roi_definition(roi))
			handles.append(hd)

			i += 1

	plt.legend(handles=handles, labels=labels)
	plt.grid(True)
	plt.xlabel('DVR')
	plt.ylabel('SUVR')
	plt.xlim(1.0, 2.4)
	if title:
		title = '%s , r = %f' %(res_linregress['equation'], res_linregress['rval'])
		# plt.title(title)
		plt.title(time_dict[str(t)])
		plt.savefig('plt_DVR_SUVr_roi_t_%s_titled.pdf' %(str(t)))
	else:
		plt.savefig('plt_DVR_SUVr_roi_t_%s.pdf' %(str(t)))
	plt.show()

def plot_DVR_VS_SUVr_together(res_linregress, pairs_AD, pairs_HC, title=True):

	# def set_ax_style(ax, labels):

	# 	# ax.set_xticks(np.arange(1, len(labels) + 1))
	# 	# ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
	# 	ax.tick_params(labelsize=8)

	# 	# ax.set_xlim(0.25, len(labels) + 0.75)
	# 	ax.set_xlabel('Group', fontproperties='Arial', fontsize=10)

	# 	ax.set_ylim( -0.1, 1.3)
	# 	ax.set_ylabel('BPnd', fontproperties='Arial', fontsize=10)

	s = 7

	# fontsize_1 = 6
	fontsize_1 = 8
	fontsize_2 = fontsize_1 + 2
	fontsize_3 = fontsize_2 + 2

	alpbet = ['A', 'B', 'C', 'D', 'E', 'F']
	# rois = ['10009', '10010', '10011', '10012', '10013', '10014']
	# rois = ['10009', '10010', '10011', '10012', '10014']
	rois = ['10011', '10009', '10014', '10010', '10012']

	# colors = ['lightsteelblue', 'darkorange', 'steelblue', 'forestgreen', 'yellowgreen']
	colors = ['steelblue', 'lightsteelblue', 'yellowgreen', 'darkorange', 'forestgreen']

	# width_cm = 14.0
	# width_cm = 20.0		# JNM
	width_cm = 19.0		# NeuroImage double column
	width_inch = width_cm / 2.54

	# for graphical abstract
	# width_cm_2 = 12.0

	# for manuscript
	width_cm_2 = 14.0
	# width_cm_2 = 24.0
	width_inch_2 = width_cm_2 / 2.54

	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(15,10))
	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(15,12))
	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(width_inch,width_inch_2), width_ratios=[0.3, 0.3, 0.3])
	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(width_inch,width_inch_2))

	fig, axss = plt.subplots(layout='constrained', nrows=2, ncols=3, figsize=(width_inch,width_inch_2), dpi=300)
	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(width_inch,30/2.54), dpi=300)

	# fig, axss = plt.subplots(layout='constrained', nrows=2, ncols=3, figsize=(width_inch,width_inch_2))
	# fig, axss = plt.subplots(nrows=3, ncols=2, figsize=(10, 15))
	# fig, axss = plt.subplots(nrows=2, ncols=3)

	for (i, t) in enumerate(res_linregress['times']):

		j = int(i%3)
		k = int( (i - int(i%3))/3 )
		# j = int(i%2)
		# k = int( (i - int(i%2))/2 )

		handles_1 = []
		labels_1 = []
		handles_2 = []
		labels_2 = []

		xs = np.linspace(1.0, 2.5, 20)
		ys = res_linregress[t]['slope'] * xs + res_linregress[t]['intercept']
		axss[k][j].plot(xs, ys, c='k', linestyle='-')

		l = 0
		# for roi in pairs_AD[res_linregress['times'][0]]['ROIs']:
		for roi in rois:

			# if (roi=='10002') or (roi=='10003') or (roi=='10004') or (roi=='10005') or (roi=='10006') or (roi=='10007'):
			if (roi=='10009') or (roi=='10010') or (roi=='10011') or (roi=='10013') or (roi=='10012') or (roi=='10014'):

				# hd1 = axss[k][j].scatter(pairs_AD[t][roi]['DVRs'], pairs_AD[t][roi]['SUVr'], c=colors_nb[l], marker='v')
				# hd2 = axss[k][j].scatter(pairs_HC[t][roi]['DVRs'], pairs_HC[t][roi]['SUVr'], c=colors_nb[l], marker='^')
				hd1 = axss[k][j].scatter(pairs_AD[t][roi]['DVRs'], pairs_AD[t][roi]['SUVr'], c=colors[l], marker='s', s=s)
				# hd2 = axss[k][j].scatter(pairs_HC[t][roi]['DVRs'], pairs_HC[t][roi]['SUVr'], c=colors_nb[l], marker='2')
				# hd2 = axss[k][j].scatter(pairs_HC[t][roi]['DVRs'], pairs_HC[t][roi]['SUVr'], c=colors_nb[l], marker='o', markerfacecolor='none')
				hd2 = axss[k][j].scatter(pairs_HC[t][roi]['DVRs'], pairs_HC[t][roi]['SUVr'], c='none', marker='o', edgecolor=colors[l], s=s)

				labels_1.append(roi_definition(roi)+', AD')
				handles_1.append(hd1)

				labels_2.append(roi_definition(roi)+', HC')
				handles_2.append(hd2)

				l += 1

		# str_tmp = '%s , r = %4.2f' %(res_linregress[t]['equation'], res_linregress[t]['rval'])
		str_tmp = '%s\nr = %4.2f' %(res_linregress[t]['equation'], res_linregress[t]['rval'])

		# if i < 2:
			# axss[k][j].text(1.2, 2.01, str_tmp, fontsize=10)
		# else:
			# axss[k][j].text(1.2, 0.81, str_tmp, fontsize=10)
		# axss[k][j].text(1.05, 2.25, alpbet[i], fontproperties='Helvetica', fontsize=12)
		if i < 2:
			# axss[k][j].text(1.2, 2.41, str_tmp, fontsize=fontsize_1)
			axss[k][j].text(1.15, 2.41, str_tmp, fontsize=fontsize_1)
		else:
			# axss[k][j].text(1.2, 1.01, str_tmp, fontsize=10)
			# axss[k][j].text(1.2, 0.8, str_tmp, fontsize=fontsize_1)
			axss[k][j].text(1.30, 0.8, str_tmp, fontsize=fontsize_1)

		# axss[k][j].text(1.05, 2.65, alpbet[i], fontproperties='Helvetica', fontsize=12)
		# axss[k][j].text(1.05, 3.05, alpbet[i], fontproperties='Helvetica', fontsize=fontsize_3)

		# axss[k][j].legend(handles=handles, labels=labels)
		# axss[k][j].grid(True)
		axss[k][j].set_xlabel('DVR', fontproperties='Arial', fontsize=fontsize_2)
		axss[k][j].set_ylabel('SUVR', fontproperties='Arial', fontsize=fontsize_2)

		axss[k][j].set_xlim( 1.0, 2.3)
		# axss[k][j].set_ylim( 0.6, 2.4)
		axss[k][j].set_ylim( 0.6, 2.8)

		axss[k][j].tick_params(labelsize=fontsize_1)

		axss[k][j].set_yticks(np.linspace(0.5, 3.0, 6))

		axss[k][j].set_title(time_dict[t], fontproperties='Arial', fontsize=fontsize_3)

	labels = []
	handles = []
	# for i in range(6):
	for i in range(len(rois)):
		# labels.append(labels_1[i])
		# labels.append(labels_2[i])
		# handles.append(handles_1[i])
		# handles.append(handles_2[i])
		labels.append(labels_2[i])
		labels.append(labels_1[i])
		handles.append(handles_2[i])
		handles.append(handles_1[i])

	# plt.legend(
	# 	handles=handles, 
	# 	labels=labels, 
	# 	ncols=3, 
	# 	# bbox_to_anchor=(1.02, 3.75)
	# 	# bbox_to_anchor=(1.02, 3.77)
	# 	bbox_to_anchor=(1.02, 2.5),
	# 	# bbox_to_anchor=(1.02, 4.07)
	# 	fontsize=fontsize_2
	# 	)

	# plt.savefig('plot_DVR_vs_SUVr.pdf')
	plt.savefig('plot_DVR_vs_SUVr.png')
	plt.show()

def plot_box_BPND_AD_HC_single(AD_quant, HC_quant):

	colors = ['lightblue', 'lightgreen']

	fig, ax = plt.subplots(nrows=1, ncols=1)

	bxplt = ax.boxplot(
		[AD_quant, HC_quant], 
		labels = ['AD', 'HC']
		)

	for patch, color in zip(bxplt['boxes'], colors):
		print(dir(patch))
		# patch.set_facecolor(color)
		patch.set_color(color)
		# patch.set_markerfacecolor(color)
		patch.set_markerfacecoloralt(color)

	ax.yaxis.grid(True)
	ax.set_xlabel('group')
	ax.set_ylabel('BPND')

	plt.show()

def plot_box_BPND_AD_HC_regional(AD_quants, HC_quants, quants_id):

	# width_cm_1 = 10.0
	width_cm_1 = 14.0
	width_inch_1 = width_cm_1 / 2.54

	# width_cm_2 = 10.0
	width_cm_2 = 12.0
	width_inch_2 = width_cm_2 / 2.54

	# colors = ['lightblue', 'lightgreen']
	colors = ['steelblue', 'darkorange']

	# s = 0.5
	s = 1.0

	fig, axss = plt.subplots(
		nrows=2, 
		ncols=3, 
		figsize=(width_inch_1, width_inch_2), 
		layout='constrained', 
		dpi=300
		)

	xs = [0, 1]

	for (i, qid) in enumerate(quants_id):

		j = int(i%3)
		k = int((i - int(i%3))/3)

		bxplt = axss[k][j].boxplot(
			# [AD_quants[i], HC_quants[i]], 
			[HC_quants[i], AD_quants[i]], 
			# labels = ['AD', 'HC'],
			labels = xs,
			widths=0.4,
			# showbars=False,
			showfliers=False
			)

		for patch in bxplt['whiskers']:
			patch.set_color('none')

		# for p in range(2):
		# 	bxplt['caps'][i].set_color(col)

		# axss[k][j].scatter(['AD' for a in AD_quants[i]], AD_quants[i], c='none', marker='o', edgecolor='k', s=5, linewidths=1)
		# axss[k][j].scatter(['HC' for a in HC_quants[i]], HC_quants[i], c='none', marker='o', edgecolor='k', s=5, linewidths=1)
		axss[k][j].scatter([2 for a in AD_quants[i]], AD_quants[i], c=colors[1], marker='o', s=s, zorder=30)
		axss[k][j].scatter([1 for a in HC_quants[i]], HC_quants[i], c=colors[0], marker='o', s=s, zorder=30)


		axss[k][j].set_title(roi_definition(qid), fontsize=12, font='arial')
		axss[k][j].set_ylabel('BPnd', fontsize=10, font='arial')
		axss[k][j].tick_params(labelsize=8)
		# axss[k][j].tick_params(labelsize=6)
		# axss[k][j].set_ylim(-0.1, 1.3)
		axss[k][j].set_ylim(-0.1, 1.4)

		axss[k][j].set_xticks([1, 2], labels=['HC','AD'])
		axss[k][j].set_yticks([0.0, 0.5, 1.0, 1.5])

		for patch in bxplt['medians']:
			patch.set_color('k')

	plt.savefig('plot_box_BPND_global_regional.png')
	plt.show()

def plot_box_BPND_AD_HC_global_regional(AD_quants, HC_quants, quants_id):

	# width_cm_1 = 14.0
	width_cm_1 = 11
	width_inch_1 = width_cm_1 / 2.54

	width_cm_2 = 8.0
	width_inch_2 = width_cm_2 / 2.54

	# colors = ['lightblue', 'lightgreen']
	colors = ['darkorange', 'sienna', 'steelblue', 'forestgreen', 'lightsteelblue', 'yellowgreen']

	xs = [1.0, 2.0]

	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(14,8))
	fig, axss = plt.subplots(nrows=1, ncols=1, figsize=(width_inch_1,width_inch_2), layout='constrained', dpi=300)

	handles = []
	labels = []

	for (i, qid) in enumerate(quants_id):

		if i == 0:
			bxplt = axss.boxplot(
				# [AD_quants[i], HC_quants[i]],
				[HC_quants[i], AD_quants[i]],
				widths=0.5,
				labels=xs,
				zorder=1
				)

			for patch in bxplt['medians']:
				patch.set_color('k')

		else:

			shift = -0.125 + (i-1) * 0.05

			# xs_AD = [ xs[0]+shift for j in range(len(AD_quants[i])) ]
			# xs_HC = [ xs[1]+shift for j in range(len(HC_quants[i])) ]
			xs_AD = [ xs[1]+shift for j in range(len(AD_quants[i])) ]
			xs_HC = [ xs[0]+shift for j in range(len(HC_quants[i])) ]

			hd = axss.scatter(xs_AD, AD_quants[i], color=colors[i-1], s=5)
			axss.scatter(xs_HC, HC_quants[i], color=colors[i-1], s=5)

			handles.append(hd)
			labels.append(roi_definition(qid))

	# axss.set_xticks(xs, labels=['AD', 'HC'])
	axss.set_xticks(xs, labels=['HC','AD'])
	axss.tick_params(labelsize=8)

	axss.set_ylabel('BPnd', fontsize=10, font='arial')

	plt.legend(
		handles=handles, 
		labels=labels, 
		ncol=2, 
		fontsize=8,
		columnspacing=0
		# reverse=True
		)

	plt.savefig('plot_box_BPND.png')
	plt.show()


	# for (i, qid) in enumerate(quants_id):

	# 	j = int(i%3)
	# 	k = int((i - int(i%3))/3)

	# 	bxplt = axss[k][j].boxplot(
	# 		[AD_quants[i], HC_quants[i]], 
	# 		labels = ['AD', 'HC']
	# 		)
	# 	axss[k][j].set_title(roi_definition(qid))

	# 	for patch, color in zip(bxplt['boxes'], colors):
	# 		patch.set_color(color)

	# for axs in axss:
	# 	for ax in axs:
	# 		# ax.yaxis.grid(True)
	# 		# ax.set_xlabel('group')
	# 		ax.set_ylabel('BPND')

	# # plt.savefig('plot_box_BPND.pdf')
	# plt.savefig('plot_box_BPND.png')
	# plt.show()

def plot_violin_BPND_AD_HC(AD_quants, HC_quants, quants_id):


	fontsize_1 = 8
	fontsize_2 = fontsize_1 + 2
	fontsize_3 = fontsize_2 + 2

	def set_ax_style(ax, labels):

		ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
		ax.tick_params(labelsize=fontsize_1)

		ax.set_xlim(0.25, len(labels) + 0.75)
		# ax.set_xlabel('Group', fontproperties='Arial', fontsize=fontsize_2)

		ax.set_ylim( -0.1, 1.3)
		ax.set_ylabel('BPnd', fontproperties='Arial', fontsize=fontsize_2)

	alpbet = ['A', 'B', 'C', 'D', 'E', 'F']
	colors_body = ['lightblue', 'lightgreen']
	colors_bars = ['blue', 'green']

	# width_cm_1 = 14.0
	width_cm_1 = 20.0
	width_inch_1 = width_cm_1 / 2.54

	width_cm_2 = 14.0
	# width_cm_2 = 24.0
	width_inch_2 = width_cm_2 / 2.54

	# fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(14,10))
	fig, axss = plt.subplots(layout='constrained', nrows=2, ncols=3, figsize=(width_inch_1,width_inch_2), dpi=300)

	for (i, qid) in enumerate(quants_id):

		j = int(i%3)
		k = int((i - int(i%3))/3)

		parts = axss[k][j].violinplot(
			[AD_quants[i], HC_quants[i]],
			showmeans=True,
			# showmedians=True,
			showextrema=True,
			# quantiles=[[0.2, 0.8], [0.2, 0.8]]
			)

		# print(parts)

		for pc in parts['bodies']:
			pc.set_facecolor('lightsteelblue')
			pc.set_alpha(1)

		parts['cmaxes'].set_color('steelblue')
		parts['cmins'].set_color('steelblue')
		parts['cmeans'].set_color('steelblue')
		parts['cbars'].set_color('steelblue')

		axss[k][j].set_title(roi_definition(qid), fontproperties='Arial', fontsize=fontsize_3)
		# axss[k][j].text(2.5, 1.2, alpbet[i], fontproperties='Helvetica', fontsize=12)
		set_ax_style(ax=axss[k][j], labels=['AD', 'HC'])

		tval, pval = t_test_calc(AD_quants[i], HC_quants[i])
		print('t-value = %f , p-value = %f' %(tval, pval))

	plt.savefig('plot_violin_BPND.png')
	plt.show()

def plot_scatter_BPND_AD_HC(AD_quants, HC_quants, quants_id):

	def set_ax_style(ax, labels):

		# ax.set_xticks(np.arange(1, len(labels) + 1))
		# ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
		ax.tick_params(labelsize=8)

		# ax.set_xlim(0.25, len(labels) + 0.75)
		ax.set_xlabel('Group', fontproperties='Arial', fontsize=10)

		ax.set_ylim( -0.1, 1.3)
		ax.set_ylabel('BPnd', fontproperties='Arial', fontsize=10)

	alpbet = ['A', 'B', 'C', 'D', 'E', 'F']
	colors_body = ['lightblue', 'lightgreen']
	colors_bars = ['blue', 'green']

	fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(14,10))

	for (i, qid) in enumerate(quants_id):

		num_AD = len(AD_quants[i])
		num_HC = len(HC_quants[i])

		ADs = ['AD' for s in range(num_AD)]
		HCs = ['HC' for s in range(num_HC)]
		# ADs = [ 0 for s in range(num_AD)]
		# HCs = [ 1 for s in range(num_HC)]

		j = int(i%3)
		k = int((i - int(i%3))/3)

		axss[k][j].scatter(ADs, AD_quants[i], c='none', marker='o', edgecolor=colors_bars[0])
		axss[k][j].scatter(HCs, HC_quants[i], c='none', marker='o', edgecolor=colors_bars[1])
		
		# parts = axss[k][j].violinplot(
		# 	[AD_quants[i], HC_quants[i]],
		# 	showmeans=True,
		# 	# showmedians=True,
		# 	showextrema=True,
		# 	# quantiles=[[0.2, 0.8], [0.2, 0.8]]
		# 	)

		axss[k][j].set_title(roi_definition(qid), fontproperties='Arial', fontsize=12)
		# axss[k][j].text(2.5, 1.2, alpbet[i], fontproperties='Helvetica', fontsize=12)
		set_ax_style(ax=axss[k][j], labels=['AD', 'HC'])

		tval, pval = t_test_calc(AD_quants[i], HC_quants[i])
		print('t-value = %f , p-value = %f' %(tval, pval))

	plt.savefig('plot_scatter_BPND.png')
	plt.show()

def plot_violin_BPND_AD_HC_w_pval(AD_quants, HC_quants, quants_id):

	def set_ax_style(ax, labels):

		ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
		ax.tick_params(labelsize=8)

		ax.set_xlim(0.25, len(labels) + 0.75)
		ax.set_xlabel('Group', fontproperties='Arial', fontsize=10)

		ax.set_ylim( -0.1, 1.5)
		ax.set_ylabel('BPND', fontproperties='Arial', fontsize=10)

	alpbet = ['A', 'B', 'C', 'D', 'E', 'F']
	colors_body = ['lightblue', 'lightgreen']
	colors_bars = ['blue', 'green']

	fig, axss = plt.subplots(nrows=2, ncols=3, figsize=(14,10))

	for (i, qid) in enumerate(quants_id):

		tval, pval = t_test_calc(AD_quants[i], HC_quants[i])
		print('t-value = %f , p-value = %f' %(tval, pval))

		j = int(i%3)
		k = int((i - int(i%3))/3)

		parts = axss[k][j].violinplot(
			[AD_quants[i], HC_quants[i]],
			showmeans=True,
			# showmedians=True,
			showextrema=True,
			# quantiles=[[0.2, 0.8], [0.2, 0.8]]
			)

		axss[k][j].set_title(roi_definition(qid), fontproperties='Arial', fontsize=12)
		axss[k][j].text(2.5, 1.35, alpbet[i], fontproperties='Helvetica', fontsize=12)
		set_ax_style(ax=axss[k][j], labels=['AD', 'HC'])

		if pval < 1E-6:
			text_tmp = 'p < 0.000001'
		elif (pval < 1E-5) and (pval >= 1E-6):
			text_tmp = 'p < 0.00001'
		elif (pval < 1E-4) and (pval >= 1E-5):
			text_tmp = 'p < 0.0001'
		axss[k][j].text(1.15, 1.4, text_tmp, fontproperties='Helvetica', fontsize=12)

		axss[k][j].vlines(1.0 , 1.29 , 1.35, color='k')
		axss[k][j].vlines(2.0 , 1.29 , 1.35, color='k')
		axss[k][j].hlines(1.35, 0.995, 2.005, color='k')


	# plt.savefig('plot_box_BPND.pdf')
	plt.savefig('plot_violin_BPND_w_pval.png')
	plt.show()

def plot_boxplot_BPND_AD_HC_w_roi_violin(AD_quants, HC_quants, quants_id):

	#----------------------- 

	xlabels = ['AD', 'HC']

	#----------------------- 

	fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 5))

	bxplt = ax.boxplot(
		[AD_quants[0], HC_quants[0]],
		widths=0.9
		)

	for patch in bxplt['boxes']:
		patch.set_color('k')

	#----------------------- 
	# plot violin for ROIs

	handles = []
	labels = []

	off_sets = [0.57, 1.57]

	for i in range(len(quants_id)-1):

		j = i + 1

		parts = ax.violinplot(
			[AD_quants[j], HC_quants[j]],
			positions=[off_sets[0]+0.125*j, off_sets[1]+0.125*j],
			widths=0.1
			)

		for pbody in parts['bodies']:
			pbody.set_color(colors_nb[i])
			pbody.set_edgecolor('none')

		parts['cmaxes'].set_color('none')
		parts['cmins'].set_color('none')
		parts['cbars'].set_color('none')

		handles.append(parts['bodies'][0])
		labels.append(roi_definition(quants_id[j]))


	#----------------------- 
	# manage axes

	ax.tick_params(labelsize=8)

	ax.set_xticks(np.arange(1, len(xlabels)+1), labels=xlabels)
	ax.set_xlim(0.25, len(xlabels)+0.75)
	ax.set_xlabel('Group', fontproperties='Arial', fontsize=10)

	ax.set_ylim( -0.1, 1.3)
	ax.set_ylabel('BPND', fontproperties='Arial', fontsize=10)

	#----------------------- 

	plt.legend(handles=handles, labels=labels)

	plt.savefig('plot_violin_BPND_w_roi.png')
	plt.show()


#=======================
#
#=======================

