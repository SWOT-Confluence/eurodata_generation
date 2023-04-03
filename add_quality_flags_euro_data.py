# imports

from netCDF4 import Dataset
import glob
import numpy as np
import random


flag_dict = {
    'node':
        {'node_q': {
        'good' : 0, 'bad': 1
    },
        'dark_frac': {
        'good' : 0.0, 'bad': 1.0
    },
        'ice_clim_f': {
        'good' : 0, 'medium':1, 'bad': 2
    },
        'ice_dyn_f': {
        'good' : 0, 'medium':1, 'bad': 2
    },
        'partial_f': {
        'good' : 0, 'bad': 1
    },
        'n_good_pix': {
        'good' : 100000, 'bad': 0
    },
        'xovr_cal_q': {
        'good' : 0, 'bad': 1
    }},


'reach' : {
        'reach_q': {
        'good' : 0, 'bad': 1
    },
        'dark_frac': {
        'good' : 0.0, 'bad': 1.0
    },
        'ice_clim_f': {
        'good' : 0, 'medium':1, 'bad': 2
    },
        'ice_dyn_f': {
        'good' : 0, 'medium':1, 'bad': 2
    },
        'partial_f': {
        'good' : 0, 'bad': 1
    },
        'n_good_nod': {
        'good' : 100, 'bad': 0
    },
        'xovr_cal_q': {
        'good' : 0, 'bad': 1
    },
        'obs_frac_n': {
        'good' : 1.0, 'bad': 0.0
    },

}}




# importing Counter module

# next do value assignment
from collections import Counter

def update_node_flags(flag_dict, nc):
    for level in ['node']:
        for flag in flag_dict[level].keys():
            try:
                nc[level].createVariable(flag, 'f8', ('nx', 'nt'))
            except:
                continue
            

            # three levels
            if len(flag_dict[level][flag]) > 2:
                ice_clim_array = np.asarray(nc[level][flag][:])
                ice_clim_array.fill(flag_dict[level][flag]['good'])

                for index, nx in enumerate(ice_clim_array):


                    index_value = random.sample(list(enumerate(nx)), round(len(nx)*0.035))

                    all_indexes_to_replace = list(list(zip(*index_value))[0])
                    # print(all_indexes_to_replace)
                    random.shuffle(all_indexes_to_replace)
                    medium_indexes_to_replace = all_indexes_to_replace[round(len(all_indexes_to_replace)/2):]
                    bad_indexes_to_replace = all_indexes_to_replace[:round(len(all_indexes_to_replace)/2)]
                    
                    # print(medium_indexes_to_replace, bad_indexes_to_replace)

                    np.put(nx, medium_indexes_to_replace, flag_dict[level][flag]['medium'])
                    np.put(nx, bad_indexes_to_replace, flag_dict[level][flag]['bad'])

                    ice_clim_array[index] = list(nx)


                # frequency = Counter(nx)
                # print(level, flag, frequency)

                nc[level][flag][:] = list(ice_clim_array)

            # two levels
            else:
                ice_clim_array = np.asarray(nc[level][flag][:])
                ice_clim_array.fill(flag_dict[level][flag]['good'])

                for index, nx in enumerate(ice_clim_array):


                    index_value = random.sample(list(enumerate(nx)), round(len(nx)*0.035))

                    all_indexes_to_replace = list(zip(*index_value))[0]

                    np.put(nx, all_indexes_to_replace, flag_dict[level][flag]['bad'])

                    ice_clim_array[index] = list(nx)
                    # print(nx)
                # frequency = Counter(nx)
                # print(level, flag, frequency)
            
                nc[level][flag][:] = list(ice_clim_array)
    # nc.close()
                        
def update_reach_flags(flag_dict, nc):
    for level in ['reach']:
        for flag in flag_dict[level].keys():
            try:
                nc[level].createVariable(flag, 'f8', ('nt'))
            except:
                continue
            

            # three levels
            if len(flag_dict[level][flag]) > 2:
                ice_clim_list = np.asarray(nc[level][flag][:])
                ice_clim_list.fill(flag_dict[level][flag]['good'])


                index_value = random.sample(list(enumerate(ice_clim_list)), round(len(ice_clim_list)*0.035))

                all_indexes_to_replace = list(list(zip(*index_value))[0])
                # print(all_indexes_to_replace)
                random.shuffle(all_indexes_to_replace)
                medium_indexes_to_replace = all_indexes_to_replace[round(len(all_indexes_to_replace)/2):]
                bad_indexes_to_replace = all_indexes_to_replace[:round(len(all_indexes_to_replace)/2)]

                # print(medium_indexes_to_replace, bad_indexes_to_replace)

                np.put(ice_clim_list, medium_indexes_to_replace, flag_dict[level][flag]['medium'])
                np.put(ice_clim_list, bad_indexes_to_replace, flag_dict[level][flag]['bad'])


                # frequency = Counter(nx)
                # print(level, flag, frequency)

                nc[level][flag][:] = list(ice_clim_list)

            # two levels
            else:
                ice_clim_list = np.asarray(nc[level][flag][:])
                ice_clim_list.fill(flag_dict[level][flag]['good'])



                index_value = random.sample(list(enumerate(ice_clim_list)), round(len(ice_clim_list)*0.035))

                all_indexes_to_replace = list(zip(*index_value))[0]

                np.put(ice_clim_list, all_indexes_to_replace, flag_dict[level][flag]['bad'])

                    #print(nx)
                # frequency = Counter(nx)
                # print(level, flag, frequency)
            
                nc[level][flag][:] = list(ice_clim_list)
    # nc.close()
    
    
all_files = glob.glob('/nas/cee-water/cjgleason/travis/out_testing/*.nc')
cnt = len(all_files)
for i in all_files:
    if str(cnt).endswith('000'):
        print(cnt, ' files left ..')
    cnt -=1
    try:
        og_file = Dataset(i, 'a')
        update_node_flags(flag_dict=flag_dict, nc=og_file)
        update_reach_flags(flag_dict=flag_dict, nc=og_file)
        og_file.close()
    except:
        print(i, 'failed')
        try:
            og_file.close()
        except:
            pass
        pass
print('complete')
