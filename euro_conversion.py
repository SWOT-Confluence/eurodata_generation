


# reach_nxs[21303000131]


# Link fix - this is the one
from netCDF4 import Dataset
import glob
import shutil
import pandas as pd
import os
import numpy as np
import math

def load_file(topology_file_path):
    matched_df = pd.read_csv(topology_file_path)
    return matched_df



from netCDF4 import Dataset
import glob
import shutil
import pandas as pd
import os
from collections import Counter
import subprocess as sp
import time
from datetime import datetime, timedelta


def get_reach_nodes(sword_path, reach_id):

    all_nodes = []

    # files = glob.glob(os.path.join(data_dir, '*'))
    # print(f'Searching across {len(files)} continents for nodes...')


    rootgrp = Dataset(sword_path, "r", format="NETCDF4")

    node_ids_indexes = np.where(rootgrp.groups['nodes'].variables['reach_id'][:].data.astype('U') == str(reach_id))

    if len(node_ids_indexes[0])!=0:
        for y in node_ids_indexes[0]:
            node_id = str(rootgrp.groups['nodes'].variables['node_id'][y].data.astype('U'))
            all_nodes.append(node_id)



            # all_nodes.extend(node_ids[0].tolist())

        rootgrp.close()

    print(f'Found {len(set(all_nodes))} nodes...')
    return list(set(all_nodes))

sword_path = '/nas/cee-water/cjgleason/travis/data/mnt/input/sword/eu_sword_v11.nc'
sword = Dataset(sword_path)
def find_reach_nx(sword):
    reach_nxs = Counter(sword['nodes']['reach_id'][:])
    return reach_nxs
reach_nxs = find_reach_nx(sword)
matched = glob.glob('/nas/cee-water/cjgleason/euro_data/euro_lisflood/Euro_sword_topology/*')

cnt = 0
for topology_file_path in matched:
    print(cnt, 'of', len(matched))
    cnt +=1
    matched_df = load_file(topology_file_path)
    # print(topology_file_path)
    # loop through links
    reach_ids = matched_df.reach_id.unique()
    
    
    #restart function
    # for sword_reach in linked_df.reach_id.unique():
    restart_paths = [os.path.exists(f'/nas/cee-water/cjgleason/travis/out_testing/{sword_reach}.nc') for sword_reach in reach_ids]
    # print(restart_paths)

    if False in restart_paths:
        

        for l in reach_ids:
            
        # these reach netcdfs seem to be messed up not sure of the reason
            if l not in [21307900361, 23303600495, 21105300021, 22738400051, 22590300211] :
                if not math.isnan(l):
                    # print('link is', l)

                    linked_df = matched_df[matched_df['reach_id'] == l]



                    #find euro reachs associated with that reach_id

                    euro_links = linked_df.link.unique()






            #         for euro_link in euro_links:


            #             # load each variable width 
            #             basin = os.path.basename(topology_file_path).split('_')[0]

            #             nc_file_path = f'/nas/cee-water/cjgleason/travis/data/euro_swot_nc/{basin}_{euro_link}_SWOT.nc'

            #             #be sure it exists in test set
            #             if os.path.exists(nc_file_path):
            #                 netcd = Dataset(nc_file_path, mode = 'r+', clobber=True)









                    agg_var_dict = {
                        'reach' : {
                            'd_x_area': [], 
                            'slope2': [], 
                            'wse' :[],
                            'width':[]
                        },

                        # 'node': {
                        #     'd_x_area': [], 
                        #     'slope2': [], 
                        #     'wse' :[],
                        #     'width':[]
                        # },
                    }

                    for euro_link in euro_links:





                        # load each variable width 
                        basin = os.path.basename(topology_file_path).split('_')[0]

                        nc_file_path = f'/nas/cee-water/cjgleason/travis/data/euro_swot_nc/{basin}_{euro_link}_SWOT.nc'

                        var_width_csv = f'/nas/cee-water/cjgleason/yuta/laplace/data/euro_data/euro_variable_width/out/wfull_lisflood_dfull_sos/{basin}_{euro_link}.csv'

                        var_width_df = pd.read_csv(var_width_csv)

                        #be sure it exists in test set
                        if os.path.exists(nc_file_path):
                            netcd = Dataset(nc_file_path, mode = 'r+', clobber=True)

                            for level in agg_var_dict.keys():

                                for agg_var in agg_var_dict[level].keys():

                                    # load variable widths for each
                                    if agg_var == 'width':
                                        if len(var_width_df.variable_width.unique()) > 1:

                                            var_values = list(var_width_df['variable_width'])
                                        else:
                                            var_values = []

                                    else: 

                                        var_values = netcd[level][agg_var][:]

                                    # check if this is the first time we are adding values
                                    if len(agg_var_dict[level][agg_var]) == 0:
                                        agg_var_dict[level][agg_var] = [[a_value] for a_value in var_values]

                                    # make a sub list elementwize
                                    else:
                                        # print(level, agg_var, var_values[0])
                                        # try:
                                        #     print(len(var_values), 'by', len(var_values[0]))
                                        # except:
                                        #     continue



                                        # therea are different nubmers of nodes that we are trying to aggregate
                                        # email out to see how we should handle it
                                        # print(nc_file_path)
                                        # print(agg_var, len(agg_var_dict[level][agg_var]), 'by', len(agg_var_dict[level][agg_var][0]))

                                        for index_v, a_value in enumerate(var_values):
                                            agg_var_dict[level][agg_var][index_v].append(a_value)



                    # average their data together and store 

                    for level in agg_var_dict.keys():
                        for agg_var in agg_var_dict[level].keys():

                            for index_v, a_value in enumerate(agg_var_dict[level][agg_var]):

                                agg_var_dict[level][agg_var][index_v] = np.ma.asarray(agg_var_dict[level][agg_var][index_v]).mean(axis=0)



                    # use last nc_file_path to do the sword reach nc write out

                    shutil.copyfile(nc_file_path, f'/nas/cee-water/cjgleason/travis/out_testing/{l}.nc')

                    # open new nc file
                    new_fp = f'/nas/cee-water/cjgleason/travis/out_testing/{l}.nc'
                    netcd = Dataset(new_fp, mode = 'r+', clobber=True)

                    print(f'/nas/cee-water/cjgleason/travis/out_testing/{l}.nc')

                    def overwrite_variable(agg_var_dict, agg_var, netcd, target_nx):

                        # resize_dims
                        starting_dim = len(netcd['node'][agg_var])
                        while starting_dim < target_nx:

                            print(starting_dim, target_nx, agg_var)
                            netcd.close()
                            extend_cmd = ['ncks', '--msa_usr_rdr', '-d', f'nx,0,{starting_dim-1}', '-d', f'nx,0,{starting_dim-1}', new_fp, new_fp, '-O']
                            sp.run(extend_cmd)
                            netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                            starting_dim = len(netcd['node'][agg_var])


                        netcd.close()
                        shrink_cmd = ['ncks', '-d', f'nx,0,{target_nx-1}', new_fp, new_fp, '-O']
                        sp.run(shrink_cmd)
                        netcd = Dataset(new_fp, mode = 'r+', clobber=True)


                        # overwrite widths of netcdf
                        # we have to check if width has values because if all the var width dfs associated with this reach had nan values
                        # we want to keep constant the ones it has now
                        if agg_var == 'width':
                            if len(agg_var_dict['reach']['width']) != 0:
                                print(agg_var, len(agg_var_dict['reach'][agg_var]))
                                # we move our width selection up by 1 becaues the starting withs are nan values
                                # in this step we are also only selecting nt number of widths
                                netcd['reach'][agg_var][:] = agg_var_dict['reach'][agg_var][1:(len(netcd['reach'][agg_var][:])+1)]

                                node_len = len(netcd['node'][agg_var][:])
                                new_nodes = [agg_var_dict['reach'][agg_var][1:(len(netcd['reach'][agg_var][:])+1)] for i in range(node_len)]
                                netcd['node'][agg_var][:] = new_nodes

                        else:
                            netcd['reach'][agg_var][:] = agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:]))]
                            new_nodes = [agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:]))] for i in range(len(netcd['node'][agg_var][:]))]
                            netcd['node'][agg_var][:] = new_nodes




                    for agg_var in agg_var_dict['reach'].keys():
                        # print(agg_var, 'ov')
                        target_nx = reach_nxs[int(l)]
                        overwrite_variable(agg_var_dict, agg_var, netcd, target_nx)
                        
                        
                    # add in missing variables
                    # netcd.close()
                    new_fp = f'/nas/cee-water/cjgleason/travis/out_testing/{l}.nc'
                    netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                    nx = netcd.dimensions['nx']
                    nt = netcd.dimensions['nt']
                    
                    
                    the_times = [] 

                    for i in range(len(nt)):
                        # input datetime
                        dt = datetime.now()
                        # epoch time
                        epoch_time = datetime(1970, 1, 1)

                        dt_1 = dt + timedelta(days = -i)

                        # subtract Datetime from epoch datetime
                        delta = (dt_1 - epoch_time)
                        the_times.append(delta.total_seconds())

                    node_ids = get_reach_nodes(sword_path, l)
                    int_node_ids = [int(i) for i in node_ids]
                    netcd['node'].createVariable('time','i8', (nx, nt))
                    netcd['node'].createVariable('node_id', 'i8', (nx))
                    netcd['reach'].createVariable('time', 'i8', (nt))
                    
                    netcd['node']['node_id'][:] = int_node_ids
                    netcd['node']['time'][:] = [the_times for i in range(len(nx))]
                    netcd['reach']['time'][:] = the_times
                    # netcd.close()
                    
                    # What is time??
                    
                    
                    