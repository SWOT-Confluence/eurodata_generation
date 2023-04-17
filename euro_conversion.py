import warnings
warnings.filterwarnings("ignore")


# reach_nxs[21303000131]


# Link fix - this is the one
from netCDF4 import Dataset
import glob
import shutil
import pandas as pd
import os
import numpy as np
import math
from collections import Counter
import subprocess as sp
import time
from datetime import datetime, timedelta
import json



global problem_ids
problem_ids = []

def find_largest_reach(upstream_ids, sword):

    all_width = []

    for i in upstream_ids:
        index = np.where(np.asarray(sword['reaches']['reach_id']) == i)
        width = sword['reaches']['max_width'][index]
        all_width.append(width.data[0])
    
    max_value = max(all_width)
    upstream_id = upstream_ids[all_width.index(max_value)]
    
    return upstream_id


def pull_upstream_data(largest_reach, outdir, agg_var):
    # check and see if this largest reach has been processed yet
    if os.path.exists(os.path.join(outdir, str(largest_reach)+'.nc')):
        print('----------------------------------------------------1')
        upstream = Dataset(os.path.join(outdir, str(largest_reach)+'.nc'))
        # print('reached processed, extracting values')
        # return upstream

        # modify agg_var_dict
        print('values found for upstream reach...')
        print(upstream['reach'][agg_var][:10])
        return upstream['reach'][agg_var]
        # raise
    else:
        return 'skipping'


def update_problem_json(reach_of_interest):
    global problem_ids

    problem_ids.append(reach_of_interest)
    problem_ids = [str(i) for i in list(set(problem_ids))]
    print('upstream reach not processed, skipping')

    if os.path.exists('retry.json'):
        with open('retry.json', 'r') as f:
            data = json.load(f)
            data['ids'] = problem_ids
            f.close()
    else:
        start_out = {'ids': problem_ids}
        with open('retry.json',"w") as f:
            json.dump(start_out,f)
            f.close()


def find_upstream_replacement(reach_of_interest, agg_var, sword, outdir):
    # global problem_ids
    # domain = int(str(reach_of_interest)[-1])-1
    index = np.where(np.asarray(sword['reaches']['reach_id']) == reach_of_interest)
    # print(sword['reaches']['rch_id_up'][:].T[index])
    # print(index)
    try:
        upstream_ids = [i for i in sword['reaches']['rch_id_up'][:].T[index][0] if i != 0]
    except:
        upstream_ids = []
    num_upstream_ids = len(upstream_ids)


    if num_upstream_ids > 1:
        print('many upstream')
        largest_reach = find_largest_reach(upstream_ids, sword=sword)
        upstream_data = pull_upstream_data(largest_reach=largest_reach, outdir= outdir, agg_var = agg_var)

    elif num_upstream_ids == 0:
        # do nothing
        print(reach_of_interest, ', is a river origin, no upstream rivers to replace ', agg_var)
        upstream_data = 'skipping'

    elif num_upstream_ids ==1:
        print('one_upstream')
        replacement_id = upstream_ids[0]
        upstream_data = pull_upstream_data(largest_reach = replacement_id, outdir = outdir, agg_var = agg_var)

    # check if upstream data is nan, if so add it to the skipping list
        
    return upstream_data


def load_file(topology_file_path):
    matched_df = pd.read_csv(topology_file_path)
    return matched_df

def get_reach_nodes(sword_path, reach_id):

    all_nodes = []

    # files = glob.glob(os.path.join(data_dir, '*'))
    # print(f'Searching across {len(files)} continents for nodes...')

    print('-------------------------------2')
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


def resize_nt(netcd, new_fp, target_nx=9362):

    # resize_dims
    # starting_dim = len(netcd['reach'][agg_var])
    # starting_dim = netcd['reach'].dimensions['nt'].size
    starting_dim = len(netcd['reach']['width'][:])
    # raise
    print('attempting to resize', netcd)
    if starting_dim != target_nx:
        print('nt resizing.....')
        while starting_dim < target_nx:

            # print(starting_dim, target_nx, agg_var)
            netcd.close()
            extend_cmd = ['ncks', '--msa_usr_rdr', '-d', f'nt,0,{starting_dim-1}', '-d', f'nt,0,{starting_dim-1}', new_fp, new_fp, '-O']
            sp.run(extend_cmd)
            print('-------------------------------3')
            netcd = Dataset(new_fp, mode = 'r+', clobber=True)
            starting_dim = starting_dim = len(netcd['reach']['width'][:])


        netcd.close()
        shrink_cmd = ['ncks', '-d', f'nt,0,{target_nx-1}', new_fp, new_fp, '-O']
        sp.run(shrink_cmd)
        print('-------------------------4')
        netcd = Dataset(new_fp, mode = 'r+', clobber=True)


def overwrite_variable(agg_var_dict, agg_var, target_nx, new_fp):


    # resize_dims
    print('------------------5')
    netcd = Dataset(new_fp, mode = 'r+', clobber=True)
    starting_dim = len(netcd['node'][agg_var])
    if starting_dim != target_nx:
        while starting_dim < target_nx:

            # print(starting_dim, target_nx, agg_var)
            netcd.close()
            extend_cmd = ['ncks', '--msa_usr_rdr', '-d', f'nx,0,{starting_dim-1}', '-d', f'nx,0,{starting_dim-1}', new_fp, new_fp, '-O']
            sp.run(extend_cmd)
            print('------------------------6')
            netcd = Dataset(new_fp, mode = 'r+', clobber=True)
            starting_dim = len(netcd['node'][agg_var][:])


        netcd.close()
        shrink_cmd = ['ncks', '-d', f'nx,0,{target_nx-1}', new_fp, new_fp, '-O']
        sp.run(shrink_cmd)
        print('---------------------------7')
        netcd = Dataset(new_fp, mode = 'r+', clobber=True)

    # print('-------------after shorten / lengthen----------------')
    # print('currently empty')
    # print(netcd['node']['wse'][:])


    # overwrite widths of netcdf
    # we have to check if width has values because if all the var width dfs associated with this reach had nan values
    # we want to keep constant the ones it has now

    if agg_var == 'width':
        # checking to see if there are varied withs
        if len(set(np.asarray(agg_var_dict['reach'][agg_var]).astype(str))) > 1:
            # print('1111111111111111111111111111111111')
            # print(agg_var, len(np.asarray(agg_var_dict['reach']['width']).astype(str)))
            # we move our width selection up by 1 becaues the starting withs are nan values
            # in this step we are also only selecting nt number of widths

            # change here
            netcd['reach'][agg_var][:] = agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:])+1)]

            node_len = len(netcd['node'][agg_var][:])
            new_nodes = [agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:])+1)] for i in range(node_len)]
            netcd['node'][agg_var][:] = new_nodes
        else:
            print(f'----------------------------------------------------keeping constant width, may be nan')
            # print(agg_var_dict['reach'][agg_var])


    else:
        # print('2222222222222222222')

        # if there is aggregated data for the reach variable
        # add it to the reach file and duplicate it for nodes
        # we cannot use the node data from the og euro nc files

        # if there is not then we have to keep the data from the previous nc files
        # this presents a problem 

        # print('--',agg_var_dict['reach'][agg_var])
        if len(set(np.asarray(agg_var_dict['reach'][agg_var]).astype(str))) > 1:
            # print('3333333333333333')
            # len(set(np.asarray(agg_var_dict['reach'][agg_var]).astype(str)))
            netcd['reach'][agg_var][:] = agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:]))]
            # print(netcd['reach'][agg_var][:])
            # print(agg_var_dict['reach'][agg_var])
            # print(agg_var_dict['reach'][agg_var][:(len(netcd['reach'][agg_var][:]))])
            # print([agg_var_dict['reach'][agg_var][:len(netcd['reach'][agg_var][:])]])
            # print(netcd['node'])
            # print(netcd['node'][agg_var][:])
            # raise
            test = len(netcd['reach'][agg_var][:])

            new_nodes = [agg_var_dict['reach'][agg_var][:test] for i in range(len(netcd['node'][agg_var][:]))]
            try:
                netcd['node'][agg_var][:] = new_nodes
            except:
                # netcd.close()
                try:
                    print('---------------7')
                    netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                except:
                    pass
                netcd['node'][agg_var][:] = new_nodes

        else:
            global failed_reaches
            failed_reaches.append(str(l)+'foo')
            # print(f'reach {l} failed ------------------------------- :(')
            print(set(failed_reaches))
            print(len(set(failed_reaches)))
            for x in range(7):
                print(str(x), sum(z.endswith(f'{x}foo') for z in set(failed_reaches)))
                # print(str(i), sum(str(i) in str(s) for s in failed_reaches)
    try:
        netcd.close()
    except:
        pass

#------------------------main-------------------------------

external_path = '/media/confluence/work/'
external_path = '/mnt/work/'

sword_path = '/nas/cee-water/cjgleason/travis/data/mnt/input/sword/eu_sword_v11.nc'
sword_path = os.path.join(external_path, 'data/euro_mnt/input/sword/eu_sword_v11.nc')

matched = glob.glob('/nas/cee-water/cjgleason/euro_data/euro_lisflood/Euro_sword_topology/*')
matched = glob.glob(os.path.join(external_path, 'data/euro_gen/Euro_sword_topology/*'))

original_eurodata_nc_dir = '/nas/cee-water/cjgleason/travis/data/euro_swot_nc/'
original_eurodata_nc_dir = os.path.join(external_path, 'data/euro_gen/euro_swot_nc/')


var_width_dir = '/nas/cee-water/cjgleason/yuta/laplace/data/euro_data/euro_variable_width/out/wfull_lisflood_dfull_sos/'
var_width_dir = os.path.join(external_path, 'data/euro_gen/wfull_lisflood_dfull_sos/')

out_path = '/nas/cee-water/cjgleason/travis/out_testing/'

out_path = os.path.join(external_path, 'data/euro_gen/euro_conversion_out/')


sword = Dataset(sword_path)
def find_reach_nx(sword):
    reach_nxs = Counter(sword['nodes']['reach_id'][:])
    return reach_nxs
reach_nxs = find_reach_nx(sword)

global failed_reaches
failed_reaches = []

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
    restart_paths = [os.path.exists(os.path.join(out_path, f'{sword_reach}.nc')) for sword_reach in reach_ids]
    # print(restart_paths)

    if False in restart_paths:
        

        for l in reach_ids:
            print(' ')
            print('Processing reach', l)
            
        # these reach netcdfs seem to be messed up not sure of the reason
            if l not in [21307900361, 23303600495, 21105300021, 22738400051, 22590300211] :
                if not math.isnan(l):
                    print('link is', l)

                    linked_df = matched_df[matched_df['reach_id'] == l]



                    #find euro reachs associated with that reach_id

                    euro_links = linked_df.link.unique()


                    # go through each of the euro nc file that are associated with a single reach_id
                    # average their variables
                    # if any of these are 
                    agg_var_dict = {
                        'reach' : {
                            'd_x_area': [], 
                            'slope2': [], 
                            'wse' :[],
                            'width':[]
                        },

                    }



                    for euro_link in list(set(euro_links)):


                        # load each variable width 
                        basin = os.path.basename(topology_file_path).split('_')[0]

                        nc_file_path = os.path.join(original_eurodata_nc_dir, f'{basin}_{euro_link}_SWOT.nc')

                        var_width_csv = os.path.join(var_width_dir, f'{basin}_{euro_link}.csv')

                        var_width_df = pd.read_csv(var_width_csv)
                        # print(var_width_df)
                        print('earlier')
                        new_fp = os.path.join(out_path, f'{l}.nc')



                        # if os.path.exists(nc_file_path):
                        #     if not os.path.exists(new_fp):

                        #         shutil.copyfile(nc_file_path, new_fp)
                        #         print('-------------------------8')
                        #         netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                        #         print('starting_dims', len(netcd['reach']['width'][:]))

                        #         resize_nt(netcd=netcd, new_fp=new_fp, target_nx=9362)

                        #         print('ending dims', len(netcd['reach']['width'][:]))

                        #         print(euro_link, 'much earlier')


                        #be sure it exists in test set
                        if os.path.exists(nc_file_path):
                            try:
                                netcd.close()
                            except:
                                pass
                            print(nc_file_path)

                            # open one of the og eurodata paths that are associated with a basin


                            # this can be moved to after averaging
                            # new_fp = os.path.join(out_path, f'{l}.nc')

                            # shutil.copyfile(nc_file_path, new_fp)
                            # print('-------------------------8')

                            # netcd = Dataset(new_fp, mode = 'r+', clobber=True)

                            netcd = Dataset(nc_file_path, mode = 'r+', clobber=True)
                            print('starting_dims', len(netcd['reach']['width'][:]))

                            resize_nt(netcd=netcd, new_fp=nc_file_path, target_nx=9362)

                            print('ending dims', len(netcd['reach']['width'][:]))


                            for level in agg_var_dict.keys():

                                for agg_var in agg_var_dict[level].keys():

                                    # load variable widths from widths dataframe
                                    if agg_var == 'width':
                                        if len(var_width_df.variable_width.unique()) > 1:
                                            # then widths are variable, update the width in the new nc
                                            var_values = list(var_width_df['variable_width'])

                                        elif len(var_width_df.variable_width.unique()) == 1:
                                            

                                            # if the one value is nans, grab the width from the og netcdf
                                            if str(list(var_width_df.variable_width.unique())[0]) == 'nan':

                                                # if the og width is greater than one unique value
                                                # then return those if not, do the if below
                                                if len(list(set([str(i) for i in netcd['reach']['width'][:]]))) > 1:
                                                       var_values = netcd['reach']['width'][:]

                                                
                                                elif list(set([str(i) for i in netcd['reach']['width'][:]]))[0] == '--':
                                                    # we allready checked to see if there is more than one value in the og widths
                                                    # if we don't, then we get here. 
                                                    # we have also allready checked to see if the variable widths is nan
                                                    # both the variable length df and og netcdf file were nans
                                                    print(f'both the variable length df and og netcdf file were nans for width')
                                                    var_values = []
                                                    
                                                
                                                
                                                else:
                                                    # variable widths df is nans, and og is constant
                                                    var_values = netcd['reach']['width']



                                            else:
                                                # variable widths is constant, and og file is nans
                                                var_values = list(var_width_df['variable_width'])
                                                # print('here is constant', var_values)



                                    # if the variable is not widths, grab the vaules from the og netcdf to start with , these will be averaged
                                    # which could be all nans
                                    else:
                                        # if list(set([str(i) for i in netcd['reach'][agg_var][:]]))[0] == '--':                                        
                                        #     var_values = netcd[level][agg_var][:]
                                        #     # print('is this nan?')
                                        #     print(list(set([str(i) for i in netcd['reach'][agg_var][:]]))[0])
                                        var_values = netcd[level][agg_var][:]


                                    # some sorces have variying nt, push them all to be 9363
                                    # var_values = var_values[:len(netcd['reach'][agg_var][:])]
                                    print('pre selection', len(var_values), agg_var)


                                    var_values = var_values[:9362]

                                    print('post selection', len(var_values), agg_var)


                                    # check if this is the first time we are adding values
                                    # if so initalize the list for the rest of the values to be included
                                    if len(agg_var_dict[level][agg_var]) == 0:
                                        agg_var_dict[level][agg_var] = [[a_value] for a_value in var_values]
                                        print('-------first through')


                                    # make a sub list elementwize
                                    else:

                                        # therea are different nubmers of nodes that we are trying to aggregate
                                        # email out to see how we should handle it
                                        # was decided to just copy reach level variables

                                        for index_v, a_value in enumerate(var_values):
                                            # print(a_value, agg_var)
                                            print(print('first pass', len(agg_var_dict[level][agg_var][:]), 'trying to update', len(var_values)), 'starting',len(netcd['reach'][agg_var][:]))

                                            agg_var_dict[level][agg_var][index_v].append(a_value)


# -----------------current dev
                    netcd.close()

                    # average their data together and store 
                    # at this point many of the variables could be masked arrays if all of the data from the combined og netcdfs was nan
                    print('averaging data...')
                    for level in agg_var_dict.keys():
                        for agg_var in agg_var_dict[level].keys():
                            print(agg_var)
                            # print('testing', agg_var_dict[level][agg_var][0])
                            for index_v, a_value in enumerate(agg_var_dict[level][agg_var]):
                                agg_var_dict[level][agg_var][index_v] = np.ma.asarray(agg_var_dict[level][agg_var][index_v]).mean(axis=0)


                    # now go through the averaged dictionary and check out if there are any empty variables
                    print('processing averaged data to check for empty values...')
                    for level in agg_var_dict.keys():
                        for agg_var in agg_var_dict[level].keys():
                            
                            # get all the unique values of the variable
                            all_string = list(set([str(i) for i in agg_var_dict[level][agg_var]]))

                            # if the set all string is less than or equal to one
                            # this means that the variable is either a valid constant (left alone), empty, or nan (replaced with upstream data)
                            if len(all_string) <= 1:

                                if len(all_string) == 0:
                                    # then the variable is empty
                                    # set to upstream reach value
                                    print(f'all empty values found in {agg_var}, attempting to find upstream replacement')
                                    upstream_data = find_upstream_replacement(reach_of_interest=l, agg_var=agg_var, sword=sword, outdir=out_path)
                                    
                                    
                                    # if there wan't an upstream value to set it too, it will be skipped
                                    # it will also be saved for later to reprocess
                                    if upstream_data != 'skipping':
                                        agg_var_dict[level][agg_var] = upstream_data

                                        # --------------------write function here that checks if a variable in the dict is gtg-------------------------

                                    else:
                                        update_problem_json(reach_of_interest=l)
                                    
                                    

                                # then the variable is constant, check to see if it is constantly nan
                                elif all_string[0] == '--':

                                    # set to upstream node value
                                    print(f'all nan values found in {agg_var}, attempting to find upstream replacement')
                                    upstream_data = find_upstream_replacement(reach_of_interest=l, agg_var=agg_var, sword=sword, outdir=out_path)
                                    if upstream_data != 'skipping':
                                        agg_var_dict[level][agg_var] = upstream_data

                                    else:
                                        update_problem_json(reach_of_interest=l)

                                else:
                                    print(agg_var, ' was valid and constant')

                            else:
                                print(agg_var, ' was valid and variable')
                                print(agg_var_dict[level][agg_var][:10])

                    # new_fp = os.path.join(out_path, f'{l}.nc')

                    # shutil.copyfile(nc_file_path, new_fp)

                    # open new nc file
                    # netcd = Dataset(new_fp, mode = 'r+', clobber=True)

                    # print(netcd)



                    new_fp = os.path.join(out_path, f'{l}.nc')

                    shutil.copyfile(nc_file_path, new_fp)
                    print('-------------------------8')
                    netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                    resize_nt(netcd=netcd,new_fp=new_fp,target_nx=9362)
                    netcd.close()

                    # netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                 

                    for agg_var in agg_var_dict['reach'].keys():
                        # print(agg_var, 'ov')
                        target_nx = reach_nxs[int(l)]
                        overwrite_variable(agg_var_dict=agg_var_dict, agg_var=agg_var, target_nx=target_nx, new_fp=new_fp)
                        
                    

                    # add in missing variables
                    # netcd.close()
                    # new_fp = os.path.join(out_path, f'{l}.nc')
                    print('----------------------8')
                    netcd = Dataset(new_fp, mode = 'r+', clobber=True)
                    nx = netcd.dimensions['nx']
                    nt = netcd.dimensions['nt']

                    # print('-------------final wse----------------')
                    # print('currently empty')
                    # print(netcd['node']['wse'][:])
                    
                    
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
                    print(' ')
                    int_node_ids = [int(i) for i in node_ids]
                    netcd['node'].createVariable('time','i8', (nx, nt))
                    netcd['node'].createVariable('node_id', 'i8', (nx))
                    netcd['reach'].createVariable('time', 'i8', (nt))
                    netcd['node']['node_id'][:] = int_node_ids
                    netcd['node']['time'][:] = [the_times for i in range(len(nx))]
                    netcd['reach']['time'][:] = the_times
                    # print('added time varaibles')
                    # netcd.close()
                    
                    # What is time??
                    
                    
                    