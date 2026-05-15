import os
import time
import netCDF4 as nc
from math import log10, sqrt, isnan
from matplotlib import pyplot as plt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
This script will start a chain of simulations with the objective
of iteratively finding suitable sliding coefficients. It is
using the same method found in the technical report Greve et
al., 2020 (DOI: 10.5281/zenodo.3971232). The only change is the
addition of a modifier, allowing for a smoother evolution of the
sliding coefficients (Dangleterre, 2023, DOI: 10.5281/zenodo.8409491).

You need to make sure the python library netCDF4 is installed
through
  python3.11 -m pip install netcdf4
on the machine. (A newer version of Python will also do.)

Change to the main SICOPOLIS directory, and make a local copy
of the template:
  cp ./tools/diverse/iterative_sliding.py ./my_iterative_sliding_xxx.py
(where 'xxx' can be any descriptor of the user's choice).

An initial 'startup' header is needed, which will be run first
('0th iteration'). Choose a suitable name for this file, and give it
as a string to the variable 'name_of_run'. Subsequent iterations
follow the naming convention 'name_of_run_iteration_modifier'.

The entry for C_SLIDE_DIMLESS must be in the line specified by the
variable 'line_number_c_slide' (to be set below).

The program will write several computed parameters in a folder
'./tmp/iterative_sliding_{name_of_run}', namely:
  rmsd_lin
  slope_lin
  rmsd_log
  slope_log
  intercept_log
in their own txt files. Another Python script is provided for plotting
out of the box.

Execution of the script:
  python3.11 my_iterative_sliding_xxx.py
(A newer version of Python will also do.)

SICOPOLIS output must be in the standard directory './sico_out';
otherwise, the script won't work!

Created by Tom Dangleterre
Last update: 2026-05-15 by Ralf Greve
'''

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------- VARIABLES TO FILL --------

name_of_run = 'grl16_bm6_spinup11_cal_100ka_iter'
# Name of the 'startup' header ('0th' iteration)

dx = '16'  # Resolution

line_number_c_slide = 1179
# Line number of the entry for C_SLIDE_DIMLESS in the 'startup' header.
# '1140' for Antarctica / '1179' for Greenland should work.

line_num = line_number_c_slide-1
# Python indexing starts with 0, therefore the '-1' is needed

kmax = 7  # Maximum number of iterations (typically 5-15)

modifier = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# Values of the modifier (relaxation factor) for the iterations;
# typically 0.5-1.0 (the smaller, the less aggressive).
# Length of array must be kmax+1 (k=0...kmax).
# The value modifier[0] for the zeroth iteration is a dummy (has no effect).

if len(modifier) != kmax+1:
    print(f'Length of \'modifier\' array does not match value of \'kmax\'!')
    exit()

tslice = '0001'
# Time-slice number for final state of iterations k=0...kmax
# (which is used for comparison with observed surface velocities)

anfdatname = 'grl16_bm6_spinup11_cal_100ka'  # Name of the initial-conditions simulation
# For example ant32_bm3_jare_aq1_spinup03_holocene_1
# NOT TO USE:
#   ant32_bm3_jare_aq1_spinup03_holocene_10002.nc
#   (no time-slice number, no extension)

targetname = 'grl16_bm6_spinup11_smooth_100a'  # Name of the target simulation for nudging

# Regions file, to change according to the ice sheet
path = './sico_in/grl/'
filename = f'grl{dx}_zwally2012_basins_with-negis_extrapolated.nc'

regions = nc.Dataset(path+filename)
reg = regions['n_basin'][:]
n_reg = max(max(sub) for sub in reg)

# Change the following depending on the ice sheet studied, dx is the resolution entered earlier
# Adjust the path to the MEaSUREs surface velocities maps 
path = './sico_in/grl/'
filename = f'SurfVel_Greenland_MEaSUREs_GridEPSG3413_{dx}km.nc'

#-------- FOLDER FOR COMPUTED PARAMETERS --------

try :
    os.mkdir(f'./tmp/iterative_sliding_{name_of_run}')
except:
    print(f'The folder \'iterative_sliding_{name_of_run}\' already exists. Please change its name or move it elsewhere and try again.')
    exit()

obs = nc.Dataset(path+filename)
vs_obs = obs['vs'][:].tolist()

#-------- USEFUL FUNCTIONS --------

def matrixmap(T):
    """Maps log10 to a square matrix"""
    ni = len(T)
    nj = len(T[0])
    result = []
    for i in range(ni):
        result.append([])
        for j in range(nj):
            if (T[i][j] == 0):
                result[i].append(float('nan'))
            else:
                result[i].append(log10(T[i][j]))
    return result 

def get_intercepts(k):
    '''Calculates the intercepts for each regions, takes in the 
    iteration number, and returns the slopes, rmsd and intercept
    for computing the sliding coefficient in each region.'''

    # dir = ''
    dir = '.'
    path = f'{dir}/sico_out/{name_of_run}_{k:02d}_{modifier[k]}/'
    filename = f'{name_of_run}_{k:02d}_{modifier[k]}{tslice}.nc'

    sim = nc.Dataset(path+filename)
    vs_sim = sim['vh_s'][:].tolist()
    mask = sim['mask'][:]
    n_cts = sim['n_cts'][:]

    if (len(vs_obs) != len(vs_sim)):
        print(len(vs_obs), len(vs_sim))
        raise Exception("Sim and Obs are not of the same size")

    imax = len(vs_obs)
    jmax = len(vs_obs[0])

    log_vs_obs = matrixmap(vs_obs)
    log_vs_sim = matrixmap(vs_sim) 

    n_data        = []
    rmsd_lin      = []
    slope_lin     = []
    slope_log     = []
    rmsd_log      = []
    intercept_log = []

    print(f'Calculating the new sliding coefficients for iteration {k+1}...')
    for n_cnt in range(1, n_reg+2):
        log_vs_obs_aux = [log_vs_obs[k][:] for k in range(imax)]
        vs_obs_aux = [vs_obs[k][:] for k in range(imax)]
        log_vs_sim_aux = [log_vs_sim[k][:] for k in range(imax)]
        vs_sim_aux = [vs_sim[k][:] for k in range(imax)]

        # Cleaning up in order to remove the ocean 
        if (n_cnt <= n_reg):
            for i in range(imax):
                for j in range(jmax):
                    if (reg[i][j] != n_cnt or mask[i][j] != 0):
                        log_vs_obs_aux[i][j] = float('nan')
                        vs_obs_aux[i][j] = float('nan')
                        log_vs_sim_aux[i][j] = float('nan')
                        vs_sim_aux[i][j] = float('nan')

        else:
            for i in range(imax):
                for j in range(jmax):
                    if (mask[i][j] != 0):
                        log_vs_obs_aux[i][j] = float('nan')
                        vs_obs_aux[i][j] = float('nan')
                        log_vs_sim_aux[i][j] = float('nan')
                        vs_sim_aux[i][j] = float('nan')

        # Calculating rmsd
        n              = 0
        rmsd_lin1      = 0
        slope_lin1     = 0
        slope_lin2     = 0
        rmsd_log1      = 0
        intercept_log1 = 0
        axis_log_min = log10(10)
        axis_log_max = log10(3000)
        for i in range(imax):
            for j in range(jmax):
                if (log_vs_obs_aux[i][j] >= axis_log_min and log_vs_obs_aux[i][j] <= axis_log_max) & \
                    (log_vs_sim_aux[i][j] >= axis_log_min and log_vs_sim_aux[i][j] <= axis_log_max):
                    n += 1
                    rmsd_lin1 = rmsd_lin1 + (vs_sim[i][j] - vs_obs[i][j])**2
                    if (isnan(rmsd_lin1)):
                        raise Exception('asdasd')
                    slope_lin1 = slope_lin1 + vs_obs[i][j] * vs_sim[i][j]
                    slope_lin2 = slope_lin2 + vs_obs[i][j]**2
                    rmsd_log1 = rmsd_log1 + (log_vs_sim[i][j] - log_vs_obs[i][j])**2
                    intercept_log1 = intercept_log1 + (log_vs_sim[i][j] - log_vs_obs[i][j])

        n_data.append(n)
        if n > 0 :
            rmsd_lin.append(sqrt(rmsd_lin1 / n))
            slope_lin.append(slope_lin1 / slope_lin2)
            rmsd_log.append(sqrt(rmsd_log1 / n))
            res = intercept_log1 / n
            intercept_log.append(res)
            slope_log.append(10**res)

    # Writing the data
    print(f'Writing data for iteration {k+1}')
    dir2 = './tmp'

    with open(f'{dir2}/iterative_sliding_{name_of_run}/rmsd_lin.txt', 'a') as f:
        g = (f'{rmsd_lin[k]},' for k in range(len(rmsd_lin)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'{dir2}/iterative_sliding_{name_of_run}/slope_lin.txt', 'a') as f:
        g = (f'{slope_lin[k]},' for k in range(len(slope_lin)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'{dir2}/iterative_sliding_{name_of_run}/rmsd_log.txt', 'a') as f:
        g = (f'{rmsd_log[k]},' for k in range(len(rmsd_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'{dir2}/iterative_sliding_{name_of_run}/intercept_log.txt', 'a') as f:
        g = (f'{intercept_log[k]},' for k in range(len(intercept_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'{dir2}/iterative_sliding_{name_of_run}/slope_log.txt', 'a') as f:
        g = (f'{slope_log[k]},' for k in range(len(slope_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')
                    
    return rmsd_lin, slope_lin, rmsd_log, intercept_log, slope_log

#-------- ZEROTH ITERATION --------

k = 0

header = open(f'./headers/sico_specs_{name_of_run}.h', 'r')
content = header.readlines()
header.close()

new_file = open(f'./headers/sico_specs_{name_of_run}_{k:02d}_{modifier[k]}.h','w')
new_file.write(''.join(content))
new_file.close()

sico_sh_opts = f'-m {name_of_run}_{k:02d}_{modifier[k]}'
sico_sh_opts = f'{sico_sh_opts} -a ${{PWD}}/sico_out/{anfdatname}'
sico_sh_opts = f'{sico_sh_opts} -t ${{PWD}}/sico_out/{targetname}'

time.sleep(5)
print(f'Running iteration {k}...')
os.system(f'(./sico.sh {sico_sh_opts}) >tmp/out_{name_of_run}_{k:02d}_{modifier[k]}.dat 2>&1')
time.sleep(15)
print(f'Iteration {k} finished...')

#-------- ITERATIONS --------

while k <= kmax:
    # Compute fit between simulated and observed surface velocities
    a, slope_lin, a, a, slope_log = get_intercepts(k)

    # Opens the previous iterations' header in order to copy it into the new header 
    header = open(f'./headers/sico_specs_{name_of_run}_{k:02d}_{modifier[k]}.h', 'r')
    content = header.readlines()
    header.close()
    a = content[line_num]
    a = a.replace('#define C_SLIDE_DIMLESS [', '')
    a = a.replace(' ', '')
    a = a.replace('\\', '')
    a = a.replace('d0', '')
    a = a.replace(']', '')
    a = a.replace('\n', '')
    a = a.split(',')
    print(a)

    slide = [float(a[k]) for k in range(len(a))]
    # For debugging
    # print('slope_log length', len(slope_log))
    # print('sliding coeff length', len(slide))
    # print([round(slope_log[j], 4) for j in range(len(slope_log))])
    for j in range(len(slope_log)):
        slope_log[j] = slope_log[j] * modifier[k] + (1 - modifier[k])
    # For debugging
    # print('for testing...')
    # print([round(slope_log[j], 4) for j in range(len(slope_log))])
    # print('slope_log length', len(slope_log))
    # print('sliding coeff length', len(slide))
    C = [max(round(slide[k]/slope_log[k],4),0.01) for k in range(n_reg)]

    values = f'#define C_SLIDE_DIMLESS [ '
    for j in range(n_reg):
        values += f'{C[j]}d0, '
    values = values[:-2]
    values +=  ']\n'

    content[line_num] = values

    #-------- CREATION OF THE NEXT ITERATION

    k = k + 1

    if k <= kmax:
        new_file = open(f'./headers/sico_specs_{name_of_run}_{k:02d}_{modifier[k]}.h','w')
        new_file.write(''.join(content))
        new_file.close()

        sico_sh_opts = f'-m {name_of_run}_{k:02d}_{modifier[k]}'
        sico_sh_opts = f'{sico_sh_opts} -a ${{PWD}}/sico_out/{anfdatname}'
        sico_sh_opts = f'{sico_sh_opts} -t ${{PWD}}/sico_out/{targetname}'

        #-------- RUNNING THE NEXT ITERATION

        time.sleep(5)
        print(f'Running iteration {k}...')
        os.system(f'(./sico.sh {sico_sh_opts}) >tmp/out_{name_of_run}_{k:02d}_{modifier[k]}.dat 2>&1')
        time.sleep(15)
        print(f'Iteration {k} finished...')

#-------- END OF SCRIPT --------

print(f'kmax = {kmax}, thus iteration {k} will not be carried out. Job done.')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
