import os
import time
import netCDF4 as nc
from math import log10, sqrt, isnan
from matplotlib import pyplot as plt

'''This script will start a chain of simulations with the
objective of iteratively finding suitable sliding coefficients. I
is using the same method found in the technical report Greve et
al., 2020 (DOI: 10.5281/zenodo.3971232). The only change is the
addition of a modifier, allowing for a smoother evolution of the
sliding coefficients (See 10.5281/zenodo.8409492 for more
details).

You need to make sure the python library netCDF4 is installed
through pip install netcdf4 on the machine.

An initial 'startup' header is needed, with all the parameters
that will remain unchanged throughout this process already filled
and checked. 
Choose a suitable name for this file, and give it as a string to
the variable 'name_of_run' as it will be carried out
for all iterations following the naming convention
name_modifier_iteration.

There should also be a file named my_multi_sico_iter_k.sh in the
/run folder.

The initial conditions for this iterative scheme requires a
previous run with sliding coefficients of 1m/a/Pa everywhere. 

All of the prefiled variables are examples, and will not work as they
are.

The 'startup' file should have the following parameters to work
correctly here

#define N_SLIDE_REGIONS 18

#define SLIDE_REGIONS_FILE 'ant32_imbie2016_basins_extrapolated.nc'

#define C_SLIDE (/ 0.5434d0, 0.4839d0, 0.7838d0, 0.5294d0, 0.4299d0, 0.377d0, \
                    0.5604d0, 0.7653d0, 1.397d0, 0.4005d0, 0.4616d0, 0.3646d0, \
                    0.7906d0, 1.0639d0, 0.7929d0, 0.6035d0, 1.1895d0, 0.3525d0 /)

With the 3 lines used by C_SLIDE at lines 1016, 1017 and 1018. If not please change in this code
the line numbers accordingly. The values used here for C_SLIDE will be
changed in the code and are arbitrary (just used as a template), as long as the correct
synthaxing is used.

The program will output in a folder /results_{name_of_run}
different parameters, namely :
rmsd_lin
slope_lin
rmsd_log
res
intercept_log
slope_log
in their own txt files. Another Python script is provided for
plotting out of the box.

Created by Tom Dangleterre
Last update: 09-10-2023
 '''


######################## VARIABLES TO FILL ############################


name_of_run = 'ant32_bm3_jare_aq1_spinup03_fixtopo_1' #name of the 'startup' file described earlier
modifier = 0.5  #value of the modifier, the change in the sliding coefficients is too agressive without it
dx = 32    #Resolution
anfdatname = 'ant32_bm3_jare_aq1_spinup03_holocene_1' #Name of the initial conditions file (previous simulation)
#For example ant32_bm3_jare_aq1_spinup03_holocene_1, 
# NOT TO USE:
#ant32_bm3_jare_aq1_spinup03_holocene_10002.nc
#No extensions
targetname = 'ant32_bm3_jare_aq1_spinup03_holocene_init100a' #name of the target for nudging.
tslice1 = '0002' #Time slice number for the initial conditions
tslice2 = '0003' #Time slice number for the final ice sheet state, after evolution with the new sliding coefficients 
#tslice1 and tslice2 are related to the initial conditions runs
#at 1m/a/Pa. 

t0 = -4500
tfinal = -4000
#Start and end times of the iterations simulations

kmax = 5 # Max number of iterations, 10-15 is usually enough

#######################################################################


######## INITIALIZATION

try :
    os.mkdir(f'./results_{name_of_run}')
except:
    print(f'the results folder with the name results_{name_of_run} already exists. Please Change it\'s name or move it elsewhere and try again')
    exit()

# regions file, to change according to the ice sheet.
path = '../sico_in/ant/'
filename = f'ant{dx}_imbie2016_basins_extrapolated.nc'

regions = nc.Dataset(path+filename)
n_reg = regions['n_basin'][:]


path = './surfvel/'
#Change the following depending on the ice sheet studied, dx is the resolution entered earlier
filename = f'SurfVel_Antarctica_MEaSUREs_GridEPSG3031_{dx}km.nc'
ob = nc.Dataset(path+filename)
vs_ob = ob['vs'][:].tolist()

################## USEFUL FUNCTIONS

def matrixmap(T):
    """Maps log10 to a square matrix"""
    n = len(T)
    result = []
    for i in range(n):
        result.append([])
        for j in range(n):
            if (T[i][j] == 0):
                result[i].append(float('nan'))
            else:
                result[i].append(log10(T[i][j]))
    return result 


def get_intercepts(k, first_iter = False):
    '''Calculates the intercepts for each regions, takes in the 
    iteration number, and returns the slopes, rmsd and intercept
    for computing the sliding coefficient in each region.'''


    if first_iter:
        path = f'../sico_out/{anfdatname}/'
        filename = f'{anfdatname}{tslice2}.nc'
    else:
        path = f'../sico_out/{name_of_run}_{modifier}_{k}/'
        filename = f'{name_of_run}_{modifier}_{k}0001.nc'


    sim = nc.Dataset(path+filename)
    vs_sim = sim['vh_s'][:].tolist()
    mask = sim['mask'][:]
    n_cts = sim['n_cts'][:]

    if (len(vs_ob) != len(vs_sim)):
        raise Exception("Sim and Observed are not of the same size")

    nmax = len(vs_ob)


    log_vs_ob = matrixmap(vs_ob)
    log_vs_sim = matrixmap(vs_sim) 


    n_data = []
    rmsd_lin     = []
    slope_lin     = []
    slope_log     = []
    rmsd_log      = []
    intercept_log = []

    print(f'Calculating the new Sliding Coefficients for iteration {k+1}...')
    for n_cnt in range(1, max(max(sub) for sub in n_reg)+2):
        log_vs_ob_aux = [log_vs_ob[k][:] for k in range(nmax)]
        vs_ob_aux = [vs_ob[k][:] for k in range(nmax)]
        log_vs_sim_aux = [log_vs_sim[k][:] for k in range(nmax)]
        vs_sim_aux = [vs_sim[k][:] for k in range(nmax)]


        #Cleaning up in order to remove the ocean 
        if (n_cnt <= max([max(sub) for sub in n_reg])):
            for i in range(nmax):
                for j in range(nmax):
                    if (n_reg[i][j] != n_cnt or mask[i][j] != 0):
                        log_vs_ob_aux[i][j] = float('nan')
                        vs_ob_aux[i][j] = float('nan')
                        log_vs_sim_aux[i][j] = float('nan')
                        vs_sim_aux[i][j] = float('nan')

        if (n_cnt > max([max(sub) for sub in n_reg])):
            for i in range(nmax):
                for j in range(nmax):
                    if (mask[i][j] != 0):
                        log_vs_ob_aux[i][j] = float('nan')
                        vs_ob_aux[i][j] = float('nan')
                        log_vs_sim_aux[i][j] = float('nan')
                        vs_sim_aux[i][j] = float('nan')

        #Calculating rmsd
        n              = 0
        rmsd_lin1      = 0
        slope_lin1     = 0
        slope_lin2     = 0
        rmsd_log1      = 0
        intercept_log1 = 0
        axis_log_min = log10(10)
        axis_log_max = log10(3000)
        for i in range(nmax):
            for j in range(nmax):
                if (log_vs_ob_aux[i][j] >= axis_log_min and log_vs_ob_aux[i][j] <= axis_log_max) & \
                    (log_vs_sim_aux[i][j] >= axis_log_min and log_vs_sim_aux[i][j] <= axis_log_max):
                    n += 1
                    rmsd_lin1 = rmsd_lin1 + (vs_sim[i][j] - vs_ob[i][j])**2
                    if (isnan(rmsd_lin1)):
                        raise Exception('asdasd')
                    slope_lin1 = slope_lin1 + vs_ob[i][j] * vs_sim[i][j]
                    slope_lin2 = slope_lin2 + vs_ob[i][j]**2
                    rmsd_log1 = rmsd_log1 + (log_vs_sim[i][j] - log_vs_ob[i][j])**2
                    intercept_log1 = intercept_log1 + (log_vs_sim[i][j] - log_vs_ob[i][j])

        n_data.append(n)
        if n > 0 :
            rmsd_lin.append(sqrt(rmsd_lin1 / n))
            slope_lin.append(slope_lin1 / slope_lin2)
            rmsd_log.append(sqrt(rmsd_log1 / n))
            res = intercept_log1 / n
            intercept_log.append(res)
            slope_log.append(10**res)

    #Writing the data
    print(f'Writing data for iteration {k+1}')
    with open(f'results_{name_of_run}/{modifier}_rmsd_lin.txt', 'a') as f:
        g = (f'{rmsd_lin[k]},' for k in range(len(rmsd_lin)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'results_{name_of_run}/{modifier}_slope_lin.txt', 'a') as f:
        g = (f'{slope_lin[k]},' for k in range(len(slope_lin)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'results_{name_of_run}/{modifier}_rmsd_log.txt', 'a') as f:
        g = (f'{rmsd_log[k]},' for k in range(len(rmsd_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'results_{name_of_run}/{modifier}_intercept_log.txt', 'a') as f:
        g = (f'{intercept_log[k]},' for k in range(len(intercept_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')

    with open(f'results_{name_of_run}/{modifier}_slope_log.txt', 'a') as f:
        g = (f'{slope_log[k]},' for k in range(len(slope_log)))
        for x in g:
            f.write(str(x))
        f.write('\n')
                
    
    return rmsd_lin, slope_lin, rmsd_log, intercept_log, slope_log


################## CREATION OF THE FIRST ITERATION

header = open(f'./headers/{name_of_run}.h', 'r')
content = header.readlines()
header.close()

rmsd_lin, slope_lin, rmsd_log, intercept_log, slope_log = get_intercepts(0, True)

a,b, c = content[1015], content[1016], content[1017]
a = a.replace('#define C_SLIDE (/ ', '')
a = a.replace(' ', '')
a = a.replace('\\', '')
a = a.replace('d0', '')
a = a.replace('\n', '')
a = a.split(',')
a.pop()
b = b.replace(' ', '')
b = b.replace('\\', '')
b = b.replace('d0', '')
b = b.replace('\n', '')
b = b.split(',')
b.pop()
c = c.replace(' ', '')
c = c.replace('\\', '')
c = c.replace('d0', '')
c = c.replace('\n', '')
c = c.replace('/)', '')
c = c.split(',')
combi = a + b + c
slide = [float(combi[k]) for k in range(len(combi))]

for j in range(len(slope_log)):
    slope_log[j] = slope_log[j] * modifier + ( 1 - modifier)
C = [max(round(slide[k]/slope_log[k],4),0.1) for k in range(18)]

content[468] = f'#define ANFDATNAME \'{anfdatname}{tslice1}.nc\'\n'
content[101] = f'#define TIME_INIT0 {t0}.0d0 \n'
content[104] = f'#define TIME_END0 {tfinal}.0d0 \n'
content[1285] = f'#define N_OUTPUT 1'
content[1290] = f'#define TIME_OUT0 {tfinal}'


content[1015] = f'#define C_SLIDE (/ {C[0]}d0, {C[1]}d0, {C[2]}d0, {C[3]}d0, {C[4]}d0, {C[5]}d0, \\\n'
content[1016] = f'                    {C[6]}d0, {C[7]}d0, {C[8]}d0, {C[9]}d0, {C[10]}d0, {C[11]}d0, \\\n'
content[1017] = f'                    {C[12]}d0, {C[13]}d0, {C[14]}d0, {C[15]}d0, {C[16]}d0, {C[17]}d0 /)\n'


new_file = open(f'headers/{name_of_run}_{modifier}_1.h','w')
new_file.write(''.join(content))
new_file.close()


################### RUNNING THE SIMULATION

run = open(f'my_multi_sico_iter_k.sh', 'r')
run_ctn = run.readlines()
run.close()

run_ctn[113] = f'   (./sico.sh ${{MULTI_OPTIONS_1}} -m {name_of_run}_{modifier}_1 \\\n'
run_ctn[114] = f'              -a ${{MULTI_OUTDIR}}/{anfdatname} \\\n'
run_ctn[115] = f'              -t ${{MULTI_OUTDIR}}/{targetname}) \\\n'
run_ctn[116] = f'              >out_multi_{name_of_run}_{modifier}_1.dat 2>&1\n'

run = open(f'my_multi_sico_iter_{name_of_run}_{modifier}_1.sh', 'w')
run.write(''.join(run_ctn))
run.close()

time.sleep(5)
os.chmod(f'my_multi_sico_iter{name_of_run}_{modifier}_1.sh', 0o755 )
time.sleep(5)
os.system(f'(./my_multi_sico_iter_{name_of_run}_{modifier}_1.sh) >out_multi_iter_{name_of_run}_{modifier}_1.dat 2>&1 &')


################ ITERATIONS

k = 1
while k < kmax:
    '''The program will wait 10 minutes before checking if the
    previous simulation is done, can be changed through time_to_wait'''

    time_to_wait = 600 #in seconds


    iterator = True
    waited_time = 0
    print(f'running simulation of iteration {k}...')
    while iterator:
        try:
            a, slope_lin, a, a, slope_log = get_intercepts(k)
            iterator = False
        except:
            time.sleep(time_to_wait)
            waited_time = waited_time + 10
            print(f'time Waited : {waited_time} minutes')

    #Opens the previous iteration' header in order to copy it into the new header 
    header = open(f'./headers/sico_specs_{name_of_run}_{modifier}_{k}.h', 'r')
    content = header.readlines()
    header.close()
    a,b, c = content[1015], content[1016], content[1017]
    a = a.replace('#define C_SLIDE (/ ', '')
    a = a.replace(' ', '')
    a = a.replace('\\', '')
    a = a.replace('d0', '')
    a = a.replace('\n', '')
    a = a.split(',')
    a.pop()
    b = b.replace(' ', '')
    b = b.replace('\\', '')
    b = b.replace('d0', '')
    b = b.replace('\n', '')
    b = b.split(',')
    b.pop()
    c = c.replace(' ', '')
    c = c.replace('\\', '')
    c = c.replace('d0', '')
    c = c.replace('\n', '')
    c = c.replace('/)', '')
    c = c.split(',')
    combi = a + b + c
    slide = [float(combi[k]) for k in range(len(combi))]
    # for debugging
    # print('slope_log length', len(slope_log))
    # print('sliding coeff length', len(slide))
    # print([round(slope_log[j], 4) for j in range(len(slope_log))])
    for j in range(len(slope_log)):
        slope_log[j] = slope_log[j] * modifier + ( 1 - modifier)
    # For debugging
    # print('for testing...')
    # print([round(slope_log[j], 4) for j in range(len(slope_log))])
    # print('slope_log length', len(slope_log))
    # print('sliding coeff length', len(slide))
    C = [max(round(slide[k]/slope_log[k],4),0.1) for k in range(18)]

    content[1015] = f'#define C_SLIDE (/ {C[0]}d0, {C[1]}d0, {C[2]}d0, {C[3]}d0, {C[4]}d0, {C[5]}d0, \\\n'
    content[1016] = f'                    {C[6]}d0, {C[7]}d0, {C[8]}d0, {C[9]}d0, {C[10]}d0, {C[11]}d0, \\\n'
    content[1017] = f'                    {C[12]}d0, {C[13]}d0, {C[14]}d0, {C[15]}d0, {C[16]}d0, {C[17]}d0 /)\n'

    k = k + 1

    #################### CREATION OF THE NEXT ITERATION
    new_file = open(f'headers/sico_specs_{name_of_run}_{modifier}_{k}.h','w')
    new_file.write(''.join(content))
    new_file.close()

    run = open('my_multi_sico_iter_k.sh', 'r')
    run_ctn = run.readlines()
    run.close()

    #Create new run.sh
    run_ctn[113] = f'   (./sico.sh ${{MULTI_OPTIONS_1}} -m {name_of_run}_{modifier}_{k} \\\n'
    run_ctn[114] = f'              -a ${{MULTI_OUTDIR}}/{anfdatname} \\\n'
    run_ctn[115] = f'              -t ${{MULTI_OUTDIR}}/{targetname}) \\\n'
    run_ctn[116] = f'              >out_multi_{name_of_run}_{modifier}_{k}.dat 2>&1\n'

    run = open(f'my_multi_sico_iter_{name_of_run}_{modifier}_{k}.sh', 'w')
    run.write(''.join(run_ctn))
    run.close()

    ################# RUNING NEW ITERATION
    time.sleep(5)
    os.chmod(f'my_multi_sico_iter_{name_of_run}_{modifier}_{k}.sh', 0o755 )
    time.sleep(5)
    os.system(f'(./my_multi_sico_iter_{name_of_run}_{modifier}_{k}.sh) >out_multi_iter_{name_of_run}_{modifier}_{k}.dat 2>&1 &')

print('Job done. It worked on the first time. Probably.')