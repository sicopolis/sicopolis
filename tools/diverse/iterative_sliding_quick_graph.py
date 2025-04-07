from statistics import mean
from matplotlib import pyplot as plt
from matplotlib import ticker, colors
import matplotlib.colors as colors
import netCDF4 as nc
import math
import numpy as np

'''
This script will produce some graphs in order to take a quick look at the results from 
the iterative_sliding.py script (for now produces some results as a .txt file)
This will plot the evolution of the sliding coefficients as well as the slopes obtained using 
the Greve et al., 2020 (DOI: 10.5281/zenodo.3971232) method.

This script should be run from the /sicopolis/tools/diverse directory, through the command
python3 quick_graph.py.

Make sure you have python 3 and the netCDF4 library installed through the command pip install netcdf4

Will save the files in the results folder created during the iterative script
'''

# NO EXTENSIONS IN THE FOLLOWING PARAMETERS
run_name = 'ant32_bm3_jare_aq1_spinup03_fixtopo_1'           # Name of the run inputted in the iterative python script
coef = 0.6                                                   # Coef used in the aforementioned script
name_topo = 'ant_bm3_32_topo'                                # Name of the topography file
region = 'ant'    					     # Name of the region (the folder for sico_in)
region_file_name = 'ant32_imbie2016_basins_extrapolated'     # Name of the region file (e.g. the imbie basins)
run_name_base_run = 'ant32_bm3_jare_aq1_spinup03_fixtopo'    # Name of the base run with 1m/Pa-1 everywhere
num_output_base_run = '0002'                                 # Output number of the base simulation with 1m/Pa-1 everywhere


def avg_sliding_slope(n):

    path = f'../../sico_in/{region}/'
    filename = f'{region_file_name}.nc'

    regions = nc.Dataset(path+filename)
    n_reg = regions['n_basin'][:]
    
    if k == 0:
       path = f'../../sico_out/{run_name_base_run}/'
       filename = f'{run_name_base_run}{num_output_base_run}.nc'  
    else:       
        path = f'../../sico_out/{run_name}_{coef}_{n}/'
        filename = f'{run_name}_{coef}_{n}0001.nc'

    sim = nc.Dataset(path+filename)
    mask = sim['mask'][:]

    nmax = len(mask)
    ntot_ground = 0
    frac_region = []
    for n_cnt in range(1, max(max(sub) for sub in n_reg)+1):

        ntot_reg = 0
        for i in range(nmax):
            for j in range(nmax):
                if (n_reg[i][j] == n_cnt and mask[i][j] == 0):
                    ntot_reg += 1
                if (n_cnt == 1 and mask[i][j] == 0):            
                    ntot_ground +=1
        frac_region.append(ntot_reg)
    return (frac_region, ntot_ground)
    

txt = open(f'../../my_results_{run_name}/{coef}_slope_log.txt', 'r')

content = txt.readlines()
txt.close()
for k in range(len(content)):
    content[k] = content[k].split(',')
    content[k].pop()
    content[k].pop()
    content[k] = [round(float(content[k][j]),4) for j in range(len(content[k]))]   

X = []
Y = []
num_reg = len(content[0])
num_iter = len(content)
slide = [[] for j in range(len(content[0]))]
print(num_iter)
print('Calculating the average slope for each iteration... \n')
for k in range(num_iter):
    avg = 0
    X.append(k)
    (frac_region, ntot_ground) = avg_sliding_slope(k)
    for j in range(len(content[0])):
        slide[j].append(content[k][j])
        avg = avg + frac_region[j] * content[k][j] / ntot_ground
    Y.append(avg)
    print('Iteration',k,':', avg)



plt.figure()
plt.title(f'Evolution of the slopes against iteration \n', fontsize = 14)
plt.xlabel('Number of iterations', fontsize = 14)
plt.ylabel('Slope', fontsize = 14)
plt.tight_layout(rect=[0, 0, 0.85, 1])
fig = plt.gcf()
for k in range(len(slide)):
    plt.plot(X, slide[k], label=f'{k}')
plt.plot(X, Y, label = 'Avg', linewidth=6)
plt.plot(X, [1 for k in range(len(Y))], 'k--')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig(f'../../my_results_{run_name}/slope_evol_{coef}_{run_name}.png', dpi=600)
plt.close('all')

C = []
for k in range(1,num_iter):
    header = open(f'../../sico_out/{run_name}_{coef}_{k}/sico_specs_{run_name}_{coef}_{k}.h', 'r')
    content = header.readlines()
    header.close()
    a = content[1117]
    a = a.replace('#define C_SLIDE (/ ', '')
    a = a.replace(' ', '')
    a = a.replace('\\', '')
    a = a.replace('d0', '')
    a = a.replace('/)', '')
    a = a.split(',')

    combi = [float(a[k]) for k in range(len(a))]
    C.append(combi)

X = []
Y = []

slide = [[] for j in range(len(C[0]))]
for k in range(len(C)):
    X.append(k)
    for j in range(len(C[0])):
        slide[j].append(C[k][j])

X.append(len(C))
plt.figure()
plt.title('Evolution of the sliding coefficient against iteration\n', fontsize = 14)
plt.xlabel('Number of iterations', fontsize = 14)
plt.ylabel('Sliding coefficients m/(a*Pa)', fontsize = 14)
plt.tight_layout(rect=[0, 0, 0.85, 1])
fig = plt.gcf()
for k in range(len(slide)):
    plt.plot(X, [1] + slide[k], label = f'{k+1}')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig(f'../../my_results_{run_name}/c_slide_evol_{run_name}_{coef}.png', dpi=600)
plt.close('all')


def create_coeff_map():
    coeff_map = []
    path = f'../../sico_in/{region}/'
    filename = f'{region_file_name}.nc'

    regions = nc.Dataset(path+filename)
    n_reg = regions['n_basin'][:]

    path = f'../../sico_out/{run_name}_{coef}_{num_iter-1}/'
    filename = f'{run_name}_{coef}_{num_iter-1}0001.nc'

    sim = nc.Dataset(path+filename)
    mask = sim['mask'][:]
    x = sim['x'][:]
    y = sim['y'][:]

    nmax = len(mask)
    coeff_map = [[math.nan for k in mask] for j in mask]
    for n_cnt in range(1, max(max(sub) for sub in n_reg)+1):
        for i in range(nmax):
            for j in range(nmax):
                if (mask[i][j] == 0 and n_reg[i][j] == n_cnt):
                    coeff_map[i][j] = slide[n_cnt-1][len(slide[0])-1]
    regions_map = []
    for i in range(nmax):
        regions_map.append([])
        for j in range(nmax):
            if (mask[i][j] == 0):
                regions_map[i].append(n_reg[i][j]*10)
            else :
                regions_map[i].append(math.nan)
    return regions_map, coeff_map, x, y

def create_plot(regions_map, coeff_map,x,y):

    x = [k/1000 for k in x]
    y = [k/1000 for k in y]

    regions = nc.Dataset(f'../../sico_in/{region}/{name_topo}.nc')
    H = regions['H'][:]
    mask = regions['mask'][:]
    H[mask != 2] +=10
    H[mask == 3] = 0

    cmap = colors.ListedColormap(['white'] + plt.cm.viridis.colors)

    plt.figure()
    plt.gca().set_aspect('equal')
    a = plt.pcolormesh(y, x, H, cmap=cmap) 
    label = f'{num_reg} regions\n'
    plt.title(label, fontsize=14)
    plt.ylabel('y (km)', fontsize= 14)
    plt.xlabel('x (km)', fontsize= 14)
    plt.xticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
    plt.yticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
    plt.contour(x,y, regions_map, levels = num_reg, colors = 'k')
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(f'../../my_results_{run_name}/regions_map.png', dpi=600)
    plt.close('all')

    levels = np.linspace(np.nanmin(regions_map), np.nanmax(regions_map), 19)
    plt.figure()

    regions_map_array = np.array(regions_map)
    regions_done = []
    for i in range(len(regions_map)):
        for j in range(len(regions_map[0])):
            if not np.isnan(regions_map[i][j]):
                if regions_map[i][j] not in regions_done:
                    coeff_indices = np.where(regions_map_array == regions_map[i][j])
                    if coeff_indices[0].size > 0:
                        coeff = coeff_map[i][j]
                        center_i = int(np.median(coeff_indices[0]))
                        center_j = int(np.median(coeff_indices[1]))
                        plt.text(y[center_j], x[center_i], f'{coeff:.2f}', ha='center', va='center', fontsize=8)
                        regions_done.append(regions_map[i][j])

    coeff_map_array = np.array(coeff_map)
    a = plt.pcolormesh(y, x, coeff_map_array, cmap='cool') 
    plt.gca().set_aspect('equal')
    label = f'Sliding coefficient in each region for \n {run_name}_{num_iter} \n (m/[a*Pa^(p-q)])'
    plt.title(label, fontsize=14)
    plt.ylabel('y (km)', fontsize= 14)
    plt.xlabel('x (km)', fontsize= 14)
    plt.xticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
    plt.yticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
    plt.contour(x,y, regions_map,levels =levels, colors = 'k', corner_mask=True)
    plt.tight_layout()
    fig = plt.gcf()
    fig.savefig(f'../../my_results_{run_name}/coeff_map_{run_name}_{coef}.png', dpi=600)
    plt.close('all')


regions, coeff_map, x, y = create_coeff_map()
create_plot(regions, coeff_map,x,y)

print('Job done')
