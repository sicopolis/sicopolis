import matplotlib
matplotlib.use('Agg')  # Forces Matplotlib to not look for a screen

from statistics import mean
from matplotlib import pyplot as plt
from matplotlib import ticker, colors
import matplotlib.colors as colors
import netCDF4 as nc
import math
import numpy as np

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
This script will produce some graphs in order to take a quick look at
the results from the iterative_sliding.py script (for now produces
some results as a .txt file).
This will plot the evolution of the sliding coefficients as well as
the slopes obtained using the Greve et al. 2020
(DOI: 10.5281/zenodo.3971232) method.

Change to the main SICOPOLIS directory, and make a local copy
of the template:
  cp ./tools/diverse/iterative_sliding_quick_graph.py \
     ./my_iterative_sliding_quick_graph_xxx.py
(where 'xxx' can be any descriptor of the user's choice).

Execution of the script:
  python3.11 my_iterative_sliding_quick_graph_xxx.py
(A newer version of Python will also do.)

Make sure you have python 3 and the netCDF4 library installed through
the command pip install netcdf4.

Will save the files in the results folder created during the iterative script.
'''

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------- Variables to fill --------

domain = 'grl'
# Name of the computational domain
# ('ant' for the Antarctic ice sheet, 'grl' for the Greenland ice sheet)

dx = '16'  # Resolution (string)

name_of_run = f'{domain}{dx}_bm6_spinup11_cal_100ka_iter'
# Name of the 'startup' header ('0th' iteration)

kmax = 7  # Maximum number of iterations (typically 5-15)

modifier = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# Values of the modifier (relaxation factor) for the iterations;
# typically 0.5-1.0 (the smaller, the less aggressive).
# Length of array must be kmax+1 (k=0...kmax).

if len(modifier) != kmax+1:
    print(f'Length of \'modifier\' array does not match value of \'kmax\'!')
    exit()

topo_file_name = f'{domain}_bm6_{dx}_topo'
# Name of the topography file (without extension '.nc')

regions_file_name = f'{domain}{dx}_zwally2012_basins_with-negis_extrapolated'
# Name of the regions file (without extension '.nc')

tslice = '0001'
# Time-slice number for final state of iterations k=0...kmax
# (which is used for comparison with observed surface velocities)

#-------- Calculating the average slope for each iteration -------

def avg_sliding_slope(n):

    path = f'./sico_in/{domain}/'
    filename = f'{regions_file_name}.nc'

    regions = nc.Dataset(path+filename)
    n_reg = regions['n_basin'][:]
    
    path = f'./sico_out/{name_of_run}_{n:02d}_{modifier[n]}/'
    filename = f'{name_of_run}_{n:02d}_{modifier[n]}{tslice}.nc'

    sim = nc.Dataset(path+filename)
    mask = sim['mask'][:]

    imax = len(mask)
    jmax = len(mask[0])

    ntot_ground = 0
    frac_region = []
    for n_cnt in range(1, max(max(sub) for sub in n_reg)+1):

        ntot_reg = 0
        for i in range(imax):
            for j in range(jmax):
                if (n_reg[i][j] == n_cnt and mask[i][j] == 0):
                    ntot_reg += 1
                if (n_cnt == 1 and mask[i][j] == 0):            
                    ntot_ground +=1
        frac_region.append(ntot_reg)
    return (frac_region, ntot_ground)
    
txt = open(f'./tmp/iterative_sliding_{name_of_run}/slope_log.txt', 'r')

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

print(f'num_reg  ={num_reg:3d}')
print(f'num_iter ={num_iter:3d}')
print()

if num_iter-1 != kmax:
    print(f'Detected value of \'num_iter\' does not match set value of \'kmax\'!')
    exit()

slide = [[] for j in range(num_reg)]

print('Calculating the average slope for each iteration... \n')

for k in range(num_iter):
    avg = 0
    X.append(k)
    (frac_region, ntot_ground) = avg_sliding_slope(k)
    for j in range(num_reg):
        slide[j].append(content[k][j])
        avg = avg + frac_region[j] * content[k][j] / ntot_ground
    Y.append(avg)
    print('Iteration',k,':', avg)

#-------- Plotting the slopes --------

print()
print('Making plots... \n')

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
fig.savefig(f'./tmp/iterative_sliding_{name_of_run}/slope_evol_{name_of_run}.png', dpi=600)
plt.close('all')

#-------- Plotting the sliding coefficients --------

C = []

for k in range(1,num_iter):
    header = open(f'./sico_out/{name_of_run}_{k:02d}_{modifier[k]}/sico_specs_{name_of_run}_{k:02d}_{modifier[k]}.h', 'r')
    content = header.readlines()
    header.close()

    line_num = None
    for num, line in enumerate(content):
        if '#define C_SLIDE_DIMLESS' in line:
            line_num = num
            break
    if line_num is None:
        print(f'Error: \'C_SLIDE_DIMLESS\' was not found in the header file!')
        exit()

    a = content[line_num]
    a = a.replace('#define C_SLIDE_DIMLESS [', '')
    a = a.replace(' ', '')
    a = a.replace('\\', '')
    a = a.replace('d0', '')
    a = a.replace(']', '')
    a = a.replace('\n', '')
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
plt.ylabel('Sliding coefficients', fontsize = 14)
plt.tight_layout(rect=[0, 0, 0.85, 1])

fig = plt.gcf()

for k in range(len(slide)):
    plt.plot(X, [1] + slide[k], label = f'{k+1}')
    # !!! This needs to be changed !!!

plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig(f'./tmp/iterative_sliding_{name_of_run}/c_slide_evol_{name_of_run}.png', dpi=600)
plt.close('all')

# def create_coeff_map():
# 
#     coeff_map = []
#     path = f'./sico_in/{domain}/'
#     filename = f'{regions_file_name}.nc'
# 
#     regions = nc.Dataset(path+filename)
#     n_reg = regions['n_basin'][:]
# 
#     path = f'./sico_out/{name_of_run}_{num_iter-1:02d}_{modifier[num_iter-1]}/'
#     filename = f'{name_of_run}_{num_iter-1:02d}_{modifier[num_iter-1]}{tslice}.nc'
# 
#     sim = nc.Dataset(path+filename)
#     mask = sim['mask'][:]
#     x = sim['x'][:]
#     y = sim['y'][:]
# 
#     imax = len(mask)
#     jmax = len(mask[0])
# 
#     coeff_map = [[math.nan for k in mask] for j in mask]
#     for n_cnt in range(1, max(max(sub) for sub in n_reg)+1):
#         for i in range(imax):
#             for j in range(jmax):
#                 if (mask[i][j] == 0 and n_reg[i][j] == n_cnt):
#                     coeff_map[i][j] = slide[n_cnt-1][len(slide[0])-1]
#     regions_map = []
#     for i in range(imax):
#         regions_map.append([])
#         for j in range(jmax):
#             if (mask[i][j] == 0):
#                 regions_map[i].append(n_reg[i][j]*10)
#             else :
#                 regions_map[i].append(math.nan)
#     return regions_map, coeff_map, x, y
# 
# def create_plot(regions_map, coeff_map,x,y):
# 
#     x = [k/1000 for k in x]
#     y = [k/1000 for k in y]
# 
#     regions = nc.Dataset(f'./sico_in/{domain}/{topo_file_name}.nc')
#     H = regions['H'][:]
#     mask = regions['mask'][:]
#     H[mask != 2] +=10
#     H[mask == 3] = 0
# 
#     cmap = colors.ListedColormap(['white'] + plt.cm.viridis.colors)
# 
#     plt.figure()
#     plt.gca().set_aspect('equal')
#     a = plt.pcolormesh(y, x, H, cmap=cmap) 
#     label = f'{num_reg} regions\n'
#     plt.title(label, fontsize=14)
#     plt.ylabel('y (km)', fontsize= 14)
#     plt.xlabel('x (km)', fontsize= 14)
#     plt.xticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
#     plt.yticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
#     plt.contour(x,y, regions_map, levels = num_reg, colors = 'k')
#     plt.tight_layout()
#     fig = plt.gcf()
#     fig.savefig(f'./tmp/iterative_sliding_{name_of_run}/regions_map.png', dpi=600)
#     plt.close('all')
# 
#     levels = np.linspace(np.nanmin(regions_map), np.nanmax(regions_map), 19)
#     plt.figure()
# 
#     regions_map_array = np.array(regions_map)
#     regions_done = []
#     for i in range(len(regions_map)):
#         for j in range(len(regions_map[0])):
#             if not np.isnan(regions_map[i][j]):
#                 if regions_map[i][j] not in regions_done:
#                     coeff_indices = np.where(regions_map_array == regions_map[i][j])
#                     if coeff_indices[0].size > 0:
#                         coeff = coeff_map[i][j]
#                         center_i = int(np.median(coeff_indices[0]))
#                         center_j = int(np.median(coeff_indices[1]))
#                         plt.text(y[center_j], x[center_i], f'{coeff:.2f}', ha='center', va='center', fontsize=8)
#                         regions_done.append(regions_map[i][j])
# 
#     coeff_map_array = np.array(coeff_map)
#     a = plt.pcolormesh(y, x, coeff_map_array, cmap='cool') 
#     plt.gca().set_aspect('equal')
#     label = f'Sliding coefficient in each region for \n {name_of_run}_{num_iter} \n (m/[a*Pa^(p-q)])'
#     plt.title(label, fontsize=14)
#     plt.ylabel('y (km)', fontsize= 14)
#     plt.xlabel('x (km)', fontsize= 14)
#     plt.xticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
#     plt.yticks(np.arange(min(x), max(x)+1, 3040), fontsize=12)
#     plt.contour(x,y, regions_map,levels =levels, colors = 'k', corner_mask=True)
#     plt.tight_layout()
#     fig = plt.gcf()
#     fig.savefig(f'./tmp/iterative_sliding_{name_of_run}/coeff_map_{name_of_run}.png', dpi=600)
#     plt.close('all')
# 
# regions, coeff_map, x, y = create_coeff_map()
# create_plot(regions, coeff_map,x,y)

#-------- End of script --------

print('Job done.')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
