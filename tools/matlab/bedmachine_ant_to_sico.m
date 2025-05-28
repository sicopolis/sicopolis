%==========================================================================
% bedmachine_ant_to_sico
%
% Description:
%   Reading of the BedMachine data (currently version 3)
%   for the Antarctic ice sheet,
%   downsampling them on the resolutions used by SICOPOLIS,
%   writing them on NetCDF files in the format required by SICOPOLIS.
%
% Author: Ralf Greve
% Date:   2025-05-02
%==========================================================================

clear variables
close all

%-------- Parameter settings --------

bm_version = '3';   % Version number of BedMachine data (string)
% bm_version = '3p8';   %%% (not officially released version)

inpath   = ['/uchi/greve/sicopolis/data/ant_bedmachine/v' bm_version];
filename = ['BedMachineAntarctica-v' bm_version '.nc'];

A = 6378137.0;
%   WGS84 semi-major axis [m]

F_INV = 298.257223563;
%   WGS84 inverse flattening

lond0 =   0.0;   % Central meridian
latd0 = -71.0;   % Standard parallel

lon0 = deg2rad(lond0);   % Central meridian
lat0 = deg2rad(latd0);   % Standard parallel

xmin = -3040e3;   % domain limits of the SICOPOLIS EPSG:3031 grid [m]
xmax =  3040e3;
ymin = -3040e3;
ymax =  3040e3;

rho_i  =  910.0;  % density of ice [kg/m3]
rho_sw = 1028.0;  % density of sea water [kg/m3]

zl_deep_ocean = -4000.0;   % sea bed of deep ocean [m]

%-------- Resolution of SICOPOLIS grid --------

disp(' ')
disp('Resolution of SICOPOLIS grid:')
disp(' ')

prompt = ['Horizontal resolution: ' ...
          '(1) 64 km, (2) 40 km, (3) 32 km, (4) 20 km, \n' ...
          '                       ' ...
          '(5) 16 km, (6) 10 km, (7)  8 km, (8)  1 km. \n' ...
          'Enter [1-8] > '];
dx_number = input(prompt);

if isempty(dx_number) ...
   || ( dx_number < 1 || dx_number > 8 )
   error('Wrong value!')
end

% disp(' ')
% disp('Treatment of Lake Vostok:')
% 
% disp(' ')
% disp('(1) grounded ice, (2) floating ice.')
% n_vostok = input('Enter [1-2] > ');
% 
% if ( (n_vostok < 1) || (n_vostok > 2) )
%    error('Wrong value!')
% end

%-------- Reading of the BedMachine data --------

disp(' ')
disp(['Reading data file ' filename ' ...'])

ncid = netcdf.open([inpath '/' filename], 'nowrite');

varid = netcdf.inqVarID(ncid, 'x');
x_bm  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'y');
y_bm  = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'bed');
zl_bm = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'surface');
zs_bm = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 'thickness');
H_bm  = netcdf.getVar(ncid, varid);

varid   = netcdf.inqVarID(ncid, 'mask');
mask_bm = netcdf.getVar(ncid, varid);

netcdf.close(ncid);

y_bm    = flip(y_bm);
zl_bm   = fliplr(zl_bm);
zs_bm   = fliplr(zs_bm);
H_bm    = fliplr(H_bm);
mask_bm = fliplr(mask_bm);

imax_bm = length(x_bm)-1;
jmax_bm = length(y_bm)-1;

xmin_bm = double(x_bm(1));
xmax_bm = double(x_bm(end));
ymin_bm = double(y_bm(1));
ymax_bm = double(y_bm(end));

dx_bm   = double(x_bm(2)-x_bm(1));

x_bm    = double(x_bm);
y_bm    = double(y_bm);
zl_bm   = double(zl_bm);
zs_bm   = double(zs_bm);
H_bm    = double(H_bm);
mask_bm = int16(mask_bm);

%  ------ Auxiliary masks

mask_ocean_bm                = zeros(size(mask_bm), 'int16');
mask_ocean_bm(mask_bm==0)    = 1;

mask_land_bm                 = zeros(size(mask_bm), 'int16');
mask_land_bm(mask_bm==1)     = 1;

mask_grounded_bm             = zeros(size(mask_bm), 'int16');
mask_grounded_bm(mask_bm==2) = 1;

mask_floating_bm             = zeros(size(mask_bm), 'int16');
mask_floating_bm(mask_bm==3) = 1;

mask_lakevstk_bm             = zeros(size(mask_bm), 'int16');
mask_lakevstk_bm(mask_bm==4) = 1;      

% if n_vostok == 1
   H_bm(mask_lakevstk_bm==1) = zs_bm(mask_lakevstk_bm==1) ...
                                  - zl_bm(mask_lakevstk_bm==1);
   mask_grounded_bm(mask_lakevstk_bm==1) = 1;
                               % treat Lake Vostok like grounded ice
% else
%    mask_floating_bm(mask_lakevstk_bm==1) = 1;
%                                % treat Lake Vostok like floating ice
% end

%-------- Definition of SICOPOLIS grid --------

disp(' ')
disp('Defining SICOPOLIS grid ...')

if (dx_number == 1)
   dx    =  64e3;
   ch_dx = '64';
elseif (dx_number == 2)
   dx    =  40e3;
   ch_dx = '40';
elseif (dx_number == 3)
   dx    =  32e3;
   ch_dx = '32';
elseif (dx_number == 4)
   dx    =  20e3;
   ch_dx = '20';
elseif (dx_number == 5)
   dx    =  16e3;
   ch_dx = '16';
elseif (dx_number == 6)
   dx    =  10e3;
   ch_dx = '10';
elseif (dx_number == 7)
   dx    =   8e3;
   ch_dx = '08';
elseif (dx_number == 8)
   dx    =   1e3;
   ch_dx = '01';
end

disp(' ')
disp(['Resolution = ' ch_dx ' km ...'])

imax = round((xmax-xmin)/dx);
jmax = round((ymax-ymin)/dx);

x = xmin:dx:xmax;
y = ymin:dx:ymax;

%-------- Downsampling of BedMachine data on SICOPOLIS grid --------

disp(' ')
disp('Downsampling BedMachine data (takes a while) ...')

zl  = zeros(imax+1,jmax+1);
zs  = zeros(imax+1,jmax+1);
H   = zeros(imax+1,jmax+1);

r_mask_grounded = zeros(imax+1,jmax+1);
r_mask_floating = zeros(imax+1,jmax+1);
r_mask_land     = zeros(imax+1,jmax+1);
r_mask_ocean    = zeros(imax+1,jmax+1);

for i=1:imax+1
for j=1:jmax+1

   if (x(i) < xmin_bm || x(i) > xmax_bm) ...
      || ...
      (y(j) < ymin_bm || y(j) > ymax_bm)

       zl(i,j)           = zl_deep_ocean;
       r_mask_ocean(i,j) = 1.0;

   else
        
       x1 = x(i) - 0.5*(dx+dx_bm);
       x2 = x(i) + 0.5*(dx+dx_bm);
       y1 = y(j) - 0.5*(dx+dx_bm);
       y2 = y(j) + 0.5*(dx+dx_bm);

       ii1 = floor((x1-xmin_bm)/dx_bm + 1);
       ii2 =  ceil((x2-xmin_bm)/dx_bm + 1);
       jj1 = floor((y1-ymin_bm)/dx_bm + 1);
       jj2 =  ceil((y2-ymin_bm)/dx_bm + 1);

       sum_weight = 0.0;
        
       for ii=ii1:ii2
       for jj=jj1:jj2

           if (ii>=1 && ii<=(imax_bm+1)) ...
              && ...
              (jj>=1 && jj<=(jmax_bm+1))

               weight = frac_area([x(i)-0.5*dx x(i)+0.5*dx], ...
                                  [y(j)-0.5*dx y(j)+0.5*dx], ...
                          [x_bm(ii)-0.5*dx_bm x_bm(ii)+0.5*dx_bm], ...
                          [y_bm(jj)-0.5*dx_bm y_bm(jj)+0.5*dx_bm]);
                           
               sum_weight = sum_weight + weight;
               zl(i,j)  = zl(i,j) + weight*zl_bm(ii,jj);
               zs(i,j)  = zs(i,j) + weight*zs_bm(ii,jj);
               H(i,j)   =  H(i,j) + weight*H_bm(ii,jj) ;
                
               r_mask_grounded(i,j) = r_mask_grounded(i,j) ...
                               + weight*double(mask_grounded_bm(ii,jj));
               r_mask_floating(i,j) = r_mask_floating(i,j) ...
                               + weight*double(mask_floating_bm(ii,jj));
               r_mask_land(i,j)     = r_mask_land(i,j) ...
                               + weight*double(mask_land_bm(ii,jj));
               r_mask_ocean(i,j)    = r_mask_ocean(i,j) ...
                               + weight*double(mask_ocean_bm(ii,jj));

           end
            
       end
       end
        
       zl(i,j) = zl(i,j)  / sum_weight;
       zs(i,j) = zs(i,j)  / sum_weight;
       H(i,j)  = H(i,j)   / sum_weight;
        
       r_mask_grounded(i,j) = r_mask_grounded(i,j) / sum_weight;
       r_mask_floating(i,j) = r_mask_floating(i,j) / sum_weight;
       r_mask_land(i,j)     = r_mask_land(i,j)     / sum_weight;
       r_mask_ocean(i,j)    = r_mask_ocean(i,j)    / sum_weight;
        
    end
        
end
end

%-------- Ice-land-ocean mask --------

disp(' ')
disp('Creating ice-land-ocean mask ...')

mask = zeros(size(zs), 'int16');

for i=1:imax+1
for j=1:jmax+1
    
    [~,n] = max([r_mask_grounded(i,j) ...
                 r_mask_floating(i,j) ...
                 r_mask_land(i,j) ...
                 r_mask_ocean(i,j)]);

    switch n
        case 1
            mask(i,j) = 0;
        case 2
            mask(i,j) = 3;
        case 3
            mask(i,j) = 1;
        case 4
            mask(i,j) = 2;
        otherwise
            error('Mask value could not be determined (NaN?)!')
    end

end
end

%  ------ Enforce connectivity of the ocean

mask_c = zeros(imax+1,jmax+1,'int16');

mask_c(:         , 1:2      ) = 1;
mask_c(:         , end-1:end) = 1;
mask_c(1:2       , :        ) = 1;
mask_c(end-1:end , :        ) = 1;

flag_change = true;

while flag_change

    mask_c_save = mask_c;

    for i=2:imax
    for j=2:jmax

        if mask_c_save(i,j) == 1
            if mask(i+1,j  ) >= 2; mask_c(i+1,j  ) = 1; end
            if mask(i-1,j  ) >= 2; mask_c(i-1,j  ) = 1; end
            if mask(i  ,j+1) >= 2; mask_c(i  ,j+1) = 1; end
            if mask(i  ,j-1) >= 2; mask_c(i  ,j-1) = 1; end
            if mask(i+1,j+1) >= 2; mask_c(i+1,j+1) = 1; end
            if mask(i-1,j+1) >= 2; mask_c(i-1,j+1) = 1; end
            if mask(i+1,j-1) >= 2; mask_c(i+1,j-1) = 1; end
            if mask(i-1,j-1) >= 2; mask_c(i-1,j-1) = 1; end
        end        

    end
    end

    mask_c_diff = abs(mask_c-mask_c_save);

    if max(mask_c_diff(:)) > 0
        flag_change = true;
    else
        flag_change = false;
    end

end

for i=1:imax+1
for j=1:jmax+1

    if (mask(i,j) == 2) && (mask_c(i,j) == 0)
        mask(i,j) = 1;
    end

    if (mask(i,j) == 3) && (mask_c(i,j) == 0)
        mask(i,j) = 0;
    end

    if (mask(i,j) <= 1) && (mask_c(i,j) == 1) && (zs(i,j) < 0.0)
        mask(i,j) = 2;
    end

end
end

%-------- Ice base, consistency corrections --------

disp(' ')
disp('Ice base, consistency corrections ...')

mask_corr = mask;

H  = max(H,  0.0);
zs = max(zs, 0.0);

H_offset = 10.0;

zb   = zeros(size(zs));

for i=1:imax+1
for j=1:jmax+1
    
    switch mask(i,j)
        case 0
            zb(i,j) = zl(i,j);
            zs(i,j) = max(zs(i,j), zb(i,j)+H_offset);
            H(i,j)  = zs(i,j) - zb(i,j);
        case 1
            zb(i,j) = zl(i,j);
            zs(i,j) = zl(i,j);
            H(i,j)  = 0.0;
        case 2
            if zl(i,j) <= 0.0
                zs(i,j) = 0.0;
                zb(i,j) = 0.0;
                H(i,j)  = 0.0;
            else   % sea bed above sea level -> correction
                mask_corr(i,j) = 1;
                zb(i,j) = zl(i,j);
                zs(i,j) = zl(i,j);
                H(i,j)  = 0.0;
            end
        case 3
            zs(i,j) = max(zs(i,j), H_offset*(rho_sw-rho_i)/rho_sw);
            zb(i,j) = -rho_i/(rho_sw-rho_i) * zs(i,j);
            H(i,j)  = zs(i,j) - zb(i,j);
            if zb(i,j) < zl(i,j)   % ice base below sea bed -> correction
                mask_corr(i,j) = 0;
                zb(i,j) = zl(i,j);
                zs(i,j) = max(zs(i,j), zb(i,j)+H_offset);
                H(i,j)  = zs(i,j) - zb(i,j);
            end
        otherwise
            error('Mask value not valid!')
    end

end
end

mask = mask_corr;

if min(H(:)) < 0.0
    error('Negative ice thickness!')
end

%  ------ Enforce connectivity of the ocean yet again

mask_c(:) = 0;

mask_c(:         , 1:2      ) = 1;
mask_c(:         , end-1:end) = 1;
mask_c(1:2       , :        ) = 1;
mask_c(end-1:end , :        ) = 1;

flag_change = true;

while flag_change

    mask_c_save = mask_c;

    for i=2:imax
    for j=2:jmax

        if mask_c_save(i,j) == 1

            if mask(i+1,j  ) >= 2; mask_c(i+1,j  ) = 1; end
            if mask(i-1,j  ) >= 2; mask_c(i-1,j  ) = 1; end
            if mask(i  ,j+1) >= 2; mask_c(i  ,j+1) = 1; end
            if mask(i  ,j-1) >= 2; mask_c(i  ,j-1) = 1; end
            if mask(i+1,j+1) >= 2; mask_c(i+1,j+1) = 1; end
            if mask(i-1,j+1) >= 2; mask_c(i-1,j+1) = 1; end
            if mask(i+1,j-1) >= 2; mask_c(i+1,j-1) = 1; end
            if mask(i-1,j-1) >= 2; mask_c(i-1,j-1) = 1; end

            if (mask(i+1,j  ) <= 1) && (zs(i+1,j  ) < 0.0)
                                    mask_c(i+1,j  ) = 1; end
            if (mask(i-1,j  ) <= 1) && (zs(i-1,j  ) < 0.0)
                                    mask_c(i-1,j  ) = 1; end
            if (mask(i  ,j+1) <= 1) && (zs(i  ,j+1) < 0.0)
                                    mask_c(i  ,j+1) = 1; end
            if (mask(i  ,j-1) <= 1) && (zs(i  ,j-1) < 0.0)
                                    mask_c(i  ,j-1) = 1; end
            if (mask(i+1,j+1) <= 1) && (zs(i+1,j+1) < 0.0)
                                    mask_c(i+1,j+1) = 1; end
            if (mask(i-1,j+1) <= 1) && (zs(i-1,j+1) < 0.0)
                                    mask_c(i-1,j+1) = 1; end
            if (mask(i+1,j-1) <= 1) && (zs(i+1,j-1) < 0.0)
                                    mask_c(i+1,j-1) = 1; end
            if (mask(i-1,j-1) <= 1) && (zs(i-1,j-1) < 0.0)
                                    mask_c(i-1,j-1) = 1; end
            
        end        

    end
    end

    mask_c_diff = abs(mask_c-mask_c_save);

    if max(mask_c_diff(:)) > 0
        flag_change = true;
    else
        flag_change = false;
    end

end

for i=1:imax+1
for j=1:jmax+1

    if (mask(i,j) == 2) && (mask_c(i,j) == 0)
        mask(i,j) = 1;
        zb(i,j)   = zl(i,j);
        zs(i,j)   = zl(i,j);
        H(i,j)    = 0.0;
    end

    if (mask(i,j) == 3) && (mask_c(i,j) == 0)
        mask(i,j) = 0;
        zb(i,j)   = zl(i,j);
        H(i,j)    = zs(i,j)-zb(i,j);
    end

    if (mask(i,j) <= 1) && (mask_c(i,j) == 1) && (zs(i,j) < 0.0)
        mask(i,j) = 2;
        zs(i,j)   = 0.0;
        zb(i,j)   = 0.0;
        H(i,j)    = 0.0;
    end

end
end

%-------- Isostatically relaxed lithosphere surface
%                               (rigid bed assumed) --------

zl0 = zl;

%-------- Maximum-extent mask --------

mask_maxextent = zeros(size(mask), 'int16');

mask_maxextent(mask==0 | mask==3) = 1;

%  ------ Connectivity of the points not allowed to glaciate

mask_c(:) = 0;

mask_c(:         , 1:2      ) = 1;
mask_c(:         , end-1:end) = 1;
mask_c(1:2       , :        ) = 1;
mask_c(end-1:end , :        ) = 1;

flag_change = true;

while flag_change

    mask_c_save = mask_c;

    for i=2:imax
    for j=2:jmax

        if mask_c_save(i,j) == 1

            if mask_maxextent(i+1,j  ) == 0; mask_c(i+1,j  ) = 1; end
            if mask_maxextent(i-1,j  ) == 0; mask_c(i-1,j  ) = 1; end
            if mask_maxextent(i  ,j+1) == 0; mask_c(i  ,j+1) = 1; end
            if mask_maxextent(i  ,j-1) == 0; mask_c(i  ,j-1) = 1; end
            if mask_maxextent(i+1,j+1) == 0; mask_c(i+1,j+1) = 1; end
            if mask_maxextent(i-1,j+1) == 0; mask_c(i-1,j+1) = 1; end
            if mask_maxextent(i+1,j-1) == 0; mask_c(i+1,j-1) = 1; end
            if mask_maxextent(i-1,j-1) == 0; mask_c(i-1,j-1) = 1; end

        end        

    end
    end

    mask_c_diff = abs(mask_c-mask_c_save);

    if max(mask_c_diff(:)) > 0
        flag_change = true;
    else
        flag_change = false;
    end

end

for i=1:imax+1
for j=1:jmax+1

    if (mask_maxextent(i,j) == 0) && (mask_c(i,j) == 0)
        mask_maxextent(i,j) = 1;
    end

end
end

%-------- Writing of data on file --------

disp(' ')
disp('Writing data on file ...')

filename = ['ant_bm' bm_version '_' ch_dx '_topo.nc'];

nccreate(filename, ...
         'mapping', ...
         'Datatype', 'int8', ...
         'Format', 'classic')

ncwrite(filename, 'mapping', int8(0))
ncwriteatt(filename, 'mapping', 'grid_mapping_name', 'polar_stereographic')
ncwriteatt(filename, 'mapping', 'reference_ellipsoid_name', 'WGS84')
ncwriteatt(filename, 'mapping', 'semi_major_axis', A)
ncwriteatt(filename, 'mapping', 'inverse_flattening', F_INV)
ncwriteatt(filename, 'mapping', 'latitude_of_projection_origin', -90.0)
ncwriteatt(filename, 'mapping', 'standard_parallel', latd0)
ncwriteatt(filename, 'mapping', ...
                     'straight_vertical_longitude_from_pole', lond0)
ncwriteatt(filename, 'mapping', 'false_easting', 0.0)
ncwriteatt(filename, 'mapping', 'false_northing', 0.0)

nccreate(filename, ...
         'x', ...
         'Dimensions', {'x', imax+1}, ...
         'Datatype', 'double', ...
         'Format', 'classic')

ncwrite(filename, 'x', x)
ncwriteatt(filename, 'x', 'standard_name', 'projection_x_coordinate')
ncwriteatt(filename, 'x', 'units', 'm')
ncwriteatt(filename, 'x', 'axis', 'x')

nccreate(filename, ...
         'y', ...
         'Dimensions', {'y', jmax+1}, ...
         'Datatype', 'double', ...
         'Format', 'classic')

ncwrite(filename, 'y', y)
ncwriteatt(filename, 'y', 'standard_name', 'projection_y_coordinate')
ncwriteatt(filename, 'y', 'units', 'm')
ncwriteatt(filename, 'y', 'axis', 'y')

% nccreate(filename, ...
%          'lon', ...
%          'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
%          'Datatype', 'single', ...
%          'Format', 'classic')
% 
% ncwrite(filename, 'lon', lond)
% ncwriteatt(filename, 'lon', 'standard_name', 'longitude')
% ncwriteatt(filename, 'lon', 'units', 'deg')
% ncwriteatt(filename, 'lon', 'grid_mapping', 'mapping')

% nccreate(filename, ...
%          'lat', ...
%          'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
%          'Datatype', 'single', ...
%          'Format', 'classic')
% 
% ncwrite(filename, 'lat', latd)
% ncwriteatt(filename, 'lat', 'standard_name', 'latitude')
% ncwriteatt(filename, 'lat', 'units', 'deg')
% ncwriteatt(filename, 'lat', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'zs', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'single', ...
         'Format', 'classic')

ncwrite(filename, 'zs', zs)
ncwriteatt(filename, 'zs', 'standard_name', 'surface_altitude')
ncwriteatt(filename, 'zs', 'units', 'm')
% ncwriteatt(filename, 'zs', '_FillValue', topo_no_value)
% ncwriteatt(filename, 'zs', 'missing_value', topo_no_value)
ncwriteatt(filename, 'zs', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'zb', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'single', ...
         'Format', 'classic')

ncwrite(filename, 'zb', zb)
ncwriteatt(filename, 'zb', 'standard_name', 'ice_base_altitude')
ncwriteatt(filename, 'zb', 'units', 'm')
% ncwriteatt(filename, 'zb', '_FillValue', topo_no_value)
% ncwriteatt(filename, 'zb', 'missing_value', topo_no_value)
ncwriteatt(filename, 'zb', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'zl', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'single', ...
         'Format', 'classic')

ncwrite(filename, 'zl', zl)
ncwriteatt(filename, 'zl', 'standard_name', 'bedrock_altitude')
ncwriteatt(filename, 'zl', 'units', 'm')
% ncwriteatt(filename, 'zl', '_FillValue', topo_no_value)
% ncwriteatt(filename, 'zl', 'missing_value', topo_no_value)
ncwriteatt(filename, 'zl', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'zl0', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'single', ...
         'Format', 'classic')

ncwrite(filename, 'zl0', zl0)
ncwriteatt(filename, 'zl0', 'standard_name', ...
                            'isostatically_relaxed_bedrock_altitude')
ncwriteatt(filename, 'zl0', 'units', 'm')
% ncwriteatt(filename, 'zl0', '_FillValue', topo_no_value)
% ncwriteatt(filename, 'zl0', 'missing_value', topo_no_value)
ncwriteatt(filename, 'zl0', 'grid_mapping', 'mapping')
ch_att = 'Rigid lithosphere assumed, thus zl and zl0 are identical';
ncwriteatt(filename, 'zl0', 'note', ch_att)

nccreate(filename, ...
         'H', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'single', ...
         'Format', 'classic')

ncwrite(filename, 'H', H)
ncwriteatt(filename, 'H', 'standard_name', ...
                            'land_ice_thickness')
ncwriteatt(filename, 'H', 'units', 'm')
% ncwriteatt(filename, 'H', '_FillValue', topo_no_value)
% ncwriteatt(filename, 'H', 'missing_value', topo_no_value)
ncwriteatt(filename, 'H', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'mask', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'int8', ...
         'Format', 'classic')

ncwrite(filename, 'mask', mask)
ncwriteatt(filename, 'mask', 'standard_name', 'ice_land_sea_mask')
ncwriteatt(filename, 'mask', 'flag_values', int8([0 1 2 3]))
ncwriteatt(filename, 'mask', 'flag_meanings', ...
                     ['grounded_ice ' 'ice_free_land '...
                      'sea ' 'floating_ice'])
ncwriteatt(filename, 'mask', 'grid_mapping', 'mapping')

nccreate(filename, ...
         'mask_maxextent', ...
         'Dimensions', {'x', imax+1, 'y', jmax+1}, ...
         'Datatype', 'int8', ...
         'Format', 'classic')

ncwrite(filename, 'mask_maxextent', mask_maxextent)
ncwriteatt(filename, 'mask_maxextent', 'standard_name', 'max_extent_mask')
ncwriteatt(filename, 'mask_maxextent', 'flag_values', int8([0 1]))
ncwriteatt(filename, 'mask_maxextent', 'flag_meanings', ...
                     ['not_allowed_to_glaciate ' 'allowed_to_glaciate'])
ncwriteatt(filename, 'mask_maxextent', 'grid_mapping', 'mapping')

%-------- End of script --------

disp(' ')
disp('Done.')
disp(' ')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
