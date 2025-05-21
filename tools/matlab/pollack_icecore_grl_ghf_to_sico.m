%==========================================================================
% pollack_icecore_grl_ghf_to_sico
%
% Description:
%   Interpolation of the spherical-harmonic representation of the
%   global heat flow by Pollack et al. (1993) to the stereographic
%   plane of SICOPOLIS for the Greenland ice sheet,
%   including prescribed values for the deep ice cores.
%
% Author: Ralf Greve
% Date:   2025-05-05
%==========================================================================

clear variables
close all

%-------- Parameter settings --------

VERSION = '2.0';
           % 1.0 - Greve (2005)
           % 1.1 - Greve and Herzfeld (2013)
           % 2.0 - Greve (2019)

GRID = 2;
       % 1 - Bamber grid
       % 2 - EPSG:3413 grid

if strcmp(VERSION, '1.0')
    ch_reference = ['Greve (2005, Ann. Glaciol. 42, 424-432, '...
                    'doi: 10.3189/172756405781812510)'];
elseif strcmp(VERSION, '1.1')
    ch_reference = ['Greve and Herzfeld (2013, Ann. Glaciol. 54(63), ' ...
                    '209-220, doi: 10.3189/2013AoG63A085)'];
elseif strcmp(VERSION, '2.0')
    ch_reference = ['Greve (2019, Polar Data J. 3, 22-36, ' ...
                    'doi: 10.20575/00000006)'];
else
   error('Wrong value of VERSION!')
end

if GRID == 1
    ch_grid = 'Bamber';
elseif GRID == 2
    ch_grid = 'EPSG3413';
else
   error('Wrong value of GRID!')
end             

eps = 1.0e-05 ;

lat_grip  = deg2rad( 72.567);
lon_grip  = deg2rad(-37.633);
lat_dye3  = deg2rad( 65.187);
lon_dye3  = deg2rad(-43.830);
lat_cc    = deg2rad( 77.167);
lon_cc    = deg2rad(-61.133);
lat_ngrip = deg2rad( 75.10 );
lon_ngrip = deg2rad(-42.32 );
lat_neem  = deg2rad( 77.5  );
lon_neem  = deg2rad(-50.9  );
lat_egrip = deg2rad( 75.63 );   % !!! to be confirmed !!!
lon_egrip = deg2rad(-35.99 );   % !!! to be confirmed !!!

lat_sass1 = deg2rad( 61.40);
lon_sass1 = deg2rad(-48.17);
lat_sass2 = deg2rad( 60.97);
lon_sass2 = deg2rad(-45.99);
lat_gap   = deg2rad( 67.15);
lon_gap   = deg2rad(-50.07);

ch_sites = {'GRIP', 'Dye3', 'CC', 'NGRIP'};

if strcmp(VERSION, '1.0')

    q_geo_grip  =  60.0;   % Prescibed heat fluxes
    q_geo_dye3  =  20.0;   % (in mW/m2)
    q_geo_cc    =  50.0;   % by Greve (2005)
    q_geo_ngrip = 135.0;

    ch_sites = {'GRIP', 'Dye3', 'CC', 'NGRIP'};

elseif strcmp(VERSION, '1.1')

    q_geo_grip  =  59.0;   % Prescibed heat fluxes
    q_geo_dye3  =  26.0;   % (in mW/m2)
    q_geo_cc    =  54.0;   % by Greve and Herzfeld (2013)
    q_geo_ngrip = 135.0;

    ch_sites = {'GRIP', 'Dye3', 'CC', 'NGRIP'};

elseif strcmp(VERSION, '2.0')

    q_geo_ave   =  60.0;   % Prescribed average heat flux for Greenland
                           % (-> scaling of the distribution by
                           %     Pollack et al. (1993)).
                           % 60 mW/m2 is what results from the
                           % GRIP-based scaling done for version 1.1.

    q_geo_grip  =  51.0;   % Prescibed heat fluxes
    q_geo_dye3  =  32.0;   % (in mW/m2)
    q_geo_cc    =  47.0;   % by Greve (2019)
    q_geo_ngrip = 135.0;
    q_geo_neem  =  70.0;

    q_geo_sass1 =  43.0;   % Borehole measurements
    q_geo_sass2 =  36.0;
    q_geo_gap   =  31.0;

    ch_sites = {'GRIP', 'Dye3', 'CC', 'NGRIP', 'NEEM', ...
                'SASS1', 'SASS2', 'GAP'};

end

if strcmp(VERSION, '2.0')

    % lat_negis_onset = deg2rad( 74.949);
    % lon_negis_onset = deg2rad(-37.711);
    lat_negis_onset = lat_egrip;
    lon_negis_onset = lon_egrip;
    
    rad_negis_onset =  25e3;   % (2-sigma radius, in m)
    
    q_geo_negis_onset_anom = -9999;
    % q_geo_negis_onset_anom = 970 - 91.3874;
    % q_geo_negis_onset_anom = 970 - 98.9137;
    % q_geo_negis_onset_anom = 2000 - 98.7222;
                             % Prescribed heat-flux anomaly (in mW/m2)
                             %%%%% Negative values will be ignored %%%%%
end                         

A = 6378137.0;
%   WGS84 semi-major axis [m]

F_INV = 298.257223563;
%   WGS84 inverse flattening

B = A*(1.0-1.0/F_INV);
%   WGS84 semi-minor axis [m]

if GRID == 1
    lond0 =   -39;   % Central meridian  Bamber [deg]
    latd0 =   +71;   % Standard parallel Bamber [deg]
    lon0  = deg2rad(lond0);
    lat0  = deg2rad(latd0);
    xmin  =  -800e3;
    xmax  =   700e3;   % Domain limits of the
    ymin  = -3400e3;   % SICOPOLIS Bamber grid [m]
    ymax  =  -600e3;
elseif GRID == 2
    lond0 =   -45;   % Central meridian  EPSG:3413 [deg]
    latd0 =   +70;   % Standard parallel EPSG:3413 [deg]
    lon0  = deg2rad(lond0);
    lat0  = deg2rad(latd0);
    xmin  =  -720e3;
    xmax  =   960e3;   % Domain limits of the
    ymin  = -3450e3;   % SICOPOLIS EPSG:3413 grid [m]
    ymax  =  -570e3;
end

domain       = 'Greenland';
domain_short = 'grl_rg_'  ;

%-------- Reading of harmonic coefficients from file --------

disp(' ')
disp('Reading harmonic coefficients...')
    
coeff_a = zeros(13,13);
coeff_b = zeros(13,13);

inpath   = '.';
filename = 'pollack_coeff.dat';

delimiterIn = ' ';
headerlinesIn = 3;
yyy = importdata([inpath '/' filename], delimiterIn, headerlinesIn);

coeff_array = yyy.data;

n = 0;

for l=0:12
   for m=0:l
      n  = n + 1;
      ll = round(coeff_array(n,1));
      mm = round(coeff_array(n,2));
      coeff_a(l+1,m+1) = coeff_array(n,3);
      coeff_b(l+1,m+1) = coeff_array(n,4);
      if (ll ~= l) || (mm ~= m)
          error ' Inconsistency in l,m values!'
      end
   end
end

%-------- Set-up of SICOPOLIS grid in the stereographic plane --------

disp(' ')
disp('Setting up SICOPOLIS grid...')

disp(' ')

disp('Horizontal resolution: (1) 40 km, (2) 20 km, (3) 10 km, (4)  5 km,')
disp('                       (5) 16 km, (6)  8 km, (7)  4 km, (8)  2 km.')
dx_number = input('Enter [1-8] > ');

if ( (dx_number < 1) || (dx_number > 8) )
   error('Wrong value!')
end

if dx_number == 1 && GRID == 2
   dx    =  40e3;
   ch_dx = '40';
elseif dx_number == 2
   dx    =  20e3;
   ch_dx = '20';
elseif dx_number == 3
   dx    =  10e3;
   ch_dx = '10';
elseif dx_number == 4
   dx    =   5e3;
   ch_dx = '05';
elseif dx_number == 5 && GRID == 2
   dx    =  16e3;
   ch_dx = '16';
elseif dx_number == 6 && GRID == 2
   dx    =   8e3;
   ch_dx = '08';
elseif dx_number == 7 && GRID == 2
   dx    =   4e3;
   ch_dx = '04';
elseif dx_number == 8 && GRID == 2
   dx    =   2e3;
   ch_dx = '02';
else
   error('Resolution not available for Bamber grid!')
end

r_imax = (xmax-xmin)/dx;
r_jmax = (ymax-ymin)/dx;

imax = round(r_imax);
jmax = round(r_jmax);

if abs(r_imax-double(imax)) > eps
   error ' No valid value of imax could be computed!'
end

if abs(r_jmax-double(jmax)) > eps
   error ' No valid value of jmax could be computed!'
end

x = xmin:dx:xmax;
y = ymin:dx:ymax;

%  ------ Corresponding latitude and longitude

lon = zeros(imax+1,jmax+1);
lat = zeros(imax+1,jmax+1);

for i=1:imax+1
for j=1:jmax+1

   [lon(i,j), lat(i,j)] ...
        = stereo_proj_m('inv_ell', {x(i), y(j), A, B, lon0, lat0});

end
end

%  ------ Read present-day topography mask

mask = zeros(imax+1,jmax+1,'int16');

if GRID == 1
    inpath   = '/uchi/greve/sicopolis/sicopolis_v33/sico_in/grl';
    filename = ['grl_b2_' ch_dx '_woem_mask.dat'];
elseif GRID == 2 && dx_number <= 4
    inpath   = '/uchi/greve/sicopolis/sicopolis/sico_in/grl';
    filename = ['grl_bm3_' ch_dx '_mask.dat'];
elseif GRID == 2 && dx_number >= 5
    inpath   = '/uchi/greve/sicopolis/sicopolis/sico_in/grl';
    filename = ['grl_bm5_' ch_dx '_topo.nc'];
end

if strcmp(filename(end-1:end), 'nc')   % NetCDF

   mask = ncread([inpath '/' filename], 'mask');

else   % ASCII

   fid = fopen([inpath '/' filename], 'r');

   for n=1:6
      fscanf(fid, '%*[^\n]\n', 1);   % skip one line
   end

   for j=jmax+1:-1:1
      for i=1:imax+1
         mask(i,j) = fscanf(fid, '%1d', 1);
      end
   end

   fclose(fid);

end

%-------- Computation of the geothermal heat flux --------

disp(' ')
disp('Computing geothermal heat flux...')

q_geo = zeros(imax+1,jmax+1);

%  ------ Global heat flow by Pollack et al. (1993)

for i=1:imax+1
for j=1:jmax+1

   for l=0:12
      for m=0:l
         q_geo(i,j) = q_geo(i,j) ...
                      + ( coeff_a(l+1,m+1)*cos(m*lon(i,j)) ...
                          + coeff_b(l+1,m+1)*sin(m*lon(i,j)) ) ...
                        * assoc_legendre(l, m, sin(lat(i,j)), true);
      end
   end

end
end

%  ------ Scaling

if strcmp(VERSION, '1.0') || strcmp(VERSION, '1.1')

%    ---- Scaling to prescribed value at GRIP

    q_geo_grip_val = 0.0;
    
    for l=0:12
        for m=0:l
            q_geo_grip_val = q_geo_grip_val ...
                + ( coeff_a(l+1,m+1)*cos(m*lon_grip) ...
                + coeff_b(l+1,m+1)*sin(m*lon_grip) ) ...
                * assoc_legendre(l, m, sin(lat_grip), true);
        end
    end
    
    q_geo = q_geo * q_geo_grip/q_geo_grip_val;

elseif strcmp(VERSION, '2.0')

%    ---- Scaling to prescribed average value

    n_mask        = 0;
    q_geo_ave_val = 0.0;

    for i=1:imax+1
    for j=1:jmax+1

        if mask(i,j) == 0 || mask(i,j) == 1
            n_mask        = n_mask + 1;
            q_geo_ave_val = q_geo_ave_val + q_geo(i,j);
        end

    end
    end

    q_geo_ave_val = q_geo_ave_val/n_mask;
    
    q_geo = q_geo * q_geo_ave/q_geo_ave_val;

end

%  ------ Interpolation between the above-computed values at the margin of
%         the Greenland domain and prescribed values at the ice-core locations

n_point_1 = 2*imax+2*jmax;   % number of margin points

if strcmp(VERSION, '1.0') || strcmp(VERSION, '1.1')
    n_point_2 = 4;   % number of ice-core locations
elseif strcmp(VERSION, '2.0')
    n_point_2 = 8;   % number of ice-core locations
end

n_point   = n_point_1 + n_point_2;

x_point     = zeros(1, n_point);
y_point     = zeros(1, n_point);
q_geo_point = zeros(1, n_point);

%    ---- Margin points

n = 0;

for i=1:imax+1
   j = 1; n = n+1;
   x_point(n) = x(i); y_point(n) = y(j); q_geo_point(n) = q_geo(i,j);
   j = jmax+1; n = n+1;
   x_point(n) = x(i); y_point(n) = y(j); q_geo_point(n) = q_geo(i,j);
end

for j=2:jmax
   i = 1; n = n+1;
   x_point(n) = x(i); y_point(n) = y(j); q_geo_point(n) = q_geo(i,j);
   i = imax+1; n = n+1;
   x_point(n) = x(i); y_point(n) = y(j); q_geo_point(n) = q_geo(i,j);
end

%    ---- Ice-core locations

n_sites = length(ch_sites);

n = n+1;
[x_point(n), y_point(n)] ...
    = stereo_proj_m('forw_ell', {lon_grip, lat_grip, A, B, lon0, lat0});
q_geo_point(n) = q_geo_grip;

n = n+1;
[x_point(n), y_point(n)] ...
    = stereo_proj_m('forw_ell', {lon_dye3, lat_dye3, A, B, lon0, lat0});
q_geo_point(n) = q_geo_dye3;

n = n+1;
[x_point(n), y_point(n)] ...
    = stereo_proj_m('forw_ell', {lon_cc, lat_cc, A, B, lon0, lat0});
q_geo_point(n) = q_geo_cc;

n = n+1;
[x_point(n), y_point(n)] ...
    = stereo_proj_m('forw_ell', {lon_ngrip, lat_ngrip, A, B, lon0, lat0});
q_geo_point(n) = q_geo_ngrip;

if strcmp(VERSION, '2.0')

    n = n+1;
    [x_point(n), y_point(n)] ...
        = stereo_proj_m('forw_ell', ...
                        {lon_neem, lat_neem, A, B, lon0, lat0});
    q_geo_point(n) = q_geo_neem;

    n = n+1;
    [x_point(n), y_point(n)] ...
        = stereo_proj_m('forw_ell', ...
                        {lon_sass1, lat_sass1, A, B, lon0, lat0});
    q_geo_point(n) = q_geo_sass1;

    n = n+1;
    [x_point(n), y_point(n)] ...
        = stereo_proj_m('forw_ell', ...
                        {lon_sass2, lat_sass2, A, B, lon0, lat0});
    q_geo_point(n) = q_geo_sass2;

    n = n+1;
    [x_point(n), y_point(n)] ...
        = stereo_proj_m('forw_ell', ...
                        {lon_gap, lat_gap, A, B, lon0, lat0});
    q_geo_point(n) = q_geo_gap;

end

if n ~= n_point
    error ' Wrong value of n_point!'
end

%    ---- Interpolation with squared-inverse-distance weighing

for i=1:imax+1
for j=1:jmax+1

   q_geo(i,j) = 0.0;
   sum_weigh  = 0.0;

   for n=1:n_point
      dist       = sqrt((x(i)-x_point(n))^2 + (y(j)-y_point(n))^2);
      weigh      = 1.0/(dist+eps)^2;
      if n <= n_point_1
         weigh = weigh/double(n_point_1);   % margin point
      else
         weigh = weigh/double(n_point_2);   % ice-core location
      end
      sum_weigh  = sum_weigh + weigh;
      q_geo(i,j) = q_geo(i,j) + weigh*q_geo_point(n);
   end

   q_geo(i,j) = q_geo(i,j)/sum_weigh;

end
end

%  ------ Adding bell-shaped anomaly at the onset of the NEGIS (optional)

if strcmp(VERSION, '2.0') ...
        && q_geo_negis_onset_anom > 0.0   % negative values ignored
        
    [x_point(1), y_point(1)] ...
        = stereo_proj_m('forw_ell', ...
                        {lon_negis_onset, lat_negis_onset, ...
                         A, B, lon0, lat0});
    sigma_negis_onset = 0.5*rad_negis_onset;
                        % 2-sigma radius -> 1-sigma radius
    
    for i=1:imax+1
    for j=1:jmax+1
        dist       = sqrt((x(i)-x_point(1))^2 + (y(j)-y_point(1))^2);
        q_geo(i,j) = q_geo(i,j) ...
                       + q_geo_negis_onset_anom ...
                           *exp(-dist^2/(2.0*sigma_negis_onset^2));
    end
    end

end

%-------- Writing of data on file --------

disp(' ')
disp('Writing data on file...')

%  ------ ASCII file

% if max(q_geo(:)) < 1000.0
%    format_q_geo = '%7.2f';
% else
%    format_q_geo = '%8.2f';
% end
%
% format_write = '';
% for i=1:imax+1; format_write = [format_write, format_q_geo]; end
% format_write = [format_write, '\n'];
%
% filename = ['GHF_' domain '_Ver' VERSION '_Grid' ch_grid ...
%             '_' ch_dx 'km'];
%
% [fid, errmess] = fopen([filename '.dat'], 'wt');
% if fid == -1;   error(errmess);   end
%
% fprintf(fid, '%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
% fprintf(fid, '%% Greenland:\n');
% fprintf(fid, '%% Geothermal heat flux, in mW/m2.\n');
% fprintf(fid, ['%% Horizontal resolution ' ch_dx ' km.\n']);
% fprintf(fid, '%% %3i records [j = %3i (-1) 0] with %3i values [i = 0 (1) %3i] in each record.\n', ...
%              jmax+1, jmax, imax+1, imax);
% fprintf(fid, '%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
%
% for j=jmax+1:-1:1
%    fprintf(fid, format_write, q_geo(:,j));
% end
%
% fclose(fid);

%  ------ NetCDF file

filename = ['GHF_' domain '_Ver' VERSION '_Grid' ch_grid ...
            '_' ch_dx 'km'];

filename_nc = [filename '.nc'];

ncid = netcdf.create(filename_nc, 'CLASSIC_MODEL');

dimid_x = netcdf.defDim(ncid, 'x', imax+1);
dimid_y = netcdf.defDim(ncid, 'y', jmax+1);

varid1 = netcdf.defVar(ncid, 'mapping', 'NC_BYTE', []);

netcdf.putAtt(ncid, varid1, 'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid, varid1, 'reference_ellipsoid_name', 'WGS84');
netcdf.putAtt(ncid, varid1, 'semi_major_axis', A);
netcdf.putAtt(ncid, varid1, 'inverse_flattening', F_INV);
netcdf.putAtt(ncid, varid1, 'latitude_of_projection_origin', 90);
netcdf.putAtt(ncid, varid1, 'standard_parallel', latd0);
netcdf.putAtt(ncid, varid1, 'straight_vertical_longitude_from_pole', lond0);
netcdf.putAtt(ncid, varid1, 'false_easting', 0);
netcdf.putAtt(ncid, varid1, 'false_northing', 0);

varid2 = netcdf.defVar(ncid, 'x', 'NC_INT', dimid_x);

netcdf.putAtt(ncid, varid2, 'units', 'm');
netcdf.putAtt(ncid, varid2, 'standard_name', 'projection_x_coordinate');
netcdf.putAtt(ncid, varid2, 'long_name', 'x-coordinate');
netcdf.putAtt(ncid, varid2, 'axis', 'x');

varid3 = netcdf.defVar(ncid, 'y', 'NC_INT', dimid_y);

netcdf.putAtt(ncid, varid3, 'units', 'm');
netcdf.putAtt(ncid, varid3, 'standard_name', 'projection_y_coordinate');
netcdf.putAtt(ncid, varid3, 'long_name', 'y-coordinate');
netcdf.putAtt(ncid, varid3, 'axis', 'y');

varid4 = netcdf.defVar(ncid, 'lat', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid4, 'units', 'degrees_N');
netcdf.putAtt(ncid, varid4, 'standard_name', 'latitude');
netcdf.putAtt(ncid, varid4, 'long_name', 'Geographical latitude');

varid5 = netcdf.defVar(ncid, 'lon', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid5, 'units', 'degrees_E');
netcdf.putAtt(ncid, varid5, 'standard_name', 'longitude');
netcdf.putAtt(ncid, varid5, 'long_name', 'Geographical longitude');

varid6 = netcdf.defVar(ncid, 'GHF', 'NC_FLOAT', [dimid_x dimid_y]);

netcdf.putAtt(ncid, varid6, 'units', 'mW m-2');
netcdf.putAtt(ncid, varid6, 'standard_name', ...
                               'upward_geothermal_heat_flux');
netcdf.putAtt(ncid, varid6, 'long_name', 'Geothermal heat flux');
netcdf.putAtt(ncid, varid6, 'grid_mapping', 'mapping');

netcdf.endDef(ncid);

mapping = 0;

netcdf.putVar(ncid, varid1, mapping);
netcdf.putVar(ncid, varid2, round(x));
netcdf.putVar(ncid, varid3, round(y));
netcdf.putVar(ncid, varid4, rad2deg(lat));
netcdf.putVar(ncid, varid5, mod(rad2deg(lon)+180, 360)-180);
netcdf.putVar(ncid, varid6, q_geo);

netcdf.close(ncid);

ncwriteatt(filename_nc, '/', 'Author'     , ...
                          'Ralf Greve <greve@lowtem.hokudai.ac.jp>');
ncwriteatt(filename_nc, '/', 'Variable'   , ...
                          'Geothermal heat flux (mW/m2) for Greenland');
ncwriteatt(filename_nc, '/', 'Version' , VERSION);
ncwriteatt(filename_nc, '/', 'Grid', ch_grid);
ncwriteatt(filename_nc, '/', 'Resolution', [ch_dx ' km']);
ncwriteatt(filename_nc, '/', 'Reference' , ch_reference);
ncwriteatt(filename_nc, '/', 'Date' , datestr(now));

%-------- Plotting --------

disp(' ')
disp('Plotting...')

x_km = x*1e-3;
y_km = y*1e-3;

x_point_km = x_point*1e-3;
y_point_km = y_point*1e-3;

figure(1)

pcolor_rg(x_km, y_km, q_geo', [25 135]);

hold on
mask_binary = 2.0.*(1.0-abs(sign(mask-0)));
[~, hh] = contour(x_km, y_km, mask_binary', [1 1]);
hh.LineColor = 'black';
hh.LineWidth = 1.5    ;
hh.LineStyle = ':'    ;

hold on
mask_binary = 2.0.*(1.0-abs(sign(mask-2)));
[~, hh] = contour(x_km, y_km, mask_binary', [1 1]);
hh.LineColor = 'black';
hh.LineWidth = 1.5    ;
hh.LineStyle = '-'    ;

NCV_banded = colormap_m('NCV_banded');
NCV_banded_bright = 1.0-NCV_banded;
NCV_banded_bright = 0.8*NCV_banded_bright;
NCV_banded_bright = 1.0-NCV_banded_bright;
colormap(NCV_banded_bright)
colorbar

axis([x_km(1) x_km(end) y_km(1) y_km(end)])

ax           = gca    ;   % axes handle
ax.Layer     = 'top'  ;
ax.FontName  = 'Arial';
ax.FontSize  = 14     ;
ax.LineWidth = 1.5    ;
ax.TickDir   = 'in'   ;
ax.XTick = -800:400:800  ;
ax.YTick = -3400:400:-600;
xlabel(ax, 'x (km)', 'FontName','Arial', 'FontSize',16)
ylabel(ax, 'y (km)', 'FontName','Arial', 'FontSize',16)
ax.PlotBoxAspectRatio = [1 (y(end)-y(1))/(x(end)-x(1)) 1];

hold on
pp = plot(x_point_km(end-n_sites+1:end), y_point_km(end-n_sites+1:end));
pp.LineStyle = 'none';
pp.Marker = 'o';
pp.MarkerSize = 8;
pp.MarkerFaceColor = 'white';
pp.MarkerEdgeColor = 'black';

hold on
pp = plot(x_point_km(end-n_sites+1:end), y_point_km(end-n_sites+1:end));
pp.LineStyle = 'none';
pp.Marker = '+';
pp.MarkerSize = 7;
pp.MarkerEdgeColor = 'black';

for n=1:n_sites

    tt = text(x_point_km(end-n_sites+n)+5, y_point_km(end-n_sites+n)+5, ...
              ch_sites{n});
    tt.Units               = 'data'  ;
    tt.FontName            = 'Arial' ;
    tt.FontSize            = 16      ;
    tt.FontWeight          = 'bold'  ;
    tt.Color               = 'white' ;
    tt.HorizontalAlignment = 'left'  ;
    tt.VerticalAlignment   = 'bottom';
end

title('Geothermal heat flux (mW/m^2)', ...
      'FontName','Arial', 'FontWeight','normal', 'FontSize',16)

% fig                   = gcf          ;   % figure handle
% fig.PaperOrientation  = 'portrait'   ;
% fig.PaperType         = 'A4'         ;
% fig.PaperUnits        = 'centimeters';
% fig.PaperPositionMode = 'manual'     ;
% fig.PaperPosition     = [0 0 20 30]  ;
% print('-dpdf', '-r300', [filename '.pdf']);
% print('-dpng', '-r600', [filename '.png']);

%-------- End of script --------

disp(' ')
disp('Done.')
disp(' ')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function assoc_legendre_r = assoc_legendre(l, m, x, norm)

%--------------------------------------------------------------------------
% Computes the associated Legendre functions P_l^m(x). Here l and m are
% integers satisfying 0 <= m <= l, while x lies in the range -1 <= x <= 1.
% Credit: Press et al., 'Numerical recipes in Fortran 77'
%         (optional normalization added by Ralf Greve)

%    norm = false: non-normalized functions
%    norm = true : normalized functions so that the square integral
%                    \int_{-1}^{1} (P_l^m(x))^2 dx is 2 (m=0) or 4 (m>0)
%--------------------------------------------------------------------------

if (m < 0) || (m > l) || (abs(x) > 1.0)
   error ' Bad argument in function assoc_legendre!'
end
 
pmm = 1.0;   % Compute P_m^m

if m > 0
   somx2 = sqrt((1.0-x)*(1.0+x));
   fact  = 1.0;
   for i=1:m
      pmm  = -pmm*fact*somx2;
      fact = fact+2.0;
   end
end

if l == m
   assoc_legendre_r = pmm;
else
   pmmp1 = x*(2*m+1)*pmm;  % Compute P_{m+1}^m
   if (l == m+1)
      assoc_legendre_r = pmmp1;
   else   % Compute P_l^m, l > m+1
        for ll=m+2:l
           pll   = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
           pmm   = pmmp1;
           pmmp1 = pll;
        end
        assoc_legendre_r = pll;
   end
end

%-------- Normalization --------

if norm

   norm_factor = (-1)^m ...
               *sqrt( double(2*l+1)*factorial(l-m)/(2.0*factorial(l+m)) );
     if m == 0
        norm_factor = norm_factor*sqrt(2.0);
     else
        norm_factor = norm_factor*2.0;
     end

     assoc_legendre_r = assoc_legendre_r*norm_factor;

end

end % function assoc_legendre

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
