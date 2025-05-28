%==========================================================================
% Module: stereo_proj_m
%
% Description:
%   Polar stereographic projection for an ellipsoidal or spherical planet.
%
% Authors: Ralf Greve, Reinhard Calov, Alex Robinson
% Date:    2024-03-13
%==========================================================================

function [ret1, ret2] = stereo_proj_m(str_function, cell_args)

if strcmp(str_function, 'forw_ell')
   [ret1, ret2] = forw_ell(cell_args);
elseif strcmp(str_function, 'forw_sph')
   [ret1, ret2] = forw_sph(cell_args);
elseif strcmp(str_function, 'inv_ell')
   [ret1, ret2] = inv_ell(cell_args);
elseif strcmp(str_function, 'inv_sph')
   [ret1, ret2] = inv_sph(cell_args);
end

end % module stereo_proj_m

%--------------------------------------------------------------------------
% Forward polar stereographic projection for an ellipsoidal planet.
%--------------------------------------------------------------------------
function [x_val, y_val] = forw_ell(cell_args)

[lon_val, lat_val, A, B, lon0, lat0] = cell_args{:};

% Input variables:
%   lon_val - geographic longitude [rad]
%   lat_val - geographic latitude [rad]
%   A - semi-major axis of planet [m]
%   B - semi-minor axis of planet [m]
%       [ B = A*(1-F), where F is the flattening ]
%   lon0 - central meridian of projection [rad]
%   lat0 - standard parallel of projection [rad]

% Output variables:
%   x_val - projection x coordinate [m]
%   y_val - projection y coordinate [m]

% Call:
%   [x_val, y_val] ...
%      = stereo_proj_m('forw_ell', {lon_val, lat_val, A, B, lon0, lat0});

%-------- Parameters, constraints --------

eps = 1.0e-05;

latd0 = rad2deg(lat0);

if lat0 > eps   % for northern hemisphere
   sign_lat0 =  1;
   latd0 = min(latd0,  (90.0-eps));
elseif lat0 < (-eps)   % for southern hemisphere
   sign_lat0 = -1;
   latd0 = max(latd0, -(90.0-eps));
else
   disp(' >>> forw_ell: lat0 must be different from zero!')
   return
end

lat0_aux = deg2rad(latd0) * sign_lat0;

e = sqrt((A^2-B^2)/(A^2));

%-------- Stereographic coordinates x,y --------

sinlat0 = sin(lat0_aux);
coslat0 = cos(lat0_aux);

lat_aux = lat_val * sign_lat0;
  
mc = coslat0/sqrt(1.0-e*e*sinlat0*sinlat0);
t = sqrt(((1.0-sin(lat_aux))/(1.0+sin(lat_aux)))* ...
        ((1.0+e*sin(lat_aux))/(1.0-e*sin(lat_aux)))^e);
tc = sqrt(((1.0-sinlat0)/(1.0+sinlat0))* ...
         ((1.0+e*sinlat0)/(1.0-e*sinlat0))^e);
rho = A*mc*t/tc;

x_val =              rho*sin(lon_val-lon0);
y_val = -sign_lat0 * rho*cos(lon_val-lon0);

%-------- End of function --------

end % function forw_ell

%--------------------------------------------------------------------------
% Forward polar stereographic projection for a spherical planet.
%--------------------------------------------------------------------------
function [x_val, y_val] = forw_sph(cell_args)

[lon_val, lat_val, R, lon0, lat0] = cell_args{:};

% Input variables:
%   lon_val - geographic longitude [rad]
%   lat_val - geographic latitude [rad]
%   R - radius of planet [m]
%   lon0 - central meridian of projection [rad]
%   lat0 - standard parallel of projection [rad]

% Output variables:
%   x_val - projection x coordinate [m]
%   y_val - projection y coordinate [m]

% Call:
%   [x_val, y_val] ...
%      = stereo_proj_m('forw_sph', {lon_val, lat_val, R, lon0, lat0});

%-------- Parameters, constraints --------

eps = 1.0e-05;

latd0 = rad2deg(lat0);

if lat0 > eps   % for northern hemisphere
   sign_lat0 =  1;
   latd0 = min(latd0,  90.0);
elseif lat0 < (-eps)   % for southern hemisphere
   sign_lat0 = -1;
   latd0 = max(latd0, -90.0);
else
   disp(' >>> forw_sph: lat0 must be different from zero!')
   return
end

lat0_aux = deg2rad(latd0) * sign_lat0;

K = (cos(0.25*pi-0.5*lat0_aux))^2;

%-------- Stereographic coordinates x,y --------

lat_aux = lat_val * sign_lat0;
  
x_val =              2.0*R*K*tan(0.25*pi-0.5*lat_aux) ...
                            *sin(lon_val-lon0);

y_val = -sign_lat0 * 2.0*R*K*tan(0.25*pi-0.5*lat_aux) ...
                            *cos(lon_val-lon0);

%-------- End of function --------

end % function forw_sph

%--------------------------------------------------------------------------
% Inverse polar stereographic projection for an ellipsoidal planet.
%--------------------------------------------------------------------------
function [lon_val, lat_val] = inv_ell(cell_args)

[x_val, y_val, A, B, lon0, lat0] = cell_args{:};

% Input variables:
%   x_val - projection x coordinate [m]
%   y_val - projection y coordinate [m]
%   A - semi-major axis of planet [m]
%   B - semi-minor axis of planet [m]
%       [ B = A*(1-F), where F is the flattening ]
%   lon0 - central meridian of projection [rad]
%   lat0 - standard parallel of projection [rad]

% Output variables:
%   lon_val - geographic longitude [rad]
%   lat_val - geographic latitude [rad]

% Call:
%   [lon_val, lat_val] ...
%      = stereo_proj_m('inv_ell', {x_val, y_val, A, B, lon0, lat0});

%-------- Parameters, constraints --------

eps = 1.0e-05;
eps_residual = 1.0e-09;

latd0 = rad2deg(lat0);

if lat0 > eps   % for northern hemisphere
   sign_lat0 =  1;
   latd0 = min(latd0,  (90.0-eps));
elseif lat0 < (-eps)   % for southern hemisphere
   sign_lat0 = -1;
   latd0 = max(latd0, -(90.0-eps));
else
   disp(' >>> inv_ell: lat0 must be different from zero!')
   return
end

lat0_aux = deg2rad(latd0) * sign_lat0;

e = sqrt((A^2-B^2)/(A^2));

%-------- Longitude --------

if (x_val ~= 0.0) || (y_val ~= 0.0)
   lon_val = lon0 + sign_lat0*atan2(y_val,x_val) + 0.5*pi;
else
   lon_val = lon0 + 0.5*pi;
end

%-------- Fix-point iteration for latitude --------

sinlat0 = sin(lat0_aux);
coslat0 = cos(lat0_aux);

tc = sqrt(((1.0-sinlat0)/(1.0+sinlat0))* ...
         ((1.0+e*sinlat0)/(1.0-e*sinlat0))^e);
mc = coslat0/sqrt(1.0-e*e*sinlat0*sinlat0);
rho = sqrt(x_val*x_val+y_val*y_val);
t = rho*tc/(A*mc);

lat_p = 0.5*pi-2.0*atan(t);
l = 0;
residual = 3600.0;

while (residual >= eps_residual)
    l = l+1;
    lat_aux = 0.5*pi-2.0*atan(t*((1.0-e*sin(lat_p))/ ...
              (1.0+e*sin(lat_p)))^(0.5*e));
    residual = abs(lat_aux-lat_p);
    lat_p = lat_aux;
end

lat_val = lat_aux * sign_lat0;

%-------- Constrain longitude to [0, 2*pi) --------

if lon_val < 0.0
   lon_val = lon_val + 2.0*pi;
elseif lon_val >= (2.0*pi)
   lon_val = lon_val - 2.0*pi;
end

%-------- End of function --------

end % function inv_ell

%--------------------------------------------------------------------------
% Inverse polar stereographic projection for a spherical planet.
%--------------------------------------------------------------------------
function [lon_val, lat_val] = inv_sph(cell_args)

[x_val, y_val, R, lon0, lat0] = cell_args{:};

% Input variables:
%   x_val - projection x coordinate [m]
%   y_val - projection y coordinate [m]
%   R - radius of planet [m]
%   lon0 - central meridian of projection [rad]
%   lat0 - standard parallel of projection [rad]

% Output variables:
%   lon_val - geographic longitude [rad]
%   lat_val - geographic latitude [rad]

% Call:
%   [lon_val, lat_val] ...
%      = stereo_proj_m('inv_sph', {x_val, y_val, R, lon0, lat0});

%-------- Parameters, constraints --------

eps = 1.0e-05;

latd0 = rad2deg(lat0);

if lat0 > eps   % for northern hemisphere
   sign_lat0 =  1;
   latd0 = min(latd0,  90.0);
elseif lat0 < (-eps)   % for southern hemisphere
   sign_lat0 = -1;
   latd0 = max(latd0, -90.0);
else
   disp(' >>> inv_sph: lat0 must be different from zero!')
   return
end

lat0_aux = deg2rad(latd0) * sign_lat0;

K = (cos(0.25*pi-0.5*lat0_aux))^2;

%-------- Longitude --------

if (x_val ~= 0.0) || (y_val ~= 0.0)
   lon_val = lon0 + sign_lat0*atan2(y_val,x_val) + 0.5*pi;
else
   lon_val = lon0 + 0.5*pi;
end

%-------- Latitude --------

lat_aux = 0.5*pi - 2.0*atan(sqrt(x_val^2+y_val^2)/(2.0*R*K));

lat_val = lat_aux * sign_lat0;

%-------- Constrain longitude to [0, 2*pi) --------

if lon_val < 0.0
   lon_val = lon_val + 2.0*pi;
elseif lon_val >= (2.0*pi)
   lon_val = lon_val - 2.0*pi;
end

%-------- End of function --------

end % function inv_sph
%
