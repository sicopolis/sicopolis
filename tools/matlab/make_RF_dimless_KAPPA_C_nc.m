%==========================================================================
% make_RF_dimless_KAPPA_C_nc
%
% Description:
%   Creating a NetCDF file (RF_KAPPA_C file) for the tabulated values
%   of the temperature-dependent rate factor (dimensionless),
%   heat conductivity and specific heat required by SICOPOLIS.
%
% Author: Ralf Greve
% Date:   2025-06-18
%==========================================================================

clear variables
close all

%-------- Parameters, settings --------

disp(' ')
disp('Rate factor:')
disp(' (1) Constant dimensionless value 1')
disp(' (2) Constant value for Vialov profile')
disp(' (3) Constant value for CalvingMIP')
disp(' (4) EISMINT Phase 2 SGE')
disp(' (5) Greve, Weis and Hutter (1998)')
disp(' (6) Cuffey and Paterson (2010)')
disp(' (7) Smith and Morland (1981)')
disp(' ')

prompt = 'Enter [1-7] ';
n_RF   = input(prompt);

if isempty(n_RF) || ( n_RF < 1 || n_RF > 7 )
   error('Wrong value of n_RF!')
end

disp(' ')
disp('Heat conductivity:')
disp(' (1) Constant value (EISMINT Phase 2 SGE)')
disp(' (2) Greve, Weis and Hutter (1998)')
disp(' ')

prompt  = 'Enter [1-2] ';
n_KAPPA = input(prompt);

if isempty(n_KAPPA) || ( n_KAPPA < 1 || n_KAPPA > 2 )
   error('Wrong value of n_KAPPA!')
end

disp(' ')
disp('Specific heat:')
disp(' (1) Constant value (EISMINT Phase 2 SGE)')
disp(' (2) Greve, Weis and Hutter (1998)')
disp(' ')

prompt = 'Enter [1-2] ';
n_C    = input(prompt);

if isempty(n_C) || ( n_C < 1 || n_C > 2 )
   error('Wrong value of n_C!')
end

temp_C2K     = 273.15;
gas_constant = 8.314;
year2sec     = 31556925.445;

chtime = datetime;
chtime.Format = 'yyyyMMdd''T''HHmmss';
filename = ['my_RF_dimless_KAPPA_C_file_' ...
            int2str(n_RF) '_' int2str(n_KAPPA) '_' int2str(n_C) '_' ...
            char(chtime) '.nc'];

%-------- Temperature --------

tempC = zeros(1,201, 'int32');
tempK = zeros(1,201);

for n=1:201
   tempC(n) = n - 191;   % -190 degC ... +10 degC
   tempK(n) = temp_C2K + double(tempC(n));
end

%-------- Rate factor RF(T)
%         [RF] = 1 (dimensionless), [T] = degC --------

RF = zeros(1,201);

if n_RF == 1
%  ------ RF(T) = const = 1 (dimensionless)

   RF(:) = 1.0;

   ch_ref_RF = 'Constant value RF = 1 (dimensionless)';

elseif n_RF == 2
%  ------ RF(T) = const, Vialov

   RF(:) = 1.0e-16/year2sec;

   ch_ref_RF = 'Constant value for Vialov profile';

elseif n_RF == 3
%  ------ RF(T) = const, CalvingMIP

   RF(:) = 9.30920838e-26;

   ch_ref_RF = ['CalvingMIP project, ' ...
                'https://github.com/JRowanJordan/CalvingMIP/wiki'];

elseif n_RF == 4
%  ------ RF(T), set-up of EISMINT Phase 2 SGE

   for n=1:180
      RF(n) = 3.61e-13 * exp( -6.0e+04/(gas_constant*tempK(n)));
   end

   for n=181:191
      RF(n) = 1.73e+03 * exp(-13.9e+04/(gas_constant*tempK(n)));
   end

   for n=192:201
      RF(n) = RF(191);   % Dummies for temperatures above 0 degC
   end

   ch_ref_RF = ['EISMINT Phase 2 Simplified Geometry Experiments; ' ...
                'Payne et al., 2000, J. Glaciol. 46, 227-238'];

elseif n_RF == 5
%  ------ RF(T) by Greve, Weis and Hutter, 1998,
%                     Paleoclimates 2 (2-3), 133-161

   for n=1:180
      RF(n) = 3.985e-13 * exp(-6.0e+04/(gas_constant*tempK(n)));
   end

   RF(181) = 4.9e-25;   % for -10 degC
   RF(186) = 1.6e-24;   % for  -5 degC
   RF(189) = 2.4e-24;   % for  -2 degC
   RF(191) = 3.2e-24;   % for   0 degC

   RF(182) = 0.8*RF(181) + 0.2*RF(186);
   RF(183) = 0.6*RF(181) + 0.4*RF(186);
   RF(184) = 0.4*RF(181) + 0.6*RF(186);
   RF(185) = 0.2*RF(181) + 0.8*RF(186);

   RF(187) = (2.0/3.0)*RF(186) + (1.0/3.0)*RF(189);
   RF(188) = (1.0/3.0)*RF(186) + (2.0/3.0)*RF(189);

   RF(190) = 0.5*RF(189) + 0.5*RF(191);

   for n=192:201
      RF(n) = RF(191);   % Dummies for temperatures above 0 degC
   end

   ch_ref_RF = 'Greve, Weis & Hutter, 1998, Paleoclimates 2, 133-161';

elseif n_RF == 6
%  ------ RF(T) by Cuffey and Paterson, "The Physics of Glaciers",
%                                      4th ed. 2010, Sect. 3.4.6

   for n=1:180
      RF(n) = 3.5e-25 ...
                 * exp( ( -6.0e+04/gas_constant) ...
                        * (1.0/tempK(n)-1.0/(temp_C2K-10.0)) );
   end

   for n=181:191
      RF(n) = 3.5e-25 ...
                 * exp( (-11.5e+04/gas_constant) ...
                        * (1.0/tempK(n)-1.0/(temp_C2K-10.0)) );
   end

   for n=192:201
      RF(n) = RF(191);   % Dummies for temperatures above 0 degC
   end

   ch_ref_RF = ['Cuffey and Paterson, The Physics of Glaciers, ' ...
                '4th ed. 2010, Sect. 3.4.6'];

elseif n_RF == 7
%  ------ RF(T) by Smith and Morland, 1981,
%                  Cold Reg. Sci. Technol. 5, 141-150

   c1 =  0.7242;
   c2 =  0.3438;
   d1 = 11.9567;
   d2 =  2.9494;
   DT = 20.0;

   for n=1:191
      T_bar = (tempK(n)-temp_C2K)/DT;
      RF(n) = c1*exp(d1*T_bar) + c2*exp(d2*T_bar);
   end

   for n=192:201
      RF(n) = RF(191);   % Dummies for temperatures above 0 degC
   end

   T10_bar = -10.0/DT;
   RF10    = c1*exp(d1*T10_bar) + c2*exp(d2*T10_bar);
   RF      = RF/RF10;   % Normalization to -10 degC

   ch_ref_RF = ['Smith and Morland, 1981, ' ...
                'Cold Reg. Sci. Technol. 5, 141-150; ' ...
                'normalized to -10 degC'];

end

%  ------ Scaling, non-dimensionalization

if n_RF < 7
   stress_dev_scale       = 1.0e+05;
   stress_dev_scale_unit  = 'Pa';
   strain_rate_scale      = 2.5e-02;
   strain_rate_scale_unit = 'a-1';
elseif n_RF == 7
   stress_dev_scale       = 1.0;
   stress_dev_scale_unit  = '-';
   strain_rate_scale      = 1.0;
   strain_rate_scale_unit = '-';
end

if n_RF < 7
   RF_scale = (strain_rate_scale/year2sec)/stress_dev_scale^3;
elseif n_RF == 7
   RF_scale = 1.0;
   %%% RF_scale = (strain_rate_scale/year2sec)/stress_dev_scale;
   %%%       Scaling needs to be checked!
end

if n_RF > 1 && n_RF < 7
RF = RF/RF_scale;
end

%-------- Heat conductivity KAPPA(T)
%         [KAPPA] = W/(m*K), [T] = degC --------

KAPPA = zeros(1,201);

if n_KAPPA == 1
%  ------ KAPPA(T) = const, by EISMINT Phase 2 SGE

   KAPPA(:) = 2.1;

   ch_ref_KAPPA = ['EISMINT Phase 2 Simplified Geometry Experiments; ' ...
                   'Payne et al., 2000, J. Glaciol. 46, 227-238'];

elseif n_KAPPA == 2
%  ------ KAPPA(T) by Greve, Weis and Hutter, 1998,
%                     Paleoclimates 2 (2-3), 133-161

   for n=1:191
      KAPPA(n) = 9.828*exp(-0.0057*tempK(n));
   end

   for n=192:201
      KAPPA(n) = KAPPA(191);   % Dummies for temperatures above 0 degC
   end

   ch_ref_KAPPA = 'Greve, Weis & Hutter, 1998, Paleoclimates 2, 133-161';

end

%-------- Specific heat C(T)
%         [C] = J/(kg*K), [T] = degC --------

C = zeros(1,201);

if n_C == 1
%  ------ C(T) = const, by EISMINT Phase 2 SGE

   C(:) = 2009.0;

   ch_ref_C = ['EISMINT Phase 2 Simplified Geometry Experiments; ' ...
               'Payne et al., 2000, J. Glaciol. 46, 227-238'];

elseif n_C == 2
%  ------ C(T) by Greve, Weis and Hutter, 1998,
%                 Paleoclimates 2 (2-3), 133-161

   for n=1:191
      C(n) = 146.3 + 7.253*tempK(n);
   end

   for n=192:201
      C(n) = C(191);   % Dummies for temperatures above 0 degC
   end

   ch_ref_C = 'Greve, Weis & Hutter, 1998, Paleoclimates 2, 133-161';

end

%-------- Write values on file --------

disp(' ')
disp(['Creating ' filename ' ...'])

nccreate(filename, ...
         'tempC', ...
         'Dimensions', {'tempC', 201}, ...
         'Datatype', 'int32', ...
         'Format', 'classic')

ncwrite(filename, 'tempC', tempC)
ncwriteatt(filename, 'tempC', 'standard_name', 'temperature')
ncwriteatt(filename, 'tempC', 'long_name', 'Temperature')
ncwriteatt(filename, 'tempC', 'units', 'degC')

nccreate(filename, ...
         'RF_dimless', ...
         'Dimensions', {'tempC', 201}, ...
         'Datatype', 'double', ...
         'Format', 'classic')

ncwrite(filename, 'RF_dimless', RF)
ncwriteatt(filename, 'RF_dimless', 'standard_name', 'ice_rate_factor')
ncwriteatt(filename, 'RF_dimless', 'long_name', 'Rate factor of ice')
ncwriteatt(filename, 'RF_dimless', 'units', '-')
ncwriteatt(filename, 'RF_dimless', 'stress_dev_scale', ...
                                    stress_dev_scale)
ncwriteatt(filename, 'RF_dimless', 'stress_dev_scale_unit', ...
                                    stress_dev_scale_unit)
ncwriteatt(filename, 'RF_dimless', 'strain_rate_scale', ...
                                    strain_rate_scale)
ncwriteatt(filename, 'RF_dimless', 'strain_rate_scale_unit', ...
                                    strain_rate_scale_unit)
ncwriteatt(filename, 'RF_dimless', 'year2sec', year2sec)
ncwriteatt(filename, 'RF_dimless', 'reference', ch_ref_RF)

nccreate(filename, ...
         'KAPPA', ...
         'Dimensions', {'tempC', 201}, ...
         'Datatype', 'double', ...
         'Format', 'classic')

ncwrite(filename, 'KAPPA', KAPPA)
ncwriteatt(filename, 'KAPPA', 'standard_name', 'ice_heat_conductivity')
ncwriteatt(filename, 'KAPPA', 'long_name', 'Heat conductivity of ice')
ncwriteatt(filename, 'KAPPA', 'units', 'W m-1 K-1')
ncwriteatt(filename, 'KAPPA', 'reference', ch_ref_KAPPA)

nccreate(filename, ...
         'C', ...
         'Dimensions', {'tempC', 201}, ...
         'Datatype', 'double', ...
         'Format', 'classic')

ncwrite(filename, 'C', C)
ncwriteatt(filename, 'C', 'standard_name', 'ice_specific_heat')
ncwriteatt(filename, 'C', 'long_name', 'Specific heat of ice')
ncwriteatt(filename, 'C', 'units', 'J kg-1 K-1')
ncwriteatt(filename, 'C', 'reference', ch_ref_C)

%-------- End of script --------

disp(' ')
disp('Done.')
disp(' ')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
