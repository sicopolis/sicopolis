.. _plotting_and_tools:

Plotting and tools
******************

.. _plotting_and_tools-plotting:

Plotting
========

The output described in ":ref:`getting_started-output`" can be visualized with any plotting tool at the user's preference. Ncview (http://meteora.ucsd.edu/~pierce/ncview_home_page.html) is a very nice browser for NetCDF files to get a quick and easy look. For more sophisticated plots, one possibility is to use MATLAB, which has an `extensive library for NetCDF files <https://www.mathworks.com/help/matlab/network-common-data-form.html>`__. For instance, the following script plots the final surface topography of the Greenland simulation ``repo_grl16_bm5_ss25ka`` (credit: Mathieu Morlighem, University of California Irvine). ::

  filename = 'repo_grl16_bm5_ss25ka0003.nc';
  x = ncread(filename,'x');
  y = ncread(filename,'y');
  surf = ncread(filename,'zs');
  % Display surface elevation
  %    (transposition needed because MATLAB is column-oriented)
  imagesc(x*1e-3,y*1e-3,surf'); axis xy equal; caxis([0 3200]); colorbar

.. _plotting_and_tools-ismip_output:

Make ISMIP output
=================

Fortran program (located in sicopolis/tools/make_ismip_output) that generates `ISMIP6 output <https://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections2300-Antarctica#A2.3.3_Table_A1:_Variable_request_for_ISMIP6>`__ from the NetCDF time-slice files produced by SICOPOLIS (see ":ref:`getting_started-output`"). For simulation ``<run_name>``, to be executed by ::

  ./tools.sh -p make_ismip_output -m <run_name>

To list further options, execute ``./tools.sh -h``.

.. _plotting_and_tools-res_dbl:

Resolution doubler
==================

Fortran program (located in sicopolis/tools/resolution_doubler) that doubles the horizontal resolution of a NetCDF time-slice output file produced by SICOPOLIS (see ":ref:`getting_started-output`"). For simulation ``<run_name>``, to be executed by ::

  ./tools.sh -p resolution_doubler -m <run_name>

To list further options, execute ``./tools.sh -h``.

For example, run ``repo_grl10_b2_paleo21`` (10 km resolution) requires the resolution-doubled output of run ``repo_grl20_b2_paleo21`` (20 km resolution) for :math:`t=-9\,\mathrm{ka}` as initial condition. In order to create it, execute the resolution doubler for run ``repo_grl20_b2_paleo21`` (i.e., with the option ``-m repo_grl20_b2_paleo21``) and enter ::

  Number of time-slice file (with leading zeros, 4 digits) >
  0004

This will convert the original time-slice file ``repo_grl20_b2_paleo210004.nc`` to the resolution-doubled file ``repo_grl20_b2_paleo21_dbl_0004.nc`` that serves as initial conditions for run ``repo_grl10_b2_paleo21``.
