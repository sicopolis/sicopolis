.. _getting_started:

Getting started
***************

Requirements
============

Unix-like system (e.g., Linux).

| Fortran compiler.
| So far, the GNU GCC (gfortran) and Intel Fortran (ifort) compilers are supported. If you wish to use a different compiler, please contact help@sicopolis.net.

SICOPOLIS writes output in `NetCDF <https://doi.org/10.5065/D6H70CW6>`__ format (plus some ASCII). An installation of NetCDF version 3.6.x or newer is therefore required. For installation tips, see :ref:`Dependencies/NetCDF <dependencies-netcdf>`.

For the shallow-shelf/shelfy-stream solver, the Library of Iterative Solvers for Linear Systems (`Lis <https://www.ssisc.org/lis/>`__, version 1.4.43 or newer) is required. For installation tips, see :ref:`Dependencies/Lis <dependencies-lis>`.

Download
========

Clone the Git repository (recommended)
  The Git repository of SICOPOLIS is kindly hosted by the GitLab instance of the `Alfred Wegener Institute for Polar and Marine Research (AWI) <https://www.awi.de/>`__ in Bremerhaven, Germany. Front page: https://gitlab.awi.de/sicopolis/sicopolis/.

  Cloning the latest revision on the main branch::

    git clone https://gitlab.awi.de/sicopolis/sicopolis.git

  (Cloning with SSH instead of HTTPS is also available. See the above GitLab front page link for details.)

  You should then have a new directory ``sicopolis`` that contains the entire program package.

Download an archive
  Go to https://gitlab.awi.de/sicopolis/sicopolis/-/tree/main, click on the download symbol and choose the desired format (zip, tar.gz, tar.bz2, tar).

  Tagged versions of SICOPOLIS are also available from `Zenodo <https://doi.org/10.5281/zenodo.3687337>`__.

Initial configuration
=====================

1. Change to the new directory ``sicopolis`` and execute the bash script ``copy_templates.sh``::

      ./copy_templates.sh

   It copies several scripts from ``runs/templates`` to ``runs`` and the run-specs header files from ``runs/headers/templates`` to ``runs/headers``. This allows modifying the scripts and headers suitably if needed, while the original files are always stored in the respective templates subdirectories for reference. 

2. Execute the bash script ``get_input_files.sh``::

      ./get_input_files.sh

  It downloads the input data files for the several model domains (Antarctica, Greenland, etc.) These files are stored on a server (`Zenodo archive <https://doi.org/10.5281/zenodo.6371122>`__) and needed for various inputs such as topography, precipitation, geothermal heat flux, etc. The script can be configured before execution if the input files are only needed for selected domains (default is downloading everything). To do so, open it with a text editor and change the flag variables according to the instructions in the script.

3. Locate the file ``sico_configs.sh`` in the directory ``runs``, and open it with a text editor.

   Set the flags ``LIS_FLAG``, ``OPENMP_FLAG`` and ``LARGE_DATA_FLAG`` according to your needs. 

   Default is ``"true"``/``"true"``/``"false"``, which works for all test simulations included in the SICOPOLIS package. ``LIS_FLAG`` and ``OPENMP_FLAG`` can be set to ``"false"`` for simulations with pure shallow-ice dynamics. However, ``"true"`` is required for simulations with shallow-shelf dynamics (for floating ice) or hybrid shallow-ice--shelfy-stream dynamics (for grounded ice). For high-resolution simulations (e.g., Greenland/5 km or Antarctica/8 km), ``LARGE_DATA_FLAG`` must be set to ``"true"``.

   Set ``NETCDFHOME`` to the correct path of your NetCDF installation.

   If ``LIS_FLAG="true"``, set ``LISHOME`` to the correct path of your Lis installation.

   Depending on your system, some additional settings might have to be added in ``sico_configs.sh`` (``module load`` commands for dynamic loading etc.).

4. Locate the file ``sico_environment.sh`` in the directory ``runs``, open it with a text editor, and replace the ``SICO_INSTITUTION="Default"`` entry by the name of your institution (max. 256 characters).

Directory structure
===================

Main directory
  Initialization scripts ``copy_templates.sh``, ``get_input_files.sh``.

  License file.

Directory ``runs``
  Configuration scripts ``sico_configs.sh``, ``sico_environment.sh``.

  Shell script (bash) ``sico.sh`` for running a single simulation.

  Shell scripts (bash) ``multi_sico_1.sh`` and ``multi_sico_2.sh`` for running multiple simulations by repeated calls of ``sico.sh``.

  Subdirectory ``headers``
    Run-specs header files ``sico_specs_<run_name>.h`` (see :ref:`below <getting_started-run_specs_headers>`) for the simulations to be carried out with SICOPOLIS.

    By default, it contains a number of :ref:`test simulations <test_simulations>`.

Directory ``src``
  Main program file ``sicopolis.F90``.

  Subdirectory ``subroutines/general``: general subroutines, for any modelled domain.
  
  Subdirectory ``subroutines/ant``: subroutines specific for the Antarctic ice sheet.

  Subdirectory ``subroutines/grl``: subroutines specific for the Greenland ice sheet.

  Subdirectory ``subroutines/eismint``: subroutines specific for the EISMINT simplified geometry experiments.

  Accordingly subdirectories ``subroutines/asf``, ``nhem``, ``scand``, ``tibet``, ``nmars`` and ``smars`` for Austfonna, the northern hemisphere, Scandinavia, Tibet and the north and south polar caps of Mars, respectively.

  Subdirectory ``subroutines/tapenade``: AD-specific subroutines and files.

  Subdirectory ``subroutines/xyz``: For :ref:`creating new domains <new_domain>`.

Directory ``sico_in``
  Input data files for SICOPOLIS.

  Subdirectory ``general``: general input files, for any modelled domain.

  Subdirectory ``ant``: input files specific for the Antarctic ice sheet. 

  Subdirectory ``grl``: input files specific for the Greenland ice sheet.

  Subdirectory ``eismint``: input files specific for the EISMINT simplified geometry experiments.

  Accordingly subdirectories ``asf``, ``nhem``, ``scand``, ``tibet``, ``nmars`` and ``smars`` for Austfonna, the northern hemisphere, Scandinavia, Tibet and the north and south polar caps of Mars, respectively.

  Subdirectory ``xyz``: For :ref:`creating new domains <new_domain>`.

  NOTE: These subdirectories also contain README files that describe the input data and provide the corresponding references.

Directory ``sico_out``
  Directory into which output files of SICOPOLIS simulations are written by default.

Directory ``docs``
  Documentation for SICOPOLIS.

Directory ``tools``
  Some useful tools and a shell script (``tools.sh``) to execute them (see ":ref:`plotting_and_tools`").

Directory ``test_ad``
  AD-specific utilities and CI testing framework.

.. _getting_started-run_specs_headers:

Run-specs header files
======================

Each simulation (run) must be specified by a run-specs header file (or "header" for short). If the name of the simulation is supposed to be ``<run_name>``, then the name of the header must be ``sico_specs_<run_name>.h``. SICOPOLIS actually extracts the name of the simulation from the name of the header according to this pattern.

A header consists of a pretty large number of preprocessor directives of the form

.. code-block:: fortran

  #define PARAMETER value

These allow specifying many aspects of a simulation and are documented in the headers themselves. See also the ":ref:`modelling_choices`" section.

For a number of :ref:`test simulations <test_simulations>`, the run-specs header files are contained in the SICOPOLIS repository. Further examples can be found in the several paper-accompanying datasets on `Zenodo <https://zenodo.org/communities/sicopolis/>`__.

.. _getting_started-phys_para:

Physical-parameter files
========================

In these files, a number of physical parameters (densities, acceleration due to gravity, heat conductivity, specific heat, latent heat, etc.) are defined. SICOPOLIS expects them in the respective directory for the input files (``sico_in/ant`` for Antarctica, ``sico_in/grl`` for Greenland, etc.) If the name of the file is ``phys_para_xxx.dat``, it must be specified in the run-specs header file as

.. code-block:: fortran

  #define PHYS_PARA_FILE 'phys_para_xxx.dat'

The physical-parameter files can be provided in either ASCII or NetCDF format. The file type is recognized automatically by the extension (``*.nc`` for NetCDF, otherwise ASCII is assumed).

How to run a simulation
=======================

For example, to run the EISMINT Phase 2 Simplified Geometry Experiment A (Payne et al. :cite:`payne_etal_2000`), named ``repo_emtp2sge25_expA``, use the script ``sico.sh`` as follows::

  (./sico.sh -m repo_emtp2sge25_expA) >out_001.dat 2>&1

(from directory ``runs``, bash required). Accordingly for any other simulations.

For further options, try ``./sico.sh -h``.

WARNING: Do not use out\_\<run\_name\>.dat for the redirected output of ``sico.sh``. This name is reserved for the runtime output of SICOPOLIS itself. (Both are very useful in case of compilation or runtime errors!)

Alternatively, if you prefer to run :ref:`all EISMINT, Antarctica and Greenland simulations <test_simulations>` consecutively, execute the script ``multi_sico_1.sh``::

  (./multi_sico_1.sh) >out_multi_100.dat 2>&1 &

For further options, try ``./multi_sico_1.sh -h``.

Approximate computing times are listed in the ":ref:`Test simulations <test_simulations>`" section.

.. _getting_started-output:

Output files
============

Output files are written by default to the directory ``sico_out/<run_name>`` (this can be changed with the ``-d /path/to/output/directory`` option). Four types are produced:

``<run_name>.log``:
  ASCII file that lists the main specifications of simulation ``<run_name>``.

``<run_name>.ser``, ``<run_name>_ser.nc``:
  Time-series files (ASCII, NetCDF) that contain scalar variables:

  * Time, t
  * Surface temperature anomaly, D\_Ts, or glacial index, glac\_ind (forcing)
  * Sea level, z\_sl (forcing)
  * Total ice volume, V
  * Volume of grounded ice, V\_g
  * Volume of floating ice, V\_f
  * Total ice area, A
  * Area of grounded ice, A\_g
  * Area of floating ice, A\_f
  * Ice volume above flotation in sea level equivalent, V\_sle
  * Volume of temperate ice, V\_t
  * Area of temperate-based grounded ice, A\_t
  * Maximum ice thickness, H\_max
  * Maximum thickness of temperate ice, H\_t\_max
  * Maximum surface elevation, zs\_max
  * Maximum surface speed, vs\_max
  * Maximum basal temperature (relative to pmp), Tbh\_max
  * (Some more in the NetCDF file, try ``ncdump -h <run_name>_ser.nc``)

``<run_name>.core``, ``<run_name>_core.nc``:
  Time-series files (ASCII, NetCDF) that contain for selected locations xxx:

  * Time, t
  * Surface temperature anomaly, D\_Ts, or glacial index, glac\_ind (forcing)
  * Sea level, z\_sl (forcing)
  * Thickness, H\_xxx
  * Surface velocity, v\_xxx
  * Basal temperature, T\_xxx
  * (Some more in the NetCDF file, try ``ncdump -h <run_name>_core.nc``)

  | For the Greenland ice sheet, these data are written for seven locations:
  | GRIP (xxx=GR), GISP2 (xxx=G2), Dye 3 (xxx=D3), Camp Century (xxx=CC), NorthGRIP (xxx=NG), NEEM (xxx=NE), EastGRIP (xxx=EG).

  | For the Antarctic ice sheet, these data are written for six locations:
  | Vostok (xxx=Vo), Dome A (xxx=DA), Dome C (xxx=DC), Dome F (xxx=DF), Kohnen (xxx=Ko), Byrd (xxx=By).

``<run_name>0001.nc``, ``<run_name>0002.nc``, ...:
  Complete set of fields (topography, velocity, temperature etc., written in NetCDF format) for selected time slices defined in the run-specs header file.

  For example, simulation ``repo_emtp2sge25_expA`` produces three files ``repo_emtp2sge25_expA0001.nc``, ``repo_emtp2sge25_expA0002.nc`` and ``repo_emtp2sge25_expA0003.nc``, which correspond to the times :math:`t=5\,\mathrm{ka}`, :math:`50\,\mathrm{ka}` and :math:`200\,\mathrm{ka}`, respectively.
