.. _getting_started:

Getting started
***************

Requirements
============

Unix-like system (e.g., Linux).

| Fortran compiler.
| So far, the GNU GCC (gfortran) and Intel Fortran (ifort) compilers are supported. If you wish to use a different compiler, please contact info@sicopolis.net.

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

1. Change to the main directory ``sicopolis`` and execute the bash script ``copy_templates.sh``::

      ./copy_templates.sh

   It copies several scripts from ``templates`` to ``.`` (the main directory), from ``tools/templates`` to ``tools``, and the run-specs header files from ``headers/templates`` to ``headers``. This allows modifying the scripts and headers suitably if needed, while the original files are always stored in the respective templates subdirectories for reference. 

2. Execute the bash script ``get_input_files.sh``::

      ./get_input_files.sh

  It downloads the input data files for the several model domains (Antarctica, Greenland, etc.) These files are stored on a server (`Zenodo archive <https://doi.org/10.5281/zenodo.6371122>`__) and needed for various inputs such as topography, precipitation, geothermal heat flux, etc. The script can be configured before execution if the input files are only needed for selected domains (default is downloading everything). To do so, open it with a text editor and change the flag variables according to the instructions in the script.

3. Locate the file ``sico_configs.sh`` in the main directory, and open it with a text editor.

   Set the flags ``LIS_FLAG``, ``OPENMP_FLAG`` and ``LARGE_DATA_FLAG`` according to your needs. 

   Default is ``"true"``/``"true"``/``"true"``. ``LIS_FLAG`` and ``OPENMP_FLAG`` can be set to ``"false"`` for simulations with pure shallow-ice dynamics. However, ``"true"`` is required for simulations with shallow-shelf dynamics (for floating ice) or hybrid shallow-ice--shelfy-stream dynamics (for grounded ice). The test simulations included in the repository all work with ``LARGE_DATA_FLAG`` set to ``"false"``. However, for higher-resolution simulations (e.g., Greenland/5 km or Antarctica/8 km), ``"true"`` is required.

   Set ``NETCDFHOME`` to the correct path of your NetCDF installation.

   If ``LIS_FLAG="true"``, set ``LISHOME`` to the correct path of your Lis installation.

   Depending on your system, some additional settings might have to be added in ``sico_configs.sh`` (``module load`` commands for dynamic loading etc.).

4. Locate the file ``sico_environment.sh`` in the main directory, open it with a text editor, and replace the ``SICO_INSTITUTION="Default"`` entry by the name of your institution (max. 256 characters).

Directory structure
===================

Main directory
  Initialization scripts ``copy_templates.sh``, ``get_input_files.sh``.

  Configuration scripts ``sico_configs.sh``, ``sico_environment.sh``.

  Shell script ``sico.sh`` for running a single simulation.

  Shell scripts ``multi_sico_1.sh`` and ``multi_sico_2.sh`` for running multiple simulations by repeated calls of ``sico.sh``.

  README and LICENSE files.

Directory ``headers``
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

  .. note::
    These subdirectories also contain README files that describe the input data and provide the corresponding references.
  
Directory ``sico_out``
  Directory into which output files of SICOPOLIS simulations are written by default.

Directory ``docs``
  Documentation for SICOPOLIS.

Directory ``tmp``
  Empty directory, contents ignored by Git.

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

.. _getting_started-run_simulation:

How to run a simulation
=======================

For example, to run the EISMINT Phase 2 Simplified Geometry Experiment A (Payne et al. :cite:`payne_etal_2000`), named ``repo_emtp2sge25_expA``, use the script ``sico.sh`` as follows::

  (./sico.sh -m repo_emtp2sge25_expA) >tmp/out_001.dat 2>&1 &

(from the main directory, bash required). Accordingly for any other simulation.

To list further options, execute ``./sico.sh -h``.

.. warning::
  The name ``out_<run_name>.dat`` must not be used for the redirected output of ``sico.sh``. This name is reserved for the runtime output of SICOPOLIS itself. (Both are very useful in case of compilation or runtime errors!)

If you prefer to run :ref:`all EISMINT, Antarctica and Greenland simulations <test_simulations>` consecutively, execute the script ``multi_sico_1.sh``::

  (./multi_sico_1.sh) >tmp/out_multi_100.dat 2>&1 &

To list further options, execute ``./multi_sico_1.sh -h``.

Alternatively, :ref:`all other test simulations (Austfonna etc.) <test_simulations>` can be run with the script ``multi_sico_2.sh``::

  (./multi_sico_2.sh) >tmp/out_multi_200.dat 2>&1 &

Approximate computing times are listed in the ":ref:`Test simulations <test_simulations>`" section.

.. _getting_started-output:

Output files
============

Output files are written by default to the directory ``sico_out/<run_name>``. This can be changed by executing ``sico.sh`` (or ``multi_sico_*.sh``) with the option ``-d /path/to/output/directory``. Four types are produced:

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
  * (Some more in the NetCDF file, execute ``ncdump -h <run_name>_ser.nc`` for a listing)

``<run_name>.core``, ``<run_name>_core.nc``:
  Time-series files (ASCII, NetCDF) that contain for selected locations xxx:

  * Time, t
  * Surface temperature anomaly, D\_Ts, or glacial index, glac\_ind (forcing)
  * Sea level, z\_sl (forcing)
  * Thickness, H\_xxx
  * Surface velocity, v\_xxx
  * Basal temperature, T\_xxx
  * (Some more in the NetCDF file, execute ``ncdump -h <run_name>_core.nc`` for a listing)

  | For the Greenland ice sheet, these data are written for seven locations:
  | GRIP (xxx=GR), GISP2 (xxx=G2), Dye 3 (xxx=D3), Camp Century (xxx=CC), NorthGRIP (xxx=NG), NEEM (xxx=NE), EastGRIP (xxx=EG).

  | For the Antarctic ice sheet, these data are written for six locations:
  | Vostok (xxx=Vo), Dome A (xxx=DA), Dome C (xxx=DC), Dome F (xxx=DF), Kohnen (xxx=Ko), Byrd (xxx=By).

``<run_name>0001.nc``, ``<run_name>0002.nc``, ...:
  Complete set of fields (topography, velocity, temperature etc., written in NetCDF format) for selected time slices.

Writing of output files can be controlled by the several parameters in the "Data output" section of the run-specs headers. For example, simulation ``repo_emtp2sge25_expA`` writes scalar variables into the time-series files ``repo_emtp2sge25_expA{.ser,.core,_ser.nc,_core.nc}`` every 100 years. In addition, it produces three time-slice files ``repo_emtp2sge25_expA0001.nc``, ``repo_emtp2sge25_expA0002.nc`` and ``repo_emtp2sge25_expA0003.nc``, which correspond to the times :math:`t=5\,\mathrm{ka}`, :math:`50\,\mathrm{ka}` and :math:`200\,\mathrm{ka}`, respectively.

.. note::
  By default, when trying to re-run a simulation, already existing output will not be overwritten, and an error message will be produced. However, overwriting can be enforced by executing ``sico.sh`` (or ``multi_sico_*.sh``) with the option ``-f``.

.. note::
  If a time-slice file of the initial state of a simulation is written, not all variables will already be defined (because SICOPOLIS has not done any proper computation yet). For instance, "diagnosed" 2D fields like the basal temparatures ``temp_b`` and ``temph_b`` (relative to pressure melting) or the thermal state mask ``n_cts`` will contain only default values. They will be filled with meaningful values after the first time step.
