.. _initial_conditions:

Initial conditions
******************

.. _initial_conditions_topography:

Topography
==========

The initial conditions for the ice topography can be chosen in the run-specs headers by the parameter ``ANF_DAT`` ("ANF"
being a remaining Germanicism in the code, derived from the German word "Anfang" = start, beginning, initial)\:

* ``1``: Present-day initial topography.

* ``2``: Ice-free initial topography with relaxed lithosphere.
  
* ``3``: Initial values from previous simulation.

For specifying the topography files, see Section ":ref:`Model domain, grid and time/Topography <topography>`".

.. _initial_conditions_temperature:

Temparature
===========

In case of ``ANF_DAT 1``, initial conditions for the 3D temperature field must also be prescribed. This is done by the parameter ``TEMP_INIT``\:

* ``1``: Constant value in the entire ice sheet (defined by parameter ``TEMP_INIT_VAL``).

* ``2``: In each ice column equal to the local ice surface temperature.

* ``3``: Ice temperature linearly increasing with depth.

* ``4``: Ice-temperature profiles from analytical solution for 1D steady-state advection-diffusion equation under the assumption of linearly varying :math:`v_z` (Robin :cite:`robin_1955` solution).

* ``5``: Ice temperature from previous simulation.

The settings ``ANF_DAT 3``, or ``ANF_DAT 1`` and ``TEMP_INIT 5``, require specifying an initial-value file, which is typically a time-slice output file from a previous simulation (see Section ":ref:`getting_started-output`"). This is done by assigning the parameter ``ANFDATNAME`` a text string with the file name. In addition, the path where this file is located must be given as an option ``-a /path/to/initial/value/directory`` to the script ``sico.sh`` when running SICOPOLIS (see Section ":ref:`getting_started-run_simulation`").
