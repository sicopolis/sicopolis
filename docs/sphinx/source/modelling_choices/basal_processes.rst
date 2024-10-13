.. _subglacial_processes:

Subglacial processes
********************

.. _basal_sliding:

Basal sliding
=============

Two options for basal sliding are implemented, to be chosen in the run-specs header by the parameter ``SLIDE_LAW``\:

* ``0``: No-slip conditions.

* ``1``: Weertman-Budd-type sliding law (Weertman :cite:`weertman_1957a`, Budd et al. :cite:`budd_etal_1979`, Budd and Jenssen :cite:`budd_jenssen_1987`).

Weertman-Budd-type sliding depends via the reduced pressure (ice minus water) on the basal water pressure, which can be chosen by the parameter ``BASAL_WATER_PRESSURE`` as follows\:

* ``0``: Zero everywhere.

* ``1``: Ocean pressure without cut-off (can become negative).

* ``2``: Ocean pressure with cut-off (non-negative values enforced).

Sub-melt sliding (Sato and Greve :cite:`sato_greve_2012`), water-film-enhanced sliding (requires ``BASAL_HYDROLOGY = 1``, see ":ref:`basal_hydrology`" below) and regionally varying sliding parameters can be added. The detailed settings are controlled by additional parameters as described in the run-specs headers.

.. _ghf:

Geothermal heat flux
====================

The geothermal heat flux (GHF), assumed to be time-independent, can be specified in the run-specs headers as either a constant value or a spatially varying distribution via the parameters ``Q_GEO`` and ``Q_GEO_FILE``\:

* If ``Q_GEO_FILE = 'none'`` (or undefined): Constant GHF defined by parameter ``Q_GEO``.

* Otherwise: Spatially varying GHF read from file specified by ``Q_GEO_FILE``.

If a file with gridded data is provided (second case), it must match the chosen horizontal grid (see ":ref:`spatial_grid`"). The format can either be NetCDF (``*.nc``) or ASCII (any other file extension).

A further, relevant parameter is ``Q_LITHO``\:

* ``0``: No coupled heat-conducting bedrock.

* ``1``: Coupled heat-conducting bedrock.

If set to ``0``, the GHF is imposed directly at the grounded ice base, which is suitable for steady-state simulations because it reduces the time required to reach the steady state. However, for transient simulations, ``1`` is the preferred setting. The GHF is then imposed at the base of the lithosphere layer (thickness defined by the parameter ``H_R`` in the run-specs header), so that the thermal inertia of the lithosphere is properly accounted for.

.. _basal_hydrology:

Basal hydrology
===============

Basal hydrology can be selected in the run-specs header by the parameter ``BASAL_HYDROLOGY``\:

* If set to ``0``, basal hydrology is ignored.

* If set to ``1``, it is assumed that basal water exists and moves in a thin (order of millimetres) and distributed water film. The film thickness is computed by a steady-state routing scheme for subglacial water that receives its input from the basal melting rate under grounded ice (Le Brocq et al. :cite:`lebrocq_etal_2006, lebrocq_etal_2009`, Calov et al. :cite:`calov_etal_2018`). The computations are carried out by the module ``hydro_m``.

.. _gia:

Glacial isostatic adjustment
============================

Three options are available for glacial isostatic adjustment, which can be selected in the run-specs header by the parameter ``REBOUND``\:

* ``0``: Rigid lithosphere, no adjustment.

* ``1``: Local-lithosphere--relaxing-asthenosphere (LLRA) model.

* ``2``: Elastic-lithosphere--relaxing-asthenosphere (ELRA) model.

These models are described by Le Meur and Huybrechts :cite:`lemeur_huy_1996` and Greve :cite:`greve_2001`.

The detailed settings are controlled by additional parameters (``FRAC_LLRA``, ``TIME_LAG_MOD``, ``TIME_LAG``, ``TIME_LAG_FILE``, ``FLEX_RIG_MOD``, ``FLEX_RIG``, ``FLEX_RIG_FILE``, ``DTIME_WSS0``) as described in the run-specs headers.

.. note::
  The isostatically relaxed lithosphere surface topography (parameter ``ZL0_FILE``, see ":ref:`topography`") is required for the isostasy models. A special setting for generating this topography can be enabled by

  .. code-block:: fortran

    #define EXEC_MAKE_ZL0

  It should be used together with ``ANF_DAT = 1`` (present-day topography used as initial topography), computes the isostatically relaxed lithosphere surface topography, writes it on file and then stops the simulation (irrespective of the setting for the final time :math:`t_\mathrm{final}`). The underlying assumption is that the present-day bed topography is approximately in equilibrium with the present-day ice load.
