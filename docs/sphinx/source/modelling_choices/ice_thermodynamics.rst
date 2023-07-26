.. _ice_thermodynamics:

Ice thermodynamics
******************

For modelling ice thermodynamics, five different options are available, which can be chosen in the run-specs header by the parameter ``CALCMOD``\:

* ``-1``: ISOT: isothermal method, constant temperature and age.

* ``0``: COLD: cold-ice method, resetting of temperatures above pressure melting.

* ``1``: POLY: polythermal method, separate domains for cold and temperate ice.

* ``2``: ENTC: conventional enthalpy method.

* ``3``: ENTM: melting-CTS enthalpy method.

Options ``1`` and ``3`` are the most sophisticated solvers and thus recommended for most real-world problems. The different methods are discussed by Greve and Blatter :cite:`greve_blatter_2016`.

As explained in Section ":ref:`spatial_grid`", for the polythermal method (POLY) separate vertical domains are employed for the cold-ice layer (:math:`\zeta_\mathrm{c}` domain) and the temperate-ice layer (:math:`\zeta_\mathrm{t}` domain). In all other cases, the :math:`\zeta_\mathrm{c}` domain is used for the entire ice column from the base to the surface, the :math:`\zeta_\mathrm{t}` domain is redundant, and the parameter ``KTMAX`` should be set to ``2``.

If the polythermal method is chosen (``CALCMOD = 1``), the parameter ``CTS_MELTING_FREEZING`` determines whether melting and freezing conditions are distinguished (``1``), or melting conditions are always assumed (``2``). In principle, the former is physically more adequate. However, in practice, freezing conditions occur only for a very small fraction of the areas with a temperate layer, and their treatment can cause numerical issues due to the discontinuities of the temperature gradient and the water content at the CTS (cold-temperate transition surface). Therefore, setting ``CTS_MELTING_FREEZING`` to ``2`` is safer and usually sufficient.
