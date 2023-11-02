.. _clim_ocn_forcing:

Oceanic forcing
***************

.. _sea_level:

Sea level
=========

The sea level surrounding the simulated ice sheet determines the land area available for glaciation. The level :math:`z=0` corresponds to the mean sea level at present day, and the SICOPOLIS variable ``z_sl`` denotes the sea level relative to this reference.

Two different options for prescribing the sea level are available, selected in the run-specs headers by the parameter ``SEA_LEVEL``:

* ``1``: Temporally constant sea level z_sl, specified by the parameter ``Z_SL0``.

* ``3``: Time-dependent sea level (e.g., reconstruction from data) from an input file (ASCII or NetCDF), specified by the parameter ``SEA_LEVEL_FILE``.

Both options assume a spatially constant sea level. However, the variable ``z_sl`` is actually a 2D field, so that SICOPOLIS can handle in principle a spatially variable sea level as well.

.. _ice_shelf_basal_melting:

Ice-shelf basal melting
=======================

Lorem ipsum...

.. Move calving to here?
