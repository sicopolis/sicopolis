.. _calving:

Calving
*******

.. _calving_ice_shelves:

Ice-shelf calving
=================

The options for calving of ice shelves (floating ice) can be selected in the run-specs headers by the parameter ``ICE_SHELF_CALVING``\:

* ``1``: Unlimited expansion of ice shelves, no calving.

* ``2``: Instantaneous calving of ice shelves if the thickness is less than a threshold thickness, specified by the parameter ``H_CALV``.

* ``3``: "Float-kill": Instantaneous removal of all floating ice.

For the Antarctic ice sheet, yearly ISMIP6-type ice-shelf collapse masks can be prescribed (Seroussi et al. :cite:`seroussi_etal_2020`). This requires the setting ``ICE_SHELF_COLLAPSE_MASK 1`` and additional parameters as described in the run-specs headers.

.. _calving_marine_ice:

Marine-ice calving
==================

For calving of grounded marine ice, the following options are available\:

* | Parameterization for "underwater-ice" calving (Dunse et al. :cite:`dunse_etal_2011`).
  | To be selected by the following combination of run-specs-header parameters\:
  | ``MARGIN 2``, ``MARINE_ICE_FORMATION 2``, ``MARINE_ICE_CALVING 9``.
  | Further parameters (``CALV_UW_COEFF``, ``R1_CALV_UW``, ``R2_CALV_UW``) as described in the run-specs headers.

For the Greenland ice sheet, yearly ISMIP6-type retreat masks can be prescribed (Goelzer et al. :cite:`goelzer_etal_2020`). This requires the setting ``RETREAT_MASK 1`` and additional parameters as described in the run-specs headers.
