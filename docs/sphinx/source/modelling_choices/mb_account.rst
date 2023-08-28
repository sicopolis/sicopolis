.. _mb_account:

Mass balance accounting
***********************

The global mass balance of an ice sheet reads

.. math::
  :label: global_mb

  \frac{\mathrm{d}V}{\mathrm{d}t}
    = \mathrm{SMB} + \mathrm{BMB} + \mathrm{CALV}\,,

where :math:`\mathrm{d}V/\mathrm{d}t` is the rate of volume change, SMB the total surface mass balance, BMB the total basal mass balance and CALV the total calving rate (all counted as positive for a volume gain and negative for a volume loss). SICOPOLIS attempts at closing this balance as accurately as possible by applying the so-called "hidden ablation scheme" (Calov et al. :cite:`calov_etal_2015`).

SICOPOLIS always employs a zero-ice-thickness boundary condition at the margin of the computational domain (:math:`i=0,\,i_\mathrm{max}`; :math:`j=0,\,j_\mathrm{max}`). However, for accurate accounting of calving near the margin, the "hidden ablation scheme" also requires the next lines of grid points (:math:`i=1,\,i_\mathrm{max}-1`; :math:`j=1,\,j_\mathrm{max}-1`) to be ice-free. Since this is not always desirable, it can be controlled by the run-specs-header parameter ``MB_ACCOUNT``\:

* ``0``: Glaciation of all inner points of the domain allowed (prevents accurate accounting of calving near the margin).

* ``1``: Outermost inner points of the domain (:math:`i=1,\,i_\mathrm{max}-1`; :math:`j=1,\,j_\mathrm{max}-1`) not allowed to glaciate (required for accurate accounting of calving near the margin).

For real-world problems, the setting ``MB_ACCOUNT = 1`` is usually fine. However, for some simple-geometry experiments that require the simulated ice sheet to cover the entire domain [e.g., the :ref:`test simulation <test_simulations>` ``repo_vialov3d25`` (3D Vialov profile)], ``MB_ACCOUNT = 0`` must be chosen to allow that.
