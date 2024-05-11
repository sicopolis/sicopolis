.. _clim_ocn_forcing:

Oceanic forcing
***************

.. |nbsp| unicode:: 0xA0 
   :trim:

.. _sea_level:

Sea level
=========

The sea level surrounding the simulated ice sheet determines the land area available for glaciation. The level :math:`z=0` corresponds to the mean sea level at present day, and the SICOPOLIS variable ``z_sl`` denotes the sea level relative to this reference.

Two different options for prescribing the sea level are available, selected in the run-specs headers by the parameter ``SEA_LEVEL``\:

* ``1``: Temporally constant sea level z_sl, specified by the parameter ``Z_SL0``.

* ``3``: Time-dependent sea level (e.g., reconstruction from data) from an input file (ASCII or NetCDF), specified by the parameter ``SEA_LEVEL_FILE``.

Both options assume a spatially constant sea level. However, the variable ``z_sl`` is actually a 2D field, so that SICOPOLIS can handle in principle a spatially variable sea level as well.

.. _ice_shelf_basal_melting:

Ice-shelf basal melting
=======================

The parameter ``FLOATING_ICE_BASAL_MELTING`` in the run-specs headers allows specifying the melting rate under ice shelves (floating ice), :math:`a_\mathrm{b}`. For all terrestrial ice sheets, the following options can be chosen\:

* ``1``: Constant values for the continental shelf (:math:`a_\mathrm{b}^\mathrm{c.s.}`) and the abyssal ocean (:math:`a_\mathrm{b}^\mathrm{a.o.}`), respectively\:

  .. math::
    :label: eq_ice_shelf_bas_melt_1

    a_\mathrm{b}
    = \left\{
      \begin{array}{ll}
      a_\mathrm{b}^\mathrm{c.s.} & \mbox{if}\;\; z_\mathrm{l} > z_\mathrm{abyss}\,,
      \\
      a_\mathrm{b}^\mathrm{a.o.} & \mbox{if}\;\; z_\mathrm{l} \le z_\mathrm{abyss}\,,
      \end{array}
    \right.

  where :math:`z_\mathrm{l}` is the seabed (lithosphere surface) elevation and :math:`z_\mathrm{abyss}` the threshold seabed elevation that separates the continental shelf from the abyssal ocean. The parameters :math:`a_\mathrm{b}^\mathrm{c.s.}`, :math:`a_\mathrm{b}^\mathrm{a.o.}` and :math:`z_\mathrm{abyss}` can be set in the run-specs headers (``QBM_FLOAT_1``, ``QBM_FLOAT_3`` and ``Z_ABYSS``, respectively).

* ``4``: Local parameterization as a function of the oceanic thermal forcing :math:`T_\mathrm{f}=T_\mathrm{oc}-T_\mathrm{b}` (difference between the ocean temperature :math:`T_\mathrm{oc}` and the ice-shelf basal temperature :math:`T_\mathrm{b}`)\:

  .. math::
    :label: eq_ice_shelf_bas_melt_4

    a_\mathrm{b} = \Omega\,T_\mathrm{f}^\alpha\,.

  The parameters :math:`T_\mathrm{oc}`, :math:`\Omega` and :math:`\alpha` can be set in the run-specs headers (``TEMP_OCEAN``, ``OMEGA_QBM`` and ``ALPHA_QBM``, respectively). The ice-shelf basal temperature is computed as

  .. math::
    :label: eq_ice_shelf_bas_temp

    T_\mathrm{b} = -\beta_\mathrm{sw} d 
                     - \Delta{}T_\mathrm{m,sw}\,,

  where :math:`d` is the draft (depth of the ice-shelf base below sea level), :math:`\beta_\mathrm{sw}=7.61\times{}10^{-4}\,\mathrm{K\,m^{-1}}` the Clausius-Clapeyron gradient and :math:`\Delta{}T_\mathrm{m,sw}=1.85^\circ\mathrm{C}` the melting-point lowering due to the average salinity of sea water.

For the Antarctic ice sheet, two additional options are available\:

* ``5``: Sector-wise, local parameterization as a function of the thermal forcing (Greve and Galton-Fenzi :cite:`greve_galton-fenzi_2017`). This parameterization is modified after Beckmann and Goosse :cite:`beckmann_goosse_2003`, with a linear dependence on the thermal forcing :math:`T_\mathrm{f}` and an additional power-law dependence on the draft :math:`d`\:

  .. math::
    :label: eq_ice_shelf_bas_melt_5

    a_\mathrm{b} 
      = \frac{\rho_\mathrm{sw}c_\mathrm{sw}\gamma_\mathrm{t}}{\rho L}
        \,\Omega\,\bigg(\frac{d}{d_0}\bigg)^\alpha \,T_\mathrm{f}\,,

  where :math:`\rho` and :math:`\rho_\mathrm{sw}` are the ice and sea-water densities, :math:`L` is the latent heat of melting (all defined in the run-specs header), :math:`c_\mathrm{sw}=3974\,\mathrm{J\,kg^{-1}\,K^{-1}}` is the specific heat of sea water, :math:`\gamma_\mathrm{t}=5\times{}10^{-5}\,\mathrm{m\,s^{-1}}` is the exchange velocity for temperature and :math:`d_0=200\,\mathrm{m}` is the reference draft.

  The parameters :math:`\Omega` and :math:`\alpha` result from a tuning procedure for eight different sectors, using observed present-day melt rates as a target (as explained in the main part and appendix of Greve and Galton-Fenzi :cite:`greve_galton-fenzi_2017`). For the thermal forcing :math:`T_\mathrm{f}`, :math:`T_\mathrm{oc}` is chosen for each sector as the sector-averaged temperature at 500 metres depth just outside the ice-shelf cavity (computed with data from the World Ocean Atlas 2009 :cite:`locarnini_etal_2010`), while :math:`T_\mathrm{b}` is computed by Eq. |nbsp| :eq:`eq_ice_shelf_bas_temp`.

* ``6``: "ISMIP6 standard approach": Sector-wise, non-local quadratic parameterization for the 18 IMBIE-2016 sectors (Rignot and Mouginot :cite:`rignot_mouginot_2016`, The IMBIE Team :cite:`imbie_2018`), where the two sectors feeding the Ross ice shelf and the two sectors feeding the Filchner--Ronne ice shelf are combined, leaving 16 distinct sectors (Jourdain et al. :cite:`jourdain_etal_2020`, Seroussi et al. :cite:`seroussi_etal_2020`). The parameterization depends on the local thermal forcing :math:`T_\mathrm{f}` and the sector-averaged thermal forcing :math:`\langle{}T_\mathrm{f}\rangle{}_\mathrm{sector}` as follows\:

  .. math::
    :label: eq_ice_shelf_bas_melt_6

    a_\mathrm{b} 
      = \gamma_0
        \bigg(\frac{\rho_\mathrm{sw}c_\mathrm{sw}}{\rho L}\bigg)^2
        \, (T_\mathrm{f} + \delta{}T_\mathrm{sector})
        \, |\langle{}T_\mathrm{f}\rangle{}_\mathrm{sector}
              + \delta{}T_\mathrm{sector}|\,,

  where :math:`\rho`, :math:`\rho_\mathrm{sw}`, :math:`L` and :math:`c_\mathrm{sw}` are defined as in Eq. |nbsp| :eq:`eq_ice_shelf_bas_melt_5`. The coefficient :math:`\gamma_0`, similar to an exchange velocity, and the sectorial temperature offsets :math:`\delta{}T_\mathrm{sector}` are obtained by calibrating the parameterization against observations (see Jourdain et al. :cite:`jourdain_etal_2020`).

  The thermal forcing at the ice--ocean interface is derived by extrapolating the oceanic fields from GCMs into the ice-shelf cavities. Following the ISMIP6-Antarctica protocol, it must be provided as NetCDF input files that contain for each year the mean-annual, 3D thermal forcing for the entire computational domain. Thereby, this option allows prescribing a time-dependent thermal forcing (which is currently not the case for the other options). For the detailed parameter settings, see the description in the run-specs headers.

For all cases, an additional scaling factor :math:`S_\mathrm{w}` can be applied (:math:`a_\mathrm{b}\rightarrow{}S_\mathrm{w}\,a_\mathrm{b}`), defined as

.. math::
  :label: eq_ice_shelf_bas_melt_scaling_factor

  S_\mathrm{w}
    = \mathrm{tanh}\,\bigg(\frac{H_\mathrm{w}}{H_\mathrm{w,0}}\bigg)\,.

This factor reduces the melting rate close to the grounding line where the water column :math:`H_\mathrm{w}` is thin. The parameter :math:`H_\mathrm{w,0}` can be set in the run-specs headers (``H_W_0``). A value recommended by Asay-Davis et al. :cite:`asay-davis_etal_2016` is :math:`75\,\mathrm{m}`, while Gladstone et al. :cite:`gladstone_etal_2017` used :math:`36.79\,(=100/e)\,\mathrm{m}`. Setting this parameter to zero results in :math:`S_\mathrm{w}=1` everywhere; the scaling is then switched off.

  .. _calving_ice_shelves:

Ice-shelf calving
=================

The options for calving of ice shelves (floating ice) can be selected in the run-specs headers by the parameter ``ICE_SHELF_CALVING``\:

* ``1``: Unlimited expansion of ice shelves, no calving.

* ``2``: Instantaneous calving of ice shelves if the thickness is less than a threshold thickness, specified by the parameter ``H_CALV``.

* ``3``: "Float-kill": Instantaneous removal of all floating ice.

For the Antarctic ice sheet, yearly ISMIP6-type ice-shelf collapse masks can be prescribed (Seroussi et al. :cite:`seroussi_etal_2020`). This requires the setting ``ICE_SHELF_COLLAPSE_MASK = 1`` and additional parameters as described in the run-specs headers.

.. _calving_marine_ice:

Marine-ice calving
==================

For calving of grounded marine ice, the following options are available\:

* Parameterization for "underwater-ice" calving (Dunse et al. :cite:`dunse_etal_2011`), to be selected by the following combination of run-specs-header parameters\: ``MARGIN = 2``, ``MARINE_ICE_FORMATION = 2``, ``MARINE_ICE_CALVING = 9``. This parameterization is an adaption of the law by Clarke et al. :cite:`clarke_etal_1999`, but acts here as an additional surface ablation rather than calving at a vertical front\:

  .. math::
    :label: eq_uw_ice_calving

    Q_\mathrm{c} = k_\mathrm{c} H^{r_1} D_\mathrm{w}^{r_2}\,,

  where :math:`Q_\mathrm{c}` is the calving flux, :math:`H` the ice thickness (taken to some power :math:`r_1`), :math:`D_\mathrm{w}` the water depth (taken to some power :math:`r_2`) and :math:`k_\mathrm{c}` the calving parameter (see also :numref:`uw_ice_calving`). The two exponents and the calving parameter can be set in the run-specs headers as parameters ``R1_CALV_UW``, ``R2_CALV_UW`` and ``CALV_UW_COEFF``, respectively.

  .. _uw_ice_calving:
  .. figure:: figs/UW_Ice_Calving.png
    :width: 500 px
    :alt: Underwater ice calving
    :align: center

    Schematic of underwater ice calving. The purple area marks the marine grounded ice, the white area the "underwater ice" (fulfilling the floating condition) for which the calving law (Eq. |nbsp| :eq:`eq_uw_ice_calving`) is applied.

For the Greenland ice sheet, yearly ISMIP6-type retreat masks can be prescribed (Goelzer et al. :cite:`goelzer_etal_2020`). This requires the setting ``RETREAT_MASK = 1`` and additional parameters as described in the run-specs headers.
