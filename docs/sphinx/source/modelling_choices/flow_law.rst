.. _flow_law:

Flow law
********

Several options for the flow law of polycrystalline ice are available. They can be selected by the parameter ``FLOW_LAW`` in the run-specs headers\:

* ``1``: Glen's flow law (Glen :cite:`glen_1955`, Nye :cite:`nye_1957`) with stress exponent :math:`n=3`.

* ``2``: Goldsby-Kohlstedt :cite:`goldsby_kohlstedt_1997, goldsby_kohlstedt_2001` flow law with stress exponent :math:`n=1.8` and grain-size exponent :math:`p=1.4`. Average grain size defined by the parameter ``GR_SIZE``.

* ``3``: Flow law by Durham et al. :cite:`durham_etal_1997a` with stress exponent :math:`n=4`.

* ``4``: Polynomial flow law by Smith and Morland :cite:`smith_morland_1981` (summarized by Greve and Blatter :cite:`greve_blatter_2009`, Section 4.3.3).

For the cases ``FLOW_LAW = 1, 2 or 3``, the additional parameter ``FIN_VISC`` allows choosing between the unmodified flow law with an infinite-viscosity limit for low strain rates (``FIN_VISC = 1``), or using a regularized flow law with a finite-viscosity limit (``FIN_VISC = 2``). The latter is defined by a non-vanishing residual stress :math:`\sigma_0` (parameter ``SIGMA_RES``; see Greve and Blatter :cite:`greve_blatter_2009`, Section 4.3.2).

.. note::
  The rate factor :math:`A(T')` must fit the flow law unit- and value-wise. It is defined in the :ref:`physical-parameter files <getting_started-phys_para>` as a list for integer temperature values (between which linear interpolation is applied).
