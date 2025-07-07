.. _flow_law:

Flow law
********

Several options for the flow law of polycrystalline ice are available. They can be selected by the parameter ``FLOW_LAW`` in the run-specs headers\:

* ``1``: Nye-Glen flow law (Glen :cite:`glen_1955`, Nye :cite:`nye_1957`) with stress exponent :math:`n`.

* ``4``: Polynomial flow law by Smith and Morland :cite:`smith_morland_1981` (summarized by Greve and Blatter :cite:`greve_blatter_2009`, Section 4.3.3).

For the case ``FLOW_LAW = 1``, the stress exponent :math:`n` is defined by the parameter ``N_POWER_LAW``. The additional parameter ``FIN_VISC`` allows choosing between the unmodified flow law with an infinite-viscosity limit for low strain rates (``FIN_VISC = 1``), or using a regularized flow law with a finite-viscosity limit (``FIN_VISC = 2``). The latter is defined by a non-vanishing residual stress :math:`\sigma_0` (parameter ``SIGMA_RES``; see Greve and Blatter :cite:`greve_blatter_2009`, Section 4.3.2).

.. note::
  The rate factor :math:`A(T')` is defined as a list for integer temperature values, to be read from a file specified in the run-specs header (parameter ``RF_KAPPA_C_FILE``). Between integer temperatures, linear interpolation is applied. To avoid problems with varying units and numerical values, a dimensionless formulation is chosen (Greve :cite:`greve_2025`).
