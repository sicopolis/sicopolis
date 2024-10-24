.. _ice_dynamics:

Ice dynamics
************

The dynamics (force balance) of grounded ice can be chosen in the run-specs header by the parameter ``DYNAMICS``\:

* ``0``: Ice flow velocity set to zero everywhere (static ice).

* ``1``: Shallow-ice approximation (SIA) (e.g., Greve and Blatter :cite:`greve_blatter_2009`).

* ``2``: Hybrid shallow-ice--shelfy-stream (SIA/SStA) dynamics (Bernales et al. :cite:`bernales_etal_2017a, bernales_etal_2017b`).

* ``3``: DIVA (depth-integrated viscosity approximation) dynamics (Goldberg :cite:`goldberg_2011`, Lipscomb et al. :cite:`lipscomb_etal_2019`).

The parameter ``MARGIN`` controls how the grounded ice can extend into the sea:

* ``1``: Ice extent strictly restricted to land area ("hydrophobic mode").

* ``2``: Formation of marine ice possible.

* ``3``: Formation of marine ice and ice shelves possible.

In case of ``MARGIN = 3``, ice shelves (floating ice) are included and modelled by the shallow-shelf approximation (SSA; e.g., Greve and Blatter :cite:`greve_blatter_2009`).

The detailed settings are controlled by additional parameters as described in the run-specs headers.
