.. _enhancement_factor:

Flow enhancement factor
***********************

The flow enhancement factor :math:`E` is a multiplier for the temperature-dependent rate factor :math:`A(T')` that accounts for ice softening (:math:`E>1`) or hardening (:math:`E<1`) due to impurities or an anisotropic fabric. For the setting ``ENHMOD 1``, three different values can be defined in the run-specs headers\:

* ``ENH_FACT``: :math:`E_\mathrm{SIA}`, for the shallow-ice approximation (SIA) of grounded ice.

* ``ENH_STREAM``: :math:`E_\mathrm{SStA}`, for the shelfy-stream approximation (SStA) of grounded ice.

* ``ENH_SHELF``: :math:`E_\mathrm{SSA}`, for the shallow-shelf approximation (SSA) of floating ice.

Further settings of ``ENHMOD`` allow making the enhancement factor dependent on age or time of deposition, so that different values for glacial and interglacial ice can be defined. This affects only the SIA. See the documentation in the run-specs headers for details.

For computing the velocity field with hybrid SIA/SStA dynamics, :math:`E_\mathrm{SIA}` is used for the SIA solution, :math:`E_\mathrm{SStA}` is used for the SStA solution, and then the velocity field is obtained as a weighted average of the two solutions. However, for computing the temperature field (the dissipation term depends on the ice viscosity and thus on the enhancement factor), this is not possible because no separate SIA and SStA solutions are available. Therefore, a weighted enhancement factor :math:`E_\mathrm{SIA/SStA}` is computed from :math:`E_\mathrm{SIA}` and :math:`E_\mathrm{SStA}`, using the same weighting rule used for the velocity. This weighted enhancement factor is then used for the dissipation term of the temperature equation.
