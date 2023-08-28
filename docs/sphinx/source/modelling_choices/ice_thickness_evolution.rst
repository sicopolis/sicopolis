.. _ice_thickness_evolution:

Ice-thickness evolution
***********************

Several options are available for the evolution of the ice thickness, selected by the run-specs-header parameter ``THK_EVOL``\:

* ``0``: No evolution of the ice thickness, kept fixed on the initial thickness.

* ``1``: Free evolution of the ice thickness.

* ``2``: Evolution of the ice thickness, but the ice topography (zs, zb, zl, H) is nugded towards a prescribed target with a time-dependent relaxation time read from the file ``TARGET_TOPO_TAU0_FILE`` (to be located in ``sico_in/general``).

* ``3``: Evolution of the ice thickness, but the ice topography (zs, zb, zl, H) is nugded towards a prescribed target with the constant relaxation time ``TARGET_TOPO_TAU0``.

For the cases ``THK_EVOL >= 1``, the ice-thickness equation

.. math::
  :label: ice_thickness_equation

  \frac{\partial{}H}{\partial{}t}
    = -\left( \frac{\partial{}Q_x}{\partial{}x}
              + \frac{\partial{}Q_y}{\partial{}y} \right)
      + a_\mathrm{s} + a_\mathrm{b}

is solved, where :math:`H` is the ice thickness, :math:`(Q_x,Q_y)` the volume flux (depth-integrated horizontal velocity), :math:`a_\mathrm{s}` the surface mass balance (SMB) and :math:`a_\mathrm{b}` the basal mass balance (BMB). Note that SMB and BMB are consistently counted as positive for a volume gain and negative for a volume loss.

For the cases ``THK_EVOL = 2, 3``, nudging of the ice topography towards a prescribed target topography (e.g., a slightly smoothed present-day topography) is employed. During each time step, after solving the ice-thickness equation, the relaxation equations

.. math::
  :label: topo_nudging

  \frac{\partial{}h}{\partial{}t} 
    = -\frac{h-h_\mathrm{target}}{\tau}\,,
  \quad
  \frac{\partial{}b}{\partial{}t} 
    = -\frac{b-b_\mathrm{target}}{\tau}\,,
  \quad
  \frac{\partial{}z_\mathrm{l}}{\partial{}t} 
    = -\frac{z_\mathrm{l}-z_\mathrm{l,target}}{\tau}\,,

are applied, where :math:`h` is the surface topography, :math:`b` the ice-base topography (:math:`H=h-b`), :math:`z_\mathrm{l}` the lithosphere-surface (bedrock/seabed) topography, :math:`(\cdot)_\mathrm{target}` the respective nudging target and :math:`\tau` the relaxation time (Rueckamp et al. :cite:`rueckamp_etal_2019a`). The nudging target is specified by the parameter ``TARGET_TOPO_DAT_NAME``, to be assigned a text string with the file name (typically a time-slice output file from a previous simulation, see Section ":ref:`getting_started-output`"). In addition, the path where this file is located must be given as an option ``-t /path/to/nudging/target/directory`` to the script ``sico.sh`` when running SICOPOLIS (see Section ":ref:`getting_started-run_simulation`").

Nudging is equivalent to applying an SMB correction, which is diagnosed by the model.

For the cases ``THK_EVOL >= 1``, the maximum ice extent can be constrained spatially by a prescribed mask file, specified by the parameter ``MASK_MAXEXTENT_FILE`` (file to be located in ``sico_in/ant`` for Antarctica, ``sico_in/grl`` for Greenland, etc.). If set to ``'none'`` or undefined, no spatial constraint is applied.

The numerical scheme for solving the ice-thickness equation can be chosen by the parameter ``CALCTHK``\:

  * ``1``: Explicit scheme for the diffusive SIA ice-surface equation.

  * ``2``: Over-implicit scheme for the diffusive SIA ice-surface equation, iterative built-in SOR solver.

  * ``4``: Explicit scheme for the general ice-thickness equation.

For more details on the cases ``CALCTHK = 1, 2``, see Greve and Calov :cite:`greve_calov_2002`. However, these cases are only applicable for SIA dynamics (see ":ref:`ice_dynamics`"), and they do not preserve mass well on complex bed topographies. Therefore, for most real-world problems, ``CALCTHK = 4`` is the recommended setting.
