.. _ice_thickness_evolution:

Ice thickness evolution
***********************

Lorem ipsum...

..
  #define THK_EVOL 1
                      0 : No evolution of the ice thickness, kept fixed on
                          the initial thickness
                      1 : Evolution of the ice thickness
                      2 : Evolution of the ice thickness, but between times
                          TIME_TARGET_TOPO_INIT0 and TIME_TARGET_TOPO_FINAL0
                          the ice topography (zs, zb, zl, H) is nudged
                          towards a prescribed target with the
                          relaxation time smoothly decreasing from
                          TARGET_TOPO_TAU0 to zero.
                      3 : Evolution of the ice thickness, but
                          the ice topography (zs, zb, zl, H) is nugded
                          towards a prescribed target with the
                          constant relaxation time TARGET_TOPO_TAU0.
                      4 : Evolution of the ice thickness,
                          but maximum ice extent is constrained by the
                          prescribed mask MASK_MAXEXTENT_FILE.
