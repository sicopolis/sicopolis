.. _test_simulations:

Test simulations
****************

These are a number of computationally rather inexpensive test simulations, of which the :ref:`run-specs header files <getting_started-run_specs_headers>` are contained in the SICOPOLIS repository.

Run ``repo_vialov3d25``
  | 3D version of the 2D Vialov profile (Vialov :cite:`vialov_1958`),
  | SIA, resolution 25 km, :math:`t=0\ldots{}100\,\mathrm{ka}`.
  | Similar to the EISMINT Phase 1 fixed-margin experiment (Huybrechts et al. :cite:`huybrechts_etal_1996`), but without thermodynamics. Instead, isothermal conditions with :math:`T=-10^{\circ}\mathrm{C}` everywhere are assumed.

Run ``repo_emtp2sge25_expA``
  | EISMINT Phase 2 Simplified Geometry Experiment A,
  | SIA, resolution 25 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Payne et al. :cite:`payne_etal_2000`).
  | The thermodynamics solver for this run is the one-layer melting-CTS enthalpy scheme (ENTM), while all other runs employ the polythermal two-layer scheme (POLY) (Greve and Blatter :cite:`greve_blatter_2016`).

Runs ``repo_grl16_bm5_{init100a, ss25ka}``
  | Greenland ice sheet, SIA, resolution 16 km;
  | :math:`t=-100\,\mathrm{a}\ldots{}0` for the init run without basal sliding (..._init100a),
  | :math:`t=0\,\ldots{}25\,\mathrm{ka}` for the main run (..._ss25ka),
  | steady-state run for modern climate conditions, free evolution during the first 10 ka, after that gradual nudging towards the slightly smoothed present-day topography computed by the init run (unpublished).

Run ``repo_ant40_b2_ss25ka``
  | Antarctic ice sheet without ice shelves, SIA, resolution 40 km,
  | short steady-state run (:math:`t=0\ldots{}25\,\mathrm{ka}`) for modern climate conditions (unpublished).

Run ``repo_grl20_b2_paleo21``
  | Greenland ice sheet, SIA, resolution 20 km,
  | :math:`t=-140\,\mathrm{ka}\ldots{}0`, basal sliding ramped up during the first 5 ka.
  | Modified, low-resolution version of the spin-up for ISMIP6 InitMIP (Greve et al. :cite:`greve_etal_2017a`).

Runs ``repo_grl10_b2_{paleo21, future21_ctrl, future21_asmb}``
  | Greenland ice sheet, SIA, resolution 10 km,
  | :math:`t=-9\,\mathrm{ka}\ldots{}0` for the paleo run, :math:`t=0\ldots{}100\,\mathrm{a}` for the future runs.
  | 10-km version of the spin-up and the schematic future climate runs for ISMIP6 InitMIP (Greve et al. :cite:`greve_etal_2017a`).

Runs ``repo_ant64_b2_{spinup09_init100a, spinup09_fixtopo, spinup09, future09_ctrl}``
  | Antarctic ice sheet with hybrid shallow-ice--shelfy-stream dynamics
  | (Bernales et al. :cite:`bernales_etal_2017a`) and ice shelves (SSA), resolution 64 km;
  | :math:`t=-140.1\ldots{}-140\,\mathrm{ka}` for the init run without basal sliding (..._init100a),
  | :math:`t=-140\,\mathrm{ka}\ldots{}0` for the run with almost fixed topography (..._fixtopo), basal sliding ramped up during the first 5 ka,
  | :math:`t=-0.5\,\mathrm{ka}\ldots{}0` for the final, freely-evolving-topography part of the (..._spinup09),
  | :math:`t=0\ldots{}100\,\mathrm{a}` for the constant-climate control run (..._future09_ctrl).
  | 64-km version of the spin-up and the constant-climate control run for ISMIP6 InitMIP (Greve and Galton-Fenzi, pers. comm. 2017).

Runs ``repo_asf2_steady``, ``repo_asf2_surge``
  | Austfonna, SIA, resolution 2 km, :math:`t=0\ldots{}10\,\mathrm{ka}`.
  | Similar to Dunse et al. :cite:`dunse_etal_2011`'s Exp. 2 (steady fast flow) and Exp. 5 (surging-type flow), respectively.

Runs ``repo_nmars10_steady``, ``repo_smars10_steady``
  | North-/south-polar cap of Mars, SIA, resolution 10 km, :math:`t=-10\,\mathrm{Ma}\ldots{}0`.
  | Steady-state runs by Greve :cite:`greve_2007b`.
 
Run ``repo_nhem80_nt012_new``
  | Northern hemisphere, SIA, resolution 80 km, :math:`t=-250\,\mathrm{ka}\ldots{}0`.
  | Similar to run nt012 by Greve et al. :cite:`greve_etal_1999a`.

Run ``repo_heino50_st``
  | ISMIP HEINO standard run ST, SIA, resolution 50 km, :math:`t=0\ldots{}200\,\mathrm{ka}` (Calov et al. :cite:`calov_etal_2010`).

-------------

**Model times, time steps, computing times:**

+-------------------------------------+------------+---------------------+--------------------+
| Run                                 | Model time | Time step\ :sup:`†` | CPU time\ :sup:`‡` |
+=====================================+============+=====================+====================+
| repo\_vialov3d25                    | 100 ka     | 20 a                | 1.0 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_emtp2sge25\_expA              | 200 ka     | 20 a                | 3.9 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl16\_bm5\_100a              | 100 a      | 5 a                 | 7.1 sec            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl16\_bm5\_ss25ka            | 25 ka      | 5 a                 | 9.9 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_ant40\_b2\_ss25ka             | 25 ka      | 10 a                | 5.1 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl20\_b2\_paleo21            | 140 ka     | 5 a                 | 0.8 hrs            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl10\_b2\_paleo21\ :sup:`\*` | 9 ka       | 1 a                 | 1.0 hrs            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl10\_b2\_future21\_ctrl     | 100 a      | 1 a                 | 0.9 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_grl10\_b2\_future21\_asmb     | 100 a      | 1 a                 | 0.9 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_ant64\_b2\_spinup09\_init100a | 100 a      | 2 / 10 a\ :sup:`†`  | 4.1 sec            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_ant64\_b2\_spinup09\_fixtopo  | 140 ka     | 5 / 10 a\ :sup:`†`  | 0.7 hrs            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_ant64\_b2\_spinup09           | 500 a      | 2 / 10 a\ :sup:`†`  | 0.5 min            |
+-------------------------------------+------------+---------------------+--------------------+
| repo\_ant64\_b2\_future09\_ctrl     | 100 a      | 2 / 10 a\ :sup:`†`  | 6.3 sec            |
+-------------------------------------+------------+---------------------+--------------------+

| Table 1: Model times, time steps and computing (CPU) times for the EISMINT, Greenland and Antarctica test simulations contained in the script ``multi_sico_1.sh``, run with SICOPOLIS V5 (revision 03471067a) and the Intel Fortran Compiler 19.1 for Linux (optimization options ``-xHOST -O3 -no-prec-div``) on a 12-Core Intel Xeon Gold 6256 (3.6 GHz) PC under openSUSE Leap 15.5.
| \ :sup:`†`: If one value is given, this is the common dynamic (velocity, ice thickness) and thermodynamic (temperature, water content, age) time step. If two values are given (marked by the dagger (\ :sup:`†`) symbol), the first one is the dynamic, the second one the thermodynamic time step.
| \ :sup:`‡`: All runs were done on a single core only. The ``repo_ant64_b2_xxx`` runs that include ice shelves can be done on multiple cores using OpenMP for the SSA solver. However, at the employed, low resolution of 64 km the solver does not scale well, and the gain in wall clock time by using multiple cores is very small.
| \ :sup:`\*`: For this run, see the remark in the :ref:`subsection on the resolution-doubler tool <plotting_and_tools-res_dbl>`.
