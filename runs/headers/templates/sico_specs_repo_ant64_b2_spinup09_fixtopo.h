!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Specification file sico_specs_runname.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define RUN_SPECS_HEADER_LAST_CHANGED '2023-07-26'
!                      Date of last change

!-------- Domain --------

#define ANT
!                 Simulated domain:
!                   ANT     - Antarctica
!                   ASF     - Austfonna
!                   EISMINT - EISMINT (Phase 2 SGE and modifications)
!                   GRL     - Greenland
!                   NHEM    - Northern hemisphere
!                   SCAND   - Scandinavia
!                   TIBET   - Tibet
!                   NMARS   - North polar cap of Mars
!                   SMARS   - South polar cap of Mars
!                   XYZ     - Various domains

!-------- Physical parameter file --------

#define PHYS_PARA_FILE 'phys_para_ant_cp10_02.dat'
!                       Name of the file containing the physical parameters

!-------- Type of grid, spatial resolution --------

#define GRID 0
!                       0 : Cartesian coordinates in the stereographic plane
!                           without distortion correction
!                       1 : Cartesian coordinates in the stereographic plane
!                           with distortion correction
!                       2 : Geographical coordinates (longitude/latitude)
!                           [not allowed for this application]

#define X0 -3040.0d0
!                       x coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define Y0 -3040.0d0
!                       y coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define DX 64.0d0
!                       Horizontal grid spacing in km, for GRID==0
!                       or GRID==1
!                       [64 km requires IMAX= 95 and JMAX= 95,
!                        40 km requires IMAX=152 and JMAX=152,
!                        32 km requires IMAX=190 and JMAX=190,
!                        20 km requires IMAX=304 and JMAX=304,
!                        16 km requires IMAX=380 and JMAX=380,
!                        10 km requires IMAX=608 and JMAX=608,
!                         8 km requires IMAX=760 and JMAX=760]

#define IMAX 95
!                       IMAX+1: number of grid points in x-direction
!                               (i=0...IMAX)

#define JMAX 95
!                       JMAX+1: number of grid points in y-direction
!                               (j=0...JMAX)

#define KCMAX 80
!                       KCMAX+1: number of grid points in z-direction
!                                in cold ice (kc=0...KCMAX)

#define KTMAX 2
!                       KTMAX+1: number of grid points in z-direction
!                                in temperate ice (kt=0...KTMAX)

#define KRMAX 40
!                       KRMAX+1: number of grid points in z-direction
!                                in the bedrock (kr=0...KRMAX)

#define DEFORM 2.0d0
!                       Exponential stretch parameter of the non-equidistant
!                       grid in z-direction in cold ice
!                       (0.0d0 produces an equidistant grid)

#define CHECK_RES_IMAX_JMAX 1
!                       Compatibility check between horizontal resolution
!                       and number of grid points:
!                       0 : Not carried out
!                       1 : Carried out

!-------- Initial and final times, time steps --------

#define YEAR_ZERO 1990.0d0
!                       SICOPOLIS year zero in astronomical year numbering
!                       [ = signed year CE (AD) ]

!!!!! NOTE: All time quantities below refer to the SICOPOLIS calendar. !!!!!

#define TIME_INIT0 -140000.0d0
!                       Initial time of simulation (in a)

#define TIME_END0 0.0d0
!                       Final time of simulation (in a)

#define DTIME0 5.0d0
!                       Time step (in a) for computation of velocity
!                       and topography

#define DTIME_TEMP0 10.0d0
!                       Time step (in a) for computation of
!                       temperature, water content and age of the ice

#define DTIME_WSS0 100.0d0
!                       Time step (in a) for computation of
!                       isostatic steady-state displacement of the lithosphere
!                       (only for REBOUND==2, ELRA model)

#define DTIME_SER0 100.0d0
!                       Time step (in a) for writing of data to
!                       the time-series files

!!! #define YEAR_SEC 31556926.0d0
!                       Conversion from years to seconds;
!                       only required if supposed to be different from
!                       the default value 1 a = 31556925.445 s
!                       (IUPAC-IUGS year for epoch 2000.0)

!-------- Ice sheet dynamics --------

#define DYNAMICS 2
!                         0 : Ice flow velocity set to zero everywhere
!                             (static ice)
!                         1 : SIA for grounded ice,
!                             SSA for floating ice (if existing)
!                         2 : SIA/SStA hybrid for grounded ice,
!                             SSA for floating ice (if existing)

#define HYB_MODE 0
!                         SIA/SStA hybrid (only for DYNAMICS==2):
!                         0 : Sum of weighted sliding SIA and weighted SStA
!                             (by R. Greve)
!                         1 : Sum of weighted non-sliding SIA and full SStA
!                             (by J. Bernales)
!                         2 : Pure SStA (no SIA)

#define LIS_OPTS '-i bicgsafe -maxiter 1000 -tol 1.0e-03 -p jacobi -initx_zeros false'
!                         Options string for the Lis solver for the SSA/SStA
!                         (see the Lis User Guide, www.ssisc.org/lis)
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2)

#define TOL_ITER_SSA 0.025d0
!                         Tolerance for the nonlinear iterations of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define N_ITER_SSA 25
!                         Maximum number of nonlinear iterations of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define N_ITER_SSA_MIN 2
!                         Minimum number of nonlinear iterations of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define ITER_INIT_SSA 2
!                         Initial depth-integrated viscosity for the
!                         nonlinear iterations of the SSA/SStA:
!                         1 : Constant VISC_INIT_SSA times ice thickness
!                         2 : Previous depth-averaged viscosity times
!                             ice thickness
!                         3 : Computation by subroutine calc_vis_ssa
!                             (that is always used for the further iterations)
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define VISC_INIT_SSA 1.0d+15
!                         Initial viscosity (in Pa s) for the iterations
!                         of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define N_VISC_SMOOTH 0
!                         Number of steps for the diffusive smoothing
!                         of the depth-averaged viscosity of the SSA/SStA
!                         (if negative, abs(N_VISC_SMOOTH) smoothing steps
!                         of the logarithm of the viscosity are carried out)
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define VISC_SMOOTH_DIFF 0.0d0
!                         Dimensionless diffusivity for the diffusive smoothing
!                         of the depth-averaged viscosity of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define RELAX_FACT_SSA 0.7d0
!                         Relaxation factor for the relaxation scheme
!                         of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define RATIO_SL_THRESH 0.5d0
!                         Threshold value for the slip ratio of grounded ice.
!                         If the slip ratio is larger, hybrid SIA/SStA dynamics
!                         kicks in (for DYNAMICS==2 and HYB_MOD==0).

#define SSTA_SIA_WEIGH_FCT 2
!                         SStA-SIA weighing factor as a function of the
!                         slip ratio (for DYNAMICS==2 and HYB_MOD==0):
!                         0 : Linear function (continuous transitions)
!                         1 : Cubic function (smooth transitions)
!                         2 : Quintic function (even smoother transitions)

#define HYB_REF_SPEED 30.0d0
!                         Scaling reference speed for hybrid SIA/SStA dynamics
!                         (in m/a, for DYNAMICS==2 and HYB_MOD==1).

!-------- Ice sheet thermodynamics --------

#define CALCMOD 3
!                        -1 : ISOT: isothermal method,
!                                   constant temperature and age
!                         0 : COLD: cold-ice method, resetting of temperatures
!                                   above pressure melting
!                         1 : POLY: polythermal method,
!                                   separate domains for cold and temperate ice
!                         2 : ENTC: conventional enthalpy method
!                         3 : ENTM: melting-CTS enthalpy method
!
!                             For CALCMOD == -1, 0, 2, 3,
!                             the kt domain is redundant,
!                             therefore KTMAX==2 is recommended

#define TEMP_CONST -10.0d0
#define AGE_CONST    0.0d0
!                         Prescribed temperature and age
!                         (only for ISOT; CALCMOD==-1)

#define CTS_MELTING_FREEZING 2
!                         Treatment of the transition conditions at the CTS
!                         (only for POLY; CALCMOD==1):
!                         1 : Melting and freezing conditions distinguished
!                         2 : Always melting conditions assumed

!-------- Ice margin treatment --------

#define MARGIN 3
!                         1 : Ice extent strictly restricted to land area
!                         2 : Formation of marine ice possible
!                         3 : Formation of marine ice and ice shelves possible

#define MARINE_ICE_FORMATION 1
!                         1 : No special mechanism for the formation of marine ice
!                             (only for MARGIN==2)
!                         2 : Formation of marine ice via "underwater ice"
!                             (only for MARGIN==2)

#define MARINE_ICE_CALVING 1
!                         1 : Unlimited expansion of marine ice, no calving
!                             (only for MARGIN==2)
!                         2 : Limited expansion of marine ice,
!                             constant value z_mar=Z_MAR for the minimum elevation
!                             of the isostatically relaxed bedrock
!                             (only for MARGIN==2)
!                         3 : Limited expansion of marine ice,
!                             constant value z_mar=Z_MAR for the minimum bedrock
!                             elevation (only for MARGIN==2)
!                         4 : Limited expansion of marine ice,
!                             minimum elevation of the isostatically relaxed bedrock
!                             z_mar proportional to sea-level stand
!                             (only for MARGIN==2)
!                         5 : Limited expansion of marine ice,
!                             minimum bedrock elevation z_mar proportional to
!                             sea-level stand (only for MARGIN==2)
!                         6 : Limited expansion of marine ice,
!                             minimum elevation of the isostatically relaxed bedrock
!                             z_mar related to sea-level stand by the piecewise
!                             linear relation proposed by Zweck and Huybrechts (2005)
!                             (only for MARGIN==2)
!                         7 : Limited expansion of marine ice,
!                             minimum bedrock elevation z_mar related to sea-level
!                             stand by the piecewise linear relation proposed by
!                             Zweck and Huybrechts (2005) (only for MARGIN==2)
!                         9 : Calving of marine ice by calving law for
!                             "underwater ice"
!                             (only for MARGIN==2 and MARINE_ICE_FORMATION==2)

#define Z_MAR 0.0d0
!                         Minimum elevation (in m) of the isostatically relaxed
!                         bedrock allowed to glaciate
!                         (for MARGIN==2 and MARINE_ICE_CALVING==2)
!                         Minimum bedrock elevation (in m) allowed to glaciate
!                         (for MARGIN==2 and MARINE_ICE_CALVING==3)

#define FACT_Z_MAR 2.5d0  /* suitable value for MARGIN==2 and MARINE_ICE_CALVING==4,5 */
!                         Minimum bedrock elevation or
!                         minimum elevation of the isostatically relaxed bedrock
!                         allowed to glaciate:
!                         proportionality factor to sea-level stand
!                         (for MARGIN==2 and MARINE_ICE_CALVING==4,5),
!                         modification factor for the piecewise linear
!                         relation by Zweck and Huybrechts (2005)
!                         (for MARGIN==2 and MARINE_ICE_CALVING==6,7)

#define CALV_UW_COEFF 1.0d-04
!                         Calving coefficient for "underwater ice",
!                         in m^(1-r1-r2)*a^(-1) (for MARGIN==2,
!                         MARINE_ICE_FORMATION==2 and MARINE_ICE_CALVING==9)

#define R1_CALV_UW 1.0d0
!                         Exponent r1 of the ice thickness
!                         for the calving law for "underwater ice"
!                         (for MARGIN==2,
!                         MARINE_ICE_FORMATION==2 and MARINE_ICE_CALVING==9)

#define R2_CALV_UW 1.0d0
!                         Exponent r2 of the sea depth
!                         for the calving law for "underwater ice"
!                         (for MARGIN==2,
!                         MARINE_ICE_FORMATION==2 and MARINE_ICE_CALVING==9)

#define ICE_SHELF_CALVING 2
!                         1 : Unlimited expansion of ice shelves, no calving
!                             (only for MARGIN==3)
!                         2 : Instantaneous calving of ice shelves if the
!                             thickness is less than H_CALV (only for MARGIN==3)
!                         3 : Float-kill:
!                             Instantaneous removal of all floating ice
!                             (only for MARGIN==3)

#define H_CALV 50.0d0
!                         Threshold thickness (in m) of ice shelves for calving
!                         (for MARGIN==3 and ICE_SHELF_CALVING==2)

!-------- Flow law --------

#define FLOW_LAW 1
!                         1 : Glen's flow law with stress exponent n=3
!                         2 : Goldsby-Kohlstedt flow law with stress exponent
!                             n=1.8 and grain-size exponent p=1.4
!                         3 : Durham's flow law with stress exponent n=4
!                         4 : Smith-Morland (polynomial) flow law

#define FIN_VISC 2
!                         1 : Unmodified flow law with infinite viscosity
!                             for low strain rates
!                             (only for FLOW_LAW==1, 2, 3)
!                         2 : Modified flow law with additional
!                             finite-viscosity parameter SIGMA_RES
!                             (only for FLOW_LAW==1, 2, 3)

#define GR_SIZE 1.0d-03
!                         Average grain size (in m; only for FLOW_LAW==2)

#define SIGMA_RES 1.0d+04
!                         Residual stress (finite-viscosity contribution)
!                         in the creep response function
!                         (in Pa; only for FLOW_LAW==1, 2, 3 and FIN_VISC==2)

!-------- Flow enhancement factor --------

#define ENHMOD 1
!                         1 : Flow enhancement factor enh=ENH_FACT everywhere
!                             in grounded ice
!                         2 : enh=ENH_INTG for ice younger than AGE_TRANS_0
!                             (Holocene ice),
!                             enh=ENH_FACT for ice older than AGE_TRANS_0
!                             (Pleistocene ice);
!                             for present-day steady-state simulations
!                         3 : enh=ENH_INTG for Holocene and Eemian ice,
!                             enh=ENH_FACT for Weichselian and pre-Eemian ice
!                             (as defined by the times DATE_TRANS1_0,
!                             DATE_TRANS2_0 and DATE_TRANS3_0);
!                             for transient scenarios
!                         4 : Anisotropic flow enhancement factor (quadratic
!                             function of the shear fraction) for grounded ice,
!                             between ENH_COMPR (for compression) and
!                             ENH_SHEAR (for shear)
!                         5 : Anisotropic flow enhancement factor (quadratic
!                             function of the shear fraction) for grounded and
!                             floating ice,
!                             between ENH_COMPR (for compression) and
!                             ENH_SHEAR (for shear)

#define ENH_FACT 5.0d0
!                         Flow enhancement factor (only for ENHMOD==1, 2, 3)

#define ENH_INTG 1.0d0
!                         Separate flow enhancement factor for interglacial ice
!                         (only for ENHMOD==2, 3)

#define AGE_TRANS_0 11000.0d0
!                         Age of the Holocene/Pleistocene transition
!                         (in a; only for ENHMOD==2)

#define DATE_TRANS1_0 -132000.0d0
!                         Time of the pre-Eemian/Eemian transition
!                         (in a; only for ENHMOD==3)

#define DATE_TRANS2_0 -114500.0d0
!                         Time of the Eemian/Weichselian transition
!                         (in a; only for ENHMOD==3)

#define DATE_TRANS3_0  -11000.0d0
!                         Time of the Weichselian/Holocene transition
!                         (in a; only for ENHMOD==3)

#define ENH_COMPR 3.0d0
!                         Flow enhancement factor for compression
!                         (only for ENHMOD==4, 5)

#define ENH_SHEAR 8.0d0
!                         Flow enhancement factor for shear
!                         (only for ENHMOD==4, 5)

#define ENH_STREAM -9999.9d0
!                         Separate flow enhancement factor for ice streams
!                         (SStA dynamics;
!                          only for ENHMOD==1, 2, 3 and DYNAMICS==2;
!                          ignored if negative)

#define ENH_SHELF 1.0d0
!                         Separate flow enhancement factor for floating ice
!                         (only for ENHMOD==1, 2, 3, 4 and MARGIN==3)

!-------- Initial conditions --------

#define ANF_DAT 3
!                         1 : Present initial topography
!                         2 : Ice-free initial topography with
!                             relaxed lithosphere
!                         3 : Initial values from previous
!                             simulation

#define ZS_PRESENT_FILE   'ant_b2_64_zs.dat'
!                             Name of the file containing the present-day
!                             ice-surface topography

#define ZB_PRESENT_FILE   'ant_b2_64_zb.dat'
!                             Name of the file containing the present-day
!                             ice-base topography (only for ANF_DAT==1)

#define ZL_PRESENT_FILE   'ant_b2_64_zl.dat'
!                             Name of the file containing the present-day
!                             lithosphere-surface topography
!                             (only for ANF_DAT==1)

#define ZL0_FILE          'ant_b2_64_zl0_elra.dat'
!                             Name of the file containing the topography
!                             of the relaxed lithosphere surface

#define MASK_PRESENT_FILE 'ant_b2_64_mask.dat'
!                             Name of the file containing the present-day
!                             ice-land-ocean mask

#define MASK_REGION_FILE 'none'
!                             Name of the file containing the region mask
!                             ('none' if no file is to be defined)

#define TEMP_INIT 1
!                         Initial ice temperature conditions
!                         (only for ANF_DAT==1):
!                         1 : Constant value in the entire ice sheet
!                         2 : In each ice column equal to the
!                             local ice surface temperature
!                         3 : Ice temperature linearly increasing with depth
!                         4 : Ice-temperature profiles from analytical solution
!                             for 1-d steady-state advection-diffusion equation
!                             under the assumption of linearly varying vz
!                             [Robin (1955) solution]
!                         5 : Ice temperature from previous simulation

#define TEMP_INIT_VAL -10.0d0
!                         Prescribed initial temperature (in deg C)
!                         (only for ANF_DAT==1 and TEMP_INIT==1)

#define ANFDATNAME 'repo_ant64_b2_spinup09_init100a0002.nc'
!                             Initial-value file (only for ANF_DAT==3,
!                                  or for ANF_DAT==1 and TEMP_INIT==5)

!-------- Lithosphere (bedrock) modelling --------

#define REBOUND 2
!                         0 : No bedrock adjustment
!                         1 : Isostatic bedrock adjustment with local
!                             lithosphere and relaxing asthenosphere (LLRA model)
!                         2 : Isostatic bedrock adjustment with elastic
!                             lithosphere and relaxing asthenosphere (ELRA model)

#define FRAC_LLRA 1.0d0
!                             Fraction of isostatic compensation in the LLRA
!                             model (REBOUND==1). Range 0 <= FRAC_LLRA <= 1;
!                             0: no bedrock adjustment, 1: full adjustment.

#define TIME_LAG_MOD 1
!                         1 : Constant value for the time lag of the
!                             relaxing asthenosphere (for REBOUND==1,2)
!                         2 : Spatially varying time lag of the relaxing
!                             asthenosphere read from file (for REBOUND==1,2)

#define TIME_LAG 3000.0d0
!                         Time lag of the relaxing asthenosphere (in a)
!                         (for TIME_LAG_MOD==1)

#define TIME_LAG_FILE 'none'
!                         Name of the file containing the time lag of the
!                         relaxing asthenosphere (for TIME_LAG_MOD==2)

#define FLEX_RIG_MOD 1
!                         1 : Constant value for the flexural rigidity of the
!                             lithosphere (for REBOUND==2)
!                         2 : Spatially varying flexural rigidity of the
!                             lithosphere (for REBOUND==2)

#define FLEX_RIG 1.0d+25
!                         Flexural rigidity of the lithosphere (in Nm)
!                         (for FLEX_RIG_MOD==1)

#define FLEX_RIG_FILE 'none'
!                         Name of the file containing the flexural rigidity
!                         of the lithosphere (for FLEX_RIG_MOD==2)

!!! #define EXEC_MAKE_ZL0
!                         Special setting for generating an
!                         isostatically relaxed lithosphere surface topography
!                         by the routine make_zl0.
!                         Should be used with ANF_DAT==1
!                         (present-day initial topography).
!                         For ZL0_FILE, the present-day lithosphere surface
!                         topography (ZL_PRESENT_FILE) can be used.
!                         !!! Not to be used regularly!!!

#define Q_LITHO 1
!                         0 : No coupled heat-conducting bedrock
!                         1 : Coupled heat-conducting bedrock

!-------- Evolution of the ice thickness --------

#define THK_EVOL 3
!                         0 : No evolution of the ice thickness, kept fixed on
!                             the initial thickness
!                         1 : Evolution of the ice thickness
!                         2 : Evolution of the ice thickness, but between times
!                             TIME_TARGET_TOPO_INIT0 and TIME_TARGET_TOPO_FINAL0
!                             the ice topography (zs, zb, zl, H) is nudged
!                             towards a prescribed target with the
!                             relaxation time smoothly decreasing from
!                             TARGET_TOPO_TAU0 to zero.
!                         3 : Evolution of the ice thickness, but
!                             the ice topography (zs, zb, zl, H) is nugded
!                             towards a prescribed target with the
!                             constant relaxation time TARGET_TOPO_TAU0.
!                         4 : Evolution of the ice thickness,
!                             but maximum ice extent is constrained by the
!                             prescribed mask MASK_MAXEXTENT_FILE.

#define OCEAN_CONNECTIVITY 1
!                         0 : Ocean connectivity not enforced.
!                         1 : Ocean connectivity enforced.

#define H_ISOL_MAX 1000.0d0
!                             Maximum thickness of isolated ice points (in m)
!                             (if set to 0.0d0, isolated ice points are killed).

#define TIME_TARGET_TOPO_INIT0 0.0d0
!                             Initial time for nudging towards target topography
!                             (in a; only for THK_EVOL==2)

#define TIME_TARGET_TOPO_FINAL0 0.0d0
!                             Final time for nudging towards target topography
!                             (in a; only for THK_EVOL==2)

#define TARGET_TOPO_TAU0 100.0d0
!                             Relaxation time for
!                             nudging towards target topography
!                             (in a;
!                              only for THK_EVOL==2, 3,
!                              or for ACCSURFACE==7 and ABLSURFACE==7)

#define TARGET_TOPO_DAT_NAME 'repo_ant64_b2_spinup09_init100a0002.nc'
!                             Target-topography file
!                             (only for THK_EVOL==2, 3,
!                              or for ACCSURFACE==7 and ABLSURFACE==7)

#define MASK_MAXEXTENT_FILE 'none'
!                             Maximum ice extent mask file (only for THK_EVOL==4)

#define CALCTHK 4
!                         Solution of the ice-thickness equation:
!                         1 : Explicit scheme for the diffusive
!                             SIA ice-surface equation
!                         2 : Over-implicit scheme for the diffusive
!                             SIA ice-surface equation,
!                             iterative built-in SOR solver
!                         3 : Over-implicit scheme for the diffusive
!                             SIA ice-surface equation,
!                             iterative library-based (Lis) solver
!                         4 : Explicit scheme for the general
!                             ice-thickness equation
!                         5 : Over-implicit scheme for the general
!                             ice-thickness equation,
!                             iterative built-in SOR solver
!                         6 : Over-implicit scheme for the general
!                             ice-thickness equation,
!                             iterative library-based (Lis) solver

#define OVI_WEIGHT 1.5d0
!                       Weighing parameter for the over-implicit scheme
!                       (only for CALCTHK==2, 3, 5, 6)

#define OMEGA_SOR 1.0d0
!                       Relaxation parameter for the iterative SOR solver
!                       for systems of linear equations
!                       (0 < OMEGA_SOR < 2, only for CALCTHK==2, 5)

#define ITER_MAX_SOR 1000
!                       Maximum number of iterations for the iterative
!                       SOR solver for systems of linear equations
!                       (only for CALCTHK==2, 5)

!-------- Advection treatment in the temperature and age equations --------

#define ADV_HOR 3
!                         Discretization of horizontal advection terms in the
!                         3-d temperature and age equations:
!                         1 : Not defined
!                             (central differences would be unstable!)
!                         2 : First-order upstream using
!                             velocities on the staggered grid
!                         3 : First-order upstream using
!                             interpolated velocities on the main grid

#define ADV_VERT 3
!                         Discretization of vertical advection terms in the
!                         3-d temperature and age equations:
!                         1 : Central differences (plus artificial diffusion
!                             for the age equation)
!                         2 : First-order upstream using
!                             advection terms on the staggered grid
!                         3 : First-order upstream using
!                             interpolated advection terms on the main grid

#define AGEDIFF 5.0d-08
!                         Numerical age diffusivity
!                         (in m2/s, only for ADV_VERT==1)

!-------- Discretisation of topography gradients --------

#define TOPOGRAD 1
!                         Topography gradients at grid points with
!                         0 : second-order discretisation
!                         1 : fourth-order discretisation

#define GL_SURF_GRAD 1
!                         Surface gradients at staggered-grid points at the
!                         grounding line (between grounded and floating ice)
!                         (only for MARGIN==3):
!                         1 : central differences
!                         2 : one-sided differences into either grounded or
!                             floating ice
!                             (depending on whether the staggered-grid points
!                             are grounded or floating)

!-------- Surface temperature --------

#define TSURFACE 4
!                         1 : delta_ts = DELTA_TS0, steady state
!                         3 : Sinusoidal air-temperature forcing
!                             between delta_ts = 0 C and delta_ts =
!                             -2*SINE_AMPLIT C with period
!                             SINE_PERIOD (in a)
!                         4 : delta_ts from ice-core data
!                             (e.g., GRIP, Vostok)
!                         5 : Surface temperature interpolated by using
!                             present values, LGM anomalies and a
!                             glacial index (requires ACCSURFACE==5)
!                         6 : Surface temperature climatology and anomaly
!                             read directly from NetCDF files
!                             (requires ACCSURFACE==6, ABLSURFACE==6)

#define TEMP_PRESENT_PARA 2
!                         Parameterization of the present-day mean-annual
!                         and mean-January (summer) surface temperatures by
!                         1 : Fortuin and Oerlemans (1990)
!                             for the whole ice sheet
!                         2 : Fortuin and Oerlemans (1990),
!                             separately for three different elevation ranges
!                         (for TSURFACE<=5)

#define TEMP_PRESENT_OFFSET 0.0d0
!                         Offset for the parameterization of the present-day
!                         mean-annual and mean-July surface temperatures (in C)
!                         in order to optimise the match
!                         for the chosen reference year YEAR_ZERO
!                         (for TSURFACE<=5)

#define DELTA_TS0 0.0d0
!                       Constant air-temperature deviation for steady
!                       states (only for TSURFACE==1)

#define SINE_AMPLIT 10.0d0
!                       Amplitude (in C) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define SINE_PERIOD 100000.0d0
!                       Period (in a) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define GRIP_TEMP_FILE 'vostok_temp_searise.dat'
!                       Name of the file containing the ice-core
!                       air-temperature forcing (only for TSURFACE==4)

#define GRIP_TEMP_FACT 1.0d0
!                       Modification factor for ice-core air-temperature
!                       forcing (only for TSURFACE==4)

#define GLAC_IND_FILE 'none'
!                       Name of the file containing the glacial-index
!                       forcing (only for TSURFACE==5)

#define TEMP_MA_ANOM_FILE 'none'
!                       Name of the file containing the LGM
!                       mean-annual surface-temperature-anomaly data 
!                       (difference LGM - present; for TSURFACE==5)

#define TEMP_MJ_ANOM_FILE 'none'
!                       Name of the file containing the LGM
!                       mean-July surface-temperature-anomaly data 
!                       (difference LGM - present; for TSURFACE==5)

#define TEMP_MA_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       TEMP_MA_ANOM_FILE (for TSURFACE==5)

#define TEMP_MJ_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       TEMP_MJ_ANOM_FILE (for TSURFACE==5)

!-------- Surface precipitation --------

#define ACCSURFACE 4
!                         1 : Precipitation is constant factor ACCFACT
!                             times present distribution
!                         2 : Precipitation is coupled linearly to
!                             delta_ts, coupling parameter GAMMA_S
!                         3 : Precipitation is coupled exponentially to
!                             delta_ts, coupling parameter GAMMA_S
!                         4 : Precipitation is coupled to delta_ts by the
!                             parameterisation by Huybrechts et al. (2007)
!                             [which involves the temperature above the
!                             inversion layer]
!                         5 : Precipitation interpolated by using
!                             present values, LGM anomalies and a
!                             glacial index (requires TSURFACE==5)
!                         6 : SMB climatology and anomaly
!                             read directly from NetCDF files
!                             (requires TSURFACE==6, ABLSURFACE==6)
!                         7 : Implied SMB by Calov+ (2018, Cryosphere 12)
!                             (requires ABLSURFACE==7)

#define PRECIP_PRESENT_FILE 'ant_sr_dev1.0_64_prec_a.dat'
!                       Name of the file containing the present-day
!                       precipitation data
!                       (for ACCSURFACE<=5)

#define ACCFACT 0.0d0
!                       Constant ratio between actual and present
!                       precipitation (only for ACCSURFACE==1)

#define GAMMA_S 0.0d0
!                       Parameter in the linear or exponential relation
!                       between precipitation and delta_ts
!                       (in 1/C, only for ACCSURFACE==2, 3)

#define ELEV_DESERT 1
!                         0 : No elevation desertification
!                         1 : Elevation desertification accounted for
!                             (only for ACCSURFACE==1, 2, 3)

#define GAMMA_P     -log(2.0d0)
!                       Precipitation lapse rate for elevation desertification,
!                       in km^(-1)
!                       (only for ELEV_DESERT==1 and ACCSURFACE==1, 2, 3)

#define ZS_THRESH   1500.0d0
!                       Elevation threshold for elevation desertification, in m
!                       (only for ELEV_DESERT==1 and ACCSURFACE==1, 2, 3)

#define PRECIP_ANOM_FILE 'none'
!                       Name of the file containing the
!                       LGM precipitation-anomaly data 
!                       (ratio LGM/present; only for ACCSURFACE==5)

#define PRECIP_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       PRECIP_ANOM_FILE (for ACCSURFACE==5)

#define PRECIP_ANOM_INTERPOL 2
!                         1 : Interpolation with a linear function
!                             (for ACCSURFACE==5)
!                         2 : Interpolation with an exponential function
!                             (for ACCSURFACE==5)

#define SOLID_PRECIP 1
!                         Fraction of solid precipitation:
!                         1 : Linear function
!                             of monthly mean surface temperature
!                             by Marsiat (1994)
!                         2 : Fifth-order polynomial function
!                             of monthly mean surface temperature
!                             by Bales et al. (2009)
!                         3 : Dependency on instantaneous surface temperature
!                             (statistical approach by
!                             Huybrechts and de Wolde 1999)
!                         (for ACCSURFACE<=5)

!-------- Surface ablation --------

#define ABLSURFACE 1
!                         1 : Ablation parameterized
!                             by positive-degree-day (PDD) method.
!                             Rainfall assumed to run off instantaneously.
!                             Parameters defined in physical-parameter file.
!                         2 : Ablation parameterized
!                             by positive-degree-day (PDD) method.
!                             Rainfall assumed to contribute to formation 
!                             of superimposed ice.
!                             Parameters defined in physical-parameter file.
!                         3 : Ablation parameterized
!                             by linear-temperature-index (LTI) method.
!                         6 : SMB climatology and anomaly
!                             read directly from NetCDF files
!                             (requires TSURFACE==6, ACCSURFACE==6)
!                         7 : Implied SMB by Calov+ (2018, Cryosphere 12)
!                             (requires ACCSURFACE==7)

#define LAMBDA_LTI 500.0d0
!                       Melting coefficient for the LTI method
!                       (in (mm WE)/(a*deg C), for ABLSURFACE==3)

#define TEMP_LTI -2.0d0
!                       Threshold summer temperature for the LTI method
!                       (in deg C, for ABLSURFACE==3)

#define MB_ACCOUNT 1
!                       Mass balance accounting by "hidden ablation scheme"
!                       (by R. Calov, A. Robinson):
!                         0 : Glaciation of all inner points of the
!                             domain allowed
!                             (prevents accurate accounting of calving
!                              near the margin)
!                         1 : Outermost inner points of the domain
!                             (i=1,IMAX-1, j=1,JMAX-1) not allowed to glaciate
!                             (required for accurate accounting of calving
!                              near the margin)

!-------- Surface temperature and SMB
!         (for TSURFACE==6, ACCSURFACE==6 and ABLSURFACE==6) --------

#define TEMP_SMB_CLIMATOLOGY_FILE 'none'
!                       NetCDF file containing the
!                       surface-temperature and SMB climatology

#define TEMP_SMB_ANOM_DIR 'none'
!                       Directory for the
!                       yearly surface-temperature and SMB anomalies
!                       ('none' if no directory is to be specified)

#define TEMP_ANOM_SUBDIR 'none'
!                       Subdirectory for the
!                       yearly surface-temperature anomalies
!                       ('none' if no directory is to be specified)

#define TEMP_ANOM_FILES 'none'
!                       NetCDF files containing the
!                       yearly surface-temperature anomalies
!                       (without final year number and .nc extension)
!                       ('none' if no files are to be specified)

#define dTEMPdz_SUBDIR 'none'
!                       Subdirectory for the
!                       yearly surface-temperature vertical gradients
!                       ('none' if no directory is to be specified,
!                        'value' if constant value is to be used)

#define dTEMPdz_FILES 'none'
!                       NetCDF files containing the
!                       yearly surface-temperature vertical gradients
!                       (without final year number and .nc extension)
!                       ('none' if no files are to be specified,
!                        value as string number [in K/m]
!                        if dTEMPdz_SUBDIR is set to 'value')

#define SMB_ANOM_SUBDIR 'none'
!                       Subdirectory for the
!                       yearly SMB anomalies
!                       ('none' if no directory is to be specified)

#define SMB_ANOM_FILES 'none'
!                       NetCDF files containing the
!                       yearly SMB anomalies
!                       (without final year number and .nc extension)
!                       ('none' if no files are to be specified)

#define dSMBdz_SUBDIR 'none'
!                       Subdirectory for the
!                       yearly SMB vertical gradients
!                       ('none' if no directory is to be specified,
!                        'value' if constant value is to be used)

#define dSMBdz_FILES 'none'
!                       NetCDF files containing the
!                       yearly SMB vertical gradients
!                       (without final year number and .nc extension)
!                       ('none' if no files are to be specified,
!                        value as string number [in (m/a ice equiv.)/m]
!                        if dSMBdz_SUBDIR is set to 'value')

#define TEMP_SMB_ANOM_TIME_MIN -9999
!                       Minimum time of the yearly surface-temperature
!                       and SMB anomalies (in year CE)

#define TEMP_SMB_ANOM_TIME_MAX 9999
!                       Maximum time of the yearly surface-temperature
!                       and SMB anomalies (in year CE)

!-------- Prescribed surface mass balance correction --------

#define SMB_CORR_FILE 'none'
!                       Name of the file containing the spatially dependent
!                       correction of the surface mass balance
!                       ('none' for no correction)

!-------- Special ISMIP6 InitMIP settings for the surface mass balance --------

!!! #define INITMIP_SMB_ANOM_FILE 'smb_anomaly_64km_ISMIP6.nc'
!                       Name of the file containing the surface mass balance
!                       anomaly for ISMIP6 InitMIP

!-------- Ice-shelf collapse masks --------

#define ICE_SHELF_COLLAPSE_MASK 0
!                         0 : No ice-shelf collapse masks provided
!                         1 : Ice-shelf collapse masks provided

#define ICE_SHELF_COLLAPSE_MASK_DIR 'none'
!                       Directory for the
!                       yearly ice-shelf collapse masks
!                       (for ICE_SHELF_COLLAPSE_MASK==1)

#define ICE_SHELF_COLLAPSE_MASK_FILES 'none'
!                       NetCDF files containing the
!                       yearly ice-shelf collapse masks
!                       (without final year number and .nc extension)
!                       (for ICE_SHELF_COLLAPSE_MASK==1)

#define ICE_SHELF_COLLAPSE_MASK_H_REF_FILE 'none'
!                       NetCDF file containing the
!                       reference ice thickness for the
!                       yearly ice-shelf collapse masks
!                       (for ICE_SHELF_COLLAPSE_MASK==1)

#define ICE_SHELF_COLLAPSE_MASK_TIME_MIN -9999
!                       Minimum time of the yearly ice-shelf collapse masks
!                       (in year CE)
!                       (for ICE_SHELF_COLLAPSE_MASK==1)

#define ICE_SHELF_COLLAPSE_MASK_TIME_MAX 9999
!                       Maximum time of the yearly ice-shelf collapse masks
!                       (in year CE)
!                       (for ICE_SHELF_COLLAPSE_MASK==1)

!-------- Sea level --------

#define SEA_LEVEL 1
!                         1 : Constant sea level z_sl = Z_SL0
!                         2 : Saw-tooth-shaped sea-level forcing
!                             with glacial minimum z_sl = -130;
!                             for palaeoclimatic scenarios
!                         3 : Sea-level reconstruction from data
!                             (e.g., SPECMAP); for palaeoclimatic scenarios

#define Z_SL0 0.0d0
!                       Constant sea level
!                       (in m, only for SEA_LEVEL==1)

#define SEA_LEVEL_FILE 'none'
!                       Name of the file containing the sea-level
!                       forcing (only for SEA_LEVEL==3)

!-------- Basal hydrology --------

#define BASAL_HYDROLOGY 0
!                       0 : Water film thickness under grounded ice
!                           equal to zero everywhere
!                       1 : Water film thickness under grounded ice
!                           computed by flux routing scheme

#define MELT_DRAIN 0
!                       Input for water film under grounded ice
!                       (only for BASAL_HYDROLOGY==1):
!                       0 : Basal meltwater only
!                       1 : Basal meltwater plus surface runoff

!-------- Basal sliding --------

#define SLIDE_LAW 2
!                       1 : Weertman-type sliding,
!                           full ice pressure in denominator
!                       2 : Weertman-type sliding,
!                           reduced pressure (ice minus water) in denominator,
!                           limiter RED_PRES_LIMIT_FACT applied for SIA and SStA
!                       3 : Weertman-type sliding,
!                           reduced pressure (ice minus water) in denominator,
!                           limiter RED_PRES_LIMIT_FACT applied for SIA only

#define N_SLIDE_REGIONS 1
!                       Number of regions with different sliding laws

#define SLIDE_REGIONS_FILE 'none'
!                       File defining the regions for the sliding laws
!                       (only for N_SLIDE_REGIONS > 1)

#define C_SLIDE 11.2d0
!                       Sliding coefficient, in m/[a*Pa^(p-q)]
!                       (N_SLIDE_REGIONS separate values).
!                       Set to 0.0d0 for no-slip conditions.

#define C_SLIDE_FILTER_WIDTH 0.0d0
!                       Filtering width (spatial smoothing by Gaussian filter)
!                       for the sliding coefficient, in km.
!                       Set to 0.0d0 for no smoothing.
!                       Values > 0 only make sense
!                       for constant Weertman exponents p and q!

#define GAMMA_SLIDE 1.0d0
!                       Sub-melt sliding coefficient, in K
!                       (N_SLIDE_REGIONS separate values).
!                       Set to 1.11d+11 (or any other very large value)
!                       to allow basal sliding everywhere,
!                       irrespective of the basal temperature.

#define P_WEERT 3
!                       Weertman exponent p (integer) for the basal shear stress
!                       (N_SLIDE_REGIONS separate values)

#define Q_WEERT 2
!                       Weertman exponent q (integer) for the basal pressure
!                       (N_SLIDE_REGIONS separate values)

#define TIME_RAMP_UP_SLIDE 5000.0d0
!                       Ramp-up time (in a) for basal sliding:
!                       Sliding starts gradually between the inital time
!                       (TIME_INIT0) and the initial time plus the ramp-up time
!                       (TIME_INIT0+TIME_RAMP_UP_SLIDE).
!                       Set to 0.0d0 for immediate start of basal sliding
!                       (no gradual ramp-up).

#define RED_PRES_LIMIT_FACT 0.35d0
!                       Limiter for the reduced pressure (ice minus water),
!                       ensures that the reduced pressure cannot become smaller
!                       than RED_PRES_LIMIT_FACT times the ice pressure
!                       (for SLIDE_LAW==2,3)

#define HYDRO_SLIDE_SAT_FCT 0
!                       Saturation function for water-film-enhanced basal sliding
!                       (only for BASAL_HYDROLOGY==1):
!                       0 : Exponential function
!                           by Kleiner and Humbert (2014, J. Glaciol. 60)
!                       1 : Linear function
!                       2 : Cubic S-shape function
!                       3 : Quintic S-shape function

#define C_HW_SLIDE 0.0d0
!                       Coefficient for water-film-enhanced basal sliding
!                       (only for BASAL_HYDROLOGY==1)

#define HW0_SLIDE 1.0d-03
!                       Threshold water film thickness for water-film-enhanced
!                       basal sliding (in m, only for BASAL_HYDROLOGY==1)

!-------- Geothermal heat flux --------

#define Q_GEO_MOD 2
!                         1 : Constant geothermal heat flux defined
!                             by parameter Q_GEO
!                         2 : Spatially varying geothermal heat flux
!                             read from file

#define Q_GEO 0.0d0
!                       Constant geothermal heat flux (for Q_GEO_MOD==1),
!                       in mW/m2

#define Q_GEO_FILE 'ant_pu_64_qgeo.dat'
!                       Name of the file containing the spatially varying
!                       geothermal heat flux (for Q_GEO_MOD==2)

!-------- Basal melting at the marine ice front --------

#define MARINE_ICE_BASAL_MELTING 1
!                        Basal melting rate at the marine ice front:
!                        1 : Computed by the usual energy jump condition
!                            for grounded ice
!                        2 : Prescribed by QBM_MARINE
!                        3 : Weighed average of grounded ice melting (computed)
!                            and marine ice melting (prescribed by QBM_MARINE)

#define QBM_MARINE 0.0d0
!                        Basal melting rate at the marine ice front,
!                        in m/a water equiv. (for MARINE_ICE_BASAL_MELTING==2,3)

!-------- Basal melting for floating ice (only for MARGIN==3) --------

#define FLOATING_ICE_BASAL_MELTING 5
!                       Basal melting rate for floating ice:
!                       1 : Constant values for the continental shelf
!                           and the abyssal ocean, respectively
!                       4 : Parameterization as a function of the
!                           thermal forcing
!                           (ocean temperature minus ice shelf basal
!                           temperature)
!                       5 : Sector-wise parameterization as a function of the
!                           thermal forcing, by
!                           Greve and Galton-Fenzi (doi: 10.5281/zenodo.4766982)
!                       6 : Sector-wise parameterization as a function of the
!                           thermal forcing, by ISMIP6-Antarctica

#define QBM_FLOAT_1 2.0d0
!                       Basal melting rate for the continental shelf,
!                       in m/a water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==1)

#define QBM_FLOAT_3 10.0d0
!                       Basal melting rate for the abyssal ocean,
!                       in m/a water equiv.

#define Z_ABYSS -2500.0d0
!                       Threshold seabed elevation separating
!                       the continental shelf from the abyssal ocean, in m

#define H_W_0 75.0d0
!                       Threshold water column thickness below which
!                       basal melting is reduced, in m
!                       (0.0d0 -> no reduction)

#define TEMP_OCEAN -1.5d0
!                       Ambient ocean water temperature, in degC
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define OMEGA_QBM 10.0d0
!                       Sensitivity of basal melting to thermal forcing,
!                       in m/[a*degC^alpha] water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define ALPHA_QBM 1.0d0
!                       Exponent alpha of the thermal forcing
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define N_BM_REGIONS 16
!                       Number of basal-melting regions
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define BM_REGIONS_FILE 'ismip6_imbie2_basin_numbers_64km.dat'
!                       File defining the regions
!                       (only for N_BM_REGIONS > 1,
!                        for FLOATING_ICE_BASAL_MELTING==6)

#define GAMMA0_BM 1.4477d+04
!                       Basal melting coefficient,
!                       in m/a water equiv.
!                       (N_BM_REGIONS separate values,
!                        for FLOATING_ICE_BASAL_MELTING==6)

#define DELTA_TF_BM (/ -0.1592d0,  0.5716d0,  0.1316d0,  0.5141d0,  0.6999d0, \
                        0.2711d0,  0.0760d0, -0.1193d0,  0.4328d0,  1.0665d0, \
                       -0.0054d0, -0.6644d0, -0.0553d0, -0.1235d0, -0.0588d0, \
                        0.1001d0 /)
!                       Thermal forcing offset, in K
!                       (N_BM_REGIONS separate values,
!                        for FLOATING_ICE_BASAL_MELTING==6)

#define TF_BM_PRESENT_FILE 'ismip6_obs_thermal_forcing_1995-2017_64km_60m.nc'
!                       Name of the file containing the
!                       present-day thermal forcing data of the ocean
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define TF_BM_DIR 'none'
!                       Directory for the
!                       yearly thermal forcing data of the ocean
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define TF_BM_FILES  'none'
!                       NetCDF files containing the
!                       yearly thermal forcing data of the ocean
!                       (without final year number and .nc extension)
!                       ('none' if no such files are to be specified)
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define TF_BM_TIME_MIN -9999
!                       Minimum time of the yearly thermal forcing data
!                       of the ocean (in year CE)
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define TF_BM_TIME_MAX 9999
!                       Maximum time of the yearly thermal forcing data
!                       of the ocean (in year CE)
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define ZMIN_TF_BM 30.0d0
!                       Minimum depth (positive downward) of the
!                       thermal forcing data, in m
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define NZ_TF_BM 29
!                       NZ_TF_BM+1:
!                       number of equidistant depth points of the
!                       thermal forcing data
!                       (index count 0...NZ_TF_BM)
!                       (for FLOATING_ICE_BASAL_MELTING==6)

#define DZ_TF_BM 60.0d0
!                       Depth spacing of the thermal forcing data, in m
!                       (for FLOATING_ICE_BASAL_MELTING==6)

!  ------ Special ISMIP6 InitMIP setting

!!! #define INITMIP_BMB_ANOM_FILE 'basal_melt_anomaly_64km_ISMIP6.nc'
!                       Name of the file containing the
!                       ice-shelf basal melting anomaly for ISMIP6 InitMIP
!                       (for FLOATING_ICE_BASAL_MELTING==4,5,6)

!!! #define LARMIP_REGIONS_FILE 'LARMIP_regions_initMIPgrid_64.nc'
!                       Name of the file containing the
!                       regions for ISMIP6 LARMIP
!                       (for FLOATING_ICE_BASAL_MELTING==4,5,6)

#define LARMIP_QBM_ANOM (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
!                       Ice-shelf basal melting rate anomaly
!                       for LARMIP regions 1-5
!                       (EAIS, Ross, Amundsen, Weddell, Peninsula),
!                       in m/a ice equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==4,5,6)

!-------- Data output --------

#define NETCDF4_ENABLED 0
!                         NetCDF format for time-slice output:
!                         0 : NetCDF-3 classic
!                         1 : NetCDF-4 classic (with compression)

#define OUT_TIMES 1
!                         Output of times in all files:
!                         1 : In SICOPOLIS years
!                         2 : In astronomical year numbering
!                             [ = signed year CE (AD) ],
!                             that is, SICOPOLIS years + YEAR_ZERO

#define OUTPUT_INIT 0
!                         Output of initial conditions
!                         in time-slice files '.nc'
!                         (for prescribed output time step, OUTPUT==1,3)
!                         and in time-series files '.ser' and '.core':
!                         0 : Initial conditions are not written to
!                             output files
!                         1 : Initial conditions are written to
!                             output files

#define OUTPUT 2
!                         1 : Writing of time-slice data in files
!                             '.nc' with prescribed time step
!                         2 : Writing of time-slice data in files
!                             '.nc' at arbitrarily specified times
!                         3 : Writing of time-slice data (only 2-d fields) in
!                             files '.nc' with prescribed time step
!                             plus
!                             writing of time-slice data
!                             (full set of 2-d and 3-d fields) in files
!                             '.nc' at arbitrarily specified times

#define ERGDAT 1
!                         0 : Only 2-d fields written as time-slice data
!                             (only for OUTPUT==1,2)
!                         1 : Full set of 2-d and 3-d fields written
!                             as time-slice data (only for OUTPUT==1,2)

#define OUTPUT_FLUX_VARS 1
!                         Treatment of flux-type variables (SMB etc.) in both
!                         the time-slice and time-series outputs:
!                         1 : Snapshots
!                         2 : Averaging over the time between subsequent outputs
!                             (exception:
!                              time-slice output with all 3-d fields for
!                              OUTPUT==3, in which case snapshots are written)

#define DTIME_OUT0 0.0d0
!                             Time step (in a) for writing of
!                             time-slice data (only for OUTPUT==1,3)

#define N_OUTPUT 8
!                             Number of specified times for writing of
!                             time-slice data (only for OUTPUT==2,3,
!                             not more than 100)

#define TIME_OUT0 (/ -9000.0d0, -5000.0d0, -2000.0d0, -1000.0d0, \
                      -500.0d0,  -200.0d0,  -100.0d0,     0.0d0 /)
!                             Times (in a) for writing of time-slice
!                             data (only for OUTPUT==2,3, in increasing
!                             order from #1 to #N_OUTPUT)

!-------- Limiters etc. --------

#define NUMDIFF_H_T 0.0d0
!                       Spatial smoothing parameter for computation
!                       of H_t

#define TAU_CTS 0.0d0
!                       Numerical time lag (in a) for evolution
!                       of H_t

#define VH_MAX 5.0d+03
!                       Lower (-VH_MAX) and upper (+VH_MAX) limits of
!                       horizontal velocities vx_c/t, vy_c/t (in m/a)

#define HD_MIN   0.0d0
#define HD_MAX 500.0d0
!                       Lower and upper limits of the SIA free-surface
!                       diffusivity hdiff (in m^2/s)

#define VISC_MIN 1.0d+10
#define VISC_MAX 1.0d+25
!                       Lower and upper limits of the depth-averaged viscosity
!                       of the SS(t)A, vis_ave_g (in Pa s)

#define QBM_MIN 0.0d0
#define QBM_MAX 3.0d0
!                       Lower and upper limits of the basal melting and
!                       drainage rates Q_bm, Q_tld and Q_b_tot
!                       (in m/a water equiv.)

#define AGE_MIN 0.0d0
#define AGE_MAX 2.0d+06
!                       Lower and upper limits of computed ages (in a)

#define MEAN_ACCUM 1.0d+02
!                       Mean accumulation rate over modelled ice sheet
!                       (in mm water equiv./a)
!                       [Only required in case of CALCTHK==2, 5 for
!                       the convergence criterion of the SOR method.
!                       Need not be very precise, a rough estimate is
!                       sufficient.]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
