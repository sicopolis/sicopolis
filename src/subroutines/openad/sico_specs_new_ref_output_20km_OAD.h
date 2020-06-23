!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Specification file sico_specs_runname.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define RUNNAME 'new_ref_output_20km_OAD'
!                      Name of simulation

#define SICO_VERSION '5-dev'
!                      Version number of SICOPOLIS

#define INPATH '/Will/be/inserted/automatically'
!                      Directory for SICOPOLIS input
!                      (no '/' at the end!)

#define OUTPATH '/Will/be/inserted/automatically'
!                      Directory for SICOPOLIS output
!                      (no '/' at the end!)

!-------- Domain --------

#define GRL
!                 Simulated domain:
!                   ANT      - Antarctica
!                   ASF      - Austfonna
!                   EMTP2SGE - EISMINT Phase 2 Simplified Geometry Experiment
!                   GRL      - Greenland
!                   NHEM     - Northern hemisphere
!                   SCAND    - Scandinavia
!                   TIBET    - Tibet
!                   NMARS    - North polar cap of Mars
!                   SMARS    - South polar cap of Mars
!                   XYZ      - Various domains

!-------- Physical parameter file --------

#define PHYS_PARA_FILE 'phys_para_grl_cp10_04.dat'
!                       Name of the file containing the physical parameters

!-------- Type of grid, spatial resolution --------

#define GRID 1
!                         0 : Cartesian coordinates in the stereographic plane
!                             without distortion correction
!                         1 : Cartesian coordinates in the stereographic plane
!                             with distortion correction
!                         2 : Geographical coordinates (longitude/latitude)
!                             [not allowed for this application]

#define X0 -720.0d0
!                       x coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define Y0 -3450.0d0
!                       y coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define DX 20.0d0
!                       Horizontal grid spacing in km, for GRID==0
!                       or GRID==1
!                       (40 km requires IMAX= 42 and JMAX= 72,
!                        20 km requires IMAX= 84 and JMAX=144,
!                        16 km requires IMAX=105 and JMAX=180,
!                        10 km requires IMAX=168 and JMAX=288,
!                         8 km requires IMAX=210 and JMAX=360,
!                         5 km requires IMAX=336 and JMAX=576,
!                         4 km requires IMAX=420 and JMAX=720)

#define IMAX 84
!                       IMAX+1: number of grid points in x-direction
!                               (i=0...IMAX)

#define JMAX 144
!                       JMAX+1: number of grid points in y-direction
!                               (j=0...JMAX)

#define KCMAX 24 
!                       KCMAX+1: number of grid points in z-direction
!                                in cold ice (kc=0...KCMAX)

#define KTMAX 2
!                       KTMAX+1: number of grid points in z-direction
!                                in temperate ice (kt=0...KTMAX)

#define KRMAX 4
!                       KRMAX+1: number of grid points in z-direction
!                                in the bedrock (kr=0...KRMAX)

#define DEFORM 2.0d0
!                       Exponential stretch parameter of the non-equidistant
!                       grid in z-direction in cold ice
!                       (0.0d0 produces an equidistant grid)

!-------- Initial and final times, time steps --------

#define YEAR_ZERO 1990.0d0
!                       SICOPOLIS year zero in astronomical year numbering
!                       [ = signed year CE (AD) ]

!!!!! NOTE: All time quantities below refer to the SICOPOLIS calendar. !!!!!

#define TIME_INIT0 -50.0d0
!                       Initial time of simulation (in a)

#define TIME_END0 0.0d0
!                       Final time of simulation (in a)

#define DTIME0 5.0d0
!                       Time step (in a) for computation of velocity
!                       and topography

#define DTIME_TEMP0 5.0d0
!                       Time step (in a) for computation of
!                       temperature, water content and age of the ice

#define DTIME_WSS0 100.0d0
!                       Time step (in a) for computation of
!                       isostatic steady-state displacement of the lithosphere
!                       (only for REBOUND==2, ELRA model)

#define DTIME_MAR_COA0 5.0d0
!                       Time step (in a) for computation of
!                       auxiliary fields mask_mar, cst_dist and cos_grad_tc
!                       for the ice discharge parameterization
!                       (only for DISC>0, ice discharge parameterization on)

#define DTIME_SER0 5.0d0
!                       Time step (in a) for writing of data to
!                       the time-series files

#define YEAR_SEC 31556926.0d0
!                       Conversion from years to seconds

!-------- Ice sheet dynamics --------

#define DYNAMICS 1
!                         0 : Ice flow velocity set to zero everywhere
!                             (static ice)
!                         1 : SIA for grounded ice,
!                             SSA for floating ice (if existing)
!                         2 : SIA for grounded ice,
!                             SIA/SStA hybrid for ice streams,
!                             SSA for floating ice (if existing)

#define LIS_OPTS '-i bicgsafe -maxiter 2000 -tol 1.0e-18 -p jacobi -initx_zeros false'
!                         Options string for the Lis solver for the SSA/SStA
!                         (see the Lis User Guide, www.ssisc.org/lis)
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2)

#define N_ITER_SSA 3
!                         Number of iterations for the relaxation scheme
!                         of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define RELAX_FACT_SSA 0.7d0
!                         Relaxation factor for the relaxation scheme
!                         of the SSA/SStA
!                         (for DYNAMICS==1 and MARGIN==3, or for DYNAMICS==2).

#define RATIO_SL_THRESH 0.5d0
!                         Threshold value for the slip ratio of grounded ice.
!                         If the slip ratio is larger, ice stream dynamics
!                         is employed (for DYNAMICS==2).

#define SSTA_SIA_WEIGH_FCT 1
!                         SStA-SIA weighing factor as a function of the
!                         slip ratio (for DYNAMICS==2):
!                         0 : Linear function (continuous transitions)
!                         1 : Cubic function (smooth transitions)
!                         2 : Quintic function (even smoother transitions)

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

#define MARGIN 1
!                         1 : Ice extent strictly restricted to land area
!                         2 : Formation of marine ice possible
!                         3 : Formation of marine ice and ice shelves possible

#define MARINE_ICE_FORMATION 2
!                         1 : No special mechanism for the formation of marine ice
!                             (only for MARGIN==2)
!                         2 : Formation of marine ice via "underwater ice"
!                             (only for MARGIN==2)

#define MARINE_ICE_CALVING 3
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

#define Z_MAR -300.0d0
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

#define ICE_SHELF_CALVING 1
!                         1 : Unlimited expansion of ice shelves, no calving
!                             (only for MARGIN==3)
!                         2 : Instantaneous calving of ice shelves if the
!                             thickness is less than H_CALV (only for MARGIN==3)
!                         3 : Float-kill:
!                             Instantaneous removal of all floating ice
!                             (only for MARGIN==3)

#define H_CALV 200.0d0
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

#define DATE_TRANS2_0 -114000.0d0
!                         Time of the Eemian/Weichselian transition
!                         (in a; only for ENHMOD==3)

#define DATE_TRANS3_0  -11000.0d0
!                         Time of the Weichselian/Holocene transition
!                         (in a; only for ENHMOD==3)

#define ENH_COMPR 2.0d0
!                         Flow enhancement factor for compression
!                         (only for ENHMOD==4, 5)

#define ENH_SHEAR 5.0d0
!                         Flow enhancement factor for shear
!                         (only for ENHMOD==4, 5)

#define ENH_SHELF 1.0d0
!                         Separate flow enhancement factor for floating ice
!                         (only for ENHMOD==1, 2, 3, 4 and MARGIN==3)

!-------- Initial conditions --------

#define ANF_DAT 1
!                         1 : Present initial topography
!                         2 : Ice-free initial topography with
!                             relaxed lithosphere
!                         3 : Initial values from previous
!                             simulation

#define ZS_PRESENT_FILE   'grl_b2_20_woem_zs.dat'
!                             Name of the file containing the present-day
!                             ice-surface topography

#define ZL_PRESENT_FILE   'grl_b2_20_woem_zl.dat'
!                             Name of the file containing the present-day
!                             lithosphere-surface topography
!                             (only for ANF_DAT==1)

#define ZL0_FILE          'grl_b2_20_woem_zl0_llra.dat'
!                             Name of the file containing the topography
!                             of the relaxed lithosphere surface

#define MASK_PRESENT_FILE 'grl_b2_20_woem_mask.dat'
!                             Name of the file containing the present-day
!                             ice-land-ocean mask

#define TEMP_INIT 4
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

#define ANFDATNAME 'none'
!                             Initial-value file (only for ANF_DAT==3)

#define ANFDATPATH '/Will/be/inserted/automatically/if/needed'
!                             Directory for initial-value file
!                             (no '/' at the end!)

!-------- Lithosphere (bedrock) modelling --------

#define REBOUND 0
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

#define Q_LITHO 0
!                         0 : No coupled heat-conducting bedrock
!                         1 : Coupled heat-conducting bedrock

!-------- Evolution of the ice thickness --------

#define THK_EVOL 1
!                         0 : No evolution of the ice thickness, kept fixed on
!                             the initial thickness
!                         1 : Evolution of the ice thickness
!                         2 : Evolution of the ice thickness, but between times
!                             TIME_TARGET_TOPO_INIT0 and TIME_TARGET_TOPO_FINAL0
!                             the ice topography (zs, zb, zl, H) is gradually
!                             adjusted in order to reach a prescribed target
!                             at time TIME_TARGET_TOPO_FINAL0.
!                         3 : Evolution of the ice thickness, but the
!                             ice topography (zs, zb, zl, H) is adjusted
!                             towards a prescribed target with the
!                             relaxation time TARGET_TOPO_TAU0.
!                         4 : Evolution of the ice thickness,
!                             but maximum ice extent is constrained by the
!                             prescribed mask maske_maxextent(j,i).

#define TIME_TARGET_TOPO_INIT0 0.0d0
!                             Initial time for target-topography adjustment
!                             (in a; only for THK_EVOL==2)

#define TIME_TARGET_TOPO_FINAL0 0.0d0
!                             Final time for target-topography adjustment
!                             (in a; only for THK_EVOL==2)

#define TARGET_TOPO_TAU0 100.0d0
!                             Relaxation time for the evolution of the
!                             ice topography towards a prescribed target
!                             (in a; only for THK_EVOL==3)

#define TARGET_TOPO_DAT_NAME 'none'
!                             Target-topography file (only for THK_EVOL==2, 3)

#define TARGET_TOPO_DAT_PATH '/Will/be/inserted/automatically/if/needed'
!                             Directory for target-topography file
!                             (no '/' at the end; only for THK_EVOL==2, 3)

#define MASK_MAXEXTENT_FILE 'none'
!                             Maximum ice extent mask file (only for THK_EVOL==4)

#define CALCTHK 2
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
!                       Numerical age diffusivity
!                       (in m2/s, only for ADV_VERT==1)

!-------- Discretisation of topography gradients --------

#define TOPOGRAD 1
!                             Topography gradients at grid points
!                             with
!                         0 : second-order discretisation
!                         1 : fourth-order discretisation

!-------- Surface temperature --------

#define TEMP_PRESENT_PARA 2
!                         Parameterization of the present-day mean-annual
!                         and mean-July surface temperatures by
!                         1 : Ritz et al. (1997) [no longitude dependence]
!                         2 : Fausto et al. (2009) [with longitude dependence]

#define TEMP_PRESENT_OFFSET 0.0d0
!                         Offset for the parameterization of the present-day
!                         mean-annual and mean-July surface temperatures (in C)
!                         in order to optimise the match
!                         for the chosen reference year YEAR_ZERO

#define TSURFACE 1
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
!                         6 : Surface temperature anomaly as a function of time
!                             is read directly from a NetCDF file
!                             (requires ACCSURFACE==6 and NETCDF==2)

#define DELTA_TS0 0.0d0
!                       Constant air-temperature deviation for steady
!                       states (only for TSURFACE==1)

#define SINE_AMPLIT 10.0d0
!                       Amplitude (in C) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define SINE_PERIOD 100000.0d0
!                       Period (in a) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define GRIP_TEMP_FILE 'GRIP_GICC05_105ka_GRIP-old_Delta-Ts_2.dat'
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

#define TEMP_PRECIP_ANOM_FILE 'none'
!                       Name of the NetCDF file containing the
!                       surface-temperature and precipitation anomaly data
!                       as functions of time
!                       (for TSURFACE==6 and ACCSURFACE==6)

#define TEMP_MA_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       TEMP_MA_ANOM_FILE (for TSURFACE==5)
!                       or TEMP_PRECIP_ANOM_FILE (for TSURFACE==6)

#define TEMP_MJ_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       TEMP_MJ_ANOM_FILE (for TSURFACE==5)
!                       or TEMP_PRECIP_ANOM_FILE (for TSURFACE==6)

!-------- Surface precipitation --------

#define PRECIP_PRESENT_FILE 'grl_rembo_20_precmm.dat'
!                       Name of the file containing the present-day
!                       monthly mean precipitation data
!                       ('none' if no such file is to be specified)

#define PRECIP_MA_PRESENT_FILE 'none'
!                       Name of the file containing the present-day
!                       mean annual precipitation data
!                       ('none' if no such file is to be specified)

!                       [Either PRECIP_PRESENT_FILE or PRECIP_MA_PRESENT_FILE
!                       must be specified. If both are specified,
!                       PRECIP_PRESENT_FILE will be used,
!                       while PRECIP_MA_PRESENT_FILE will be ignored.]

#define PRECIP_ZS_REF_FILE 'grl_rembo_20_zs_ref.dat'
!                       Name of the file containing the reference topography
!                       for the data in
!                       PRECIP_PRESENT_FILE or PRECIP_MA_PRESENT_FILE

#define ACCSURFACE 1
!                         1 : Precipitation is constant factor ACCFACT
!                             times present distribution
!                         2 : Precipitation is coupled linearly to
!                             delta_ts, coupling parameter GAMMA_S
!                         3 : Precipitation is coupled exponentially to
!                             delta_ts, coupling parameter GAMMA_S
!                         5 : Precipitation interpolated by using
!                             present values, LGM anomalies and a
!                             glacial index (requires TSURFACE==5)
!                         6 : Precipitation anomaly as a function of time
!                             is read directly from a NetCDF file
!                             (requires TSURFACE==6 and NETCDF==2)

#define ACCFACT 1.0d0
!                       Constant ratio between actual and present
!                       precipitation (only for ACCSURFACE==1)

#define GAMMA_S 0.070458d0 /* ACCSURFACE==3; 7.3% change per deg temp. change */
!                       Parameter in the linear or exponential relation
!                       between precipitation and delta_ts
!                       (in 1/C, only for ACCSURFACE==2, 3)

#define ELEV_DESERT 0
!                         0 : No elevation desertification
!                         1 : Elevation desertification accounted for
!                             (only for ACCSURFACE==1, 2, 3)

#define GAMMA_P     -log(2.0d0)
!                       Precipitation lapse rate for elevation desertification,
!                       in km^(-1)
!                       (only for ELEV_DESERT==1 and ACCSURFACE==1, 2, 3)

#define ZS_THRESH   2000.0d0
!                       Elevation threshold for elevation desertification, in m
!                       (only for ELEV_DESERT==1 and ACCSURFACE==1, 2, 3)

#define PRECIP_ANOM_FILE 'none'
!                       Name of the file containing the
!                       LGM precipitation-anomaly data 
!                       (ratio LGM/present; only for ACCSURFACE==5)

#define PRECIP_ANOM_FACT 1.0d0
!                       Modification factor for the anomaly data of
!                       PRECIP_ANOM_FILE (for ACCSURFACE==5)
!                       or TEMP_PRECIP_ANOM_FILE (for ACCSURFACE==6)

#define PRECIP_ANOM_INTERPOL 2
!                         1 : Interpolation with a linear function
!                             (for ACCSURFACE==5)
!                         2 : Interpolation with an exponential function
!                             (for ACCSURFACE==5)

#define SOLID_PRECIP 2
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

!-------- Surface ablation --------

#define ABLSURFACE 2
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

#define LAMBDA_LTI 500.0d0
!                       Melting coefficient for the LTI method
!                       (in (mm WE)/(a*deg C), for ABLSURFACE==3)

#define TEMP_LTI -2.0d0
!                       Threshold summer temperature for the LTI method
!                       (in deg C, for ABLSURFACE==3)

#define MB_ACCOUNT 0
!                         0 : No detailed mass balance accounting
!                         1 : Mass balance accounting by
!                             "hidden ablation scheme" (R. Calov, A. Robinson)

!-------- Prescribed surface mass balance correction --------

#define SMB_CORR_FILE 'none'
!                       Name of the file containing the spatially dependent
!                       correction of the surface mass balance
!                       ('none' for no correction)

!-------- Special ISMIP6 InitMIP settings for the surface mass balance --------

#define INITMIP_CONST_SMB 0
!                         0 or undefined : Normal computation of the SMB
!                         1 : SMB held constant at the initial distribution
!                             computed by the first call of the routine
!                             boundary

!!! #define INITMIP_SMB_ANOM_FILE 'dsmb_20_ISMIP6_v2_EPSG3413.nc'
!                       Name of the file containing the surface mass balance
!                       anomaly for ISMIP6 InitMIP

!-------- Ice discharge parameterization --------

#define DISC 1
!                         0 : Ice discharge parameterization off
!                         1 : Ice discharge parameterization on
!                         2 : Ice discharge parameterization on for
!                             interglacial, off for glacial with
!                             some transition inbetween

!!! #define EXEC_MAKE_C_DIS_0 1
!                         If defined compute c_dis_0 and stop
!                         (only for DISC>0)

#define C_DIS_0 3770.0d0
!                         Discharge parameter: scale [in m^(mD+1-mH)/s]
!                         (only for DISC>0)

#define C_DIS_FAC 1.0d0
!                         Discharge parameter: factor
!                         (only for DISC>0)

#define M_H 1.0d0
!                         Power of thickness
!                         (only for DISC>0)

#define M_D 3.0d0
!                         Power of distance
!                         (only for DISC>0)

#define R_MAR_EFF 120.0d0
!                         Width of ice marginal ring (in km)
!                         (only for DISC>0)

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

#define SEA_LEVEL_FILE 'specmap_zsl_124ka_searise.dat'
!                       Name of the file containing the sea-level
!                       forcing (only for SEA_LEVEL==3)

!-------- Basal hydrology --------

#define BASAL_HYDROLOGY 0
!                       0 : Water film thickness under grounded ice
!                           equal to zero everywhere
!                       1 : Water film thickness under grounded ice
!                           computed by flux routing scheme

!-------- Basal sliding --------

#define SLIDE_LAW 1
!                       1 : Weertman-type sliding,
!                           full ice pressure in denominator
!                       2 : Weertman-type sliding,
!                           reduced pressure (ice minus water) in denominator

#define C_SLIDE 0.0d0  /* for p=3, q=2 */
!                       Sliding coefficient, in m/[a*Pa^(p-q)]

#define SMW_COEFF 0.0d0
!                       Surface meltwater coefficient, in m^(r-s)*a^s

#define R_SMW 1
!                       Exponent r for the ice thickness in the parameterization
!                       of the surface meltwater effect (integer value)

#define S_SMW 1
!                       Exponent s for the runoff in the parameterization
!                       of the surface meltwater effect (integer value)

#define GAMMA_SLIDE 1.0d0
!                       Sub-melt sliding coefficient, in K

#define P_WEERT 3
!                       Weertman exponent p for the basal shear stress
!                       (integer value)

#define Q_WEERT 2
!                       Weertman exponent q for the basal pressure
!                       (integer value)

#define TIME_RAMP_UP_SLIDE 5000.0d0
!                       Ramp-up time (in a) for basal sliding:
!                       Sliding starts gradually between the inital time
!                       (TIME_INIT0) and the initial time plus the ramp-up time
!                       (TIME_INIT0+TIME_RAMP_UP_SLIDE).
!                       Set to 0.0d0 for immediate start of basal sliding
!                       (no gradual ramp-up).

#define RED_PRES_LIMIT_FACT 0.2d0
!                       Limiter for the reduced pressure (ice minus water),
!                       ensures that the reduced pressure cannot become smaller
!                       than RED_PRES_LIMIT_FACT times the ice pressure
!                       (for SLIDE_LAW==2)

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

#define ICE_STREAM 1
!                       1 : No special treatment of ice streams
!                       2 : Special treatment of ice streams (NEGIS)

#define MASK_SEDI_FILE 'none'
!                       Name of the file containing the
!                       hard-rock / soft-sediment (ice-stream) / ocean mask
!                       (for ICE_STREAM==2)

#define C_SLIDE_SEDI 6.72d0  /* for p=3, q=2 */
!                       Sliding coefficient, for sediment / ice stream,
!                       in m/[a*Pa^(p-q)] (for ICE_STREAM==2)

#define SMW_COEFF_SEDI 0.0d0
!                       Surface meltwater coefficient, for sediment / ice stream,
!                       in m^(r-s)*a^s (for ICE_STREAM==2)

#define R_SMW_SEDI 0
!                       Exponent r for the ice thickness in the parameterization
!                       of the surface meltwater effect, for sediment / ice stream
!                       (integer value, for ICE_STREAM==2)

#define S_SMW_SEDI 1
!                       Exponent s for the runoff in the parameterization
!                       of the surface meltwater effect, for sediment / ice stream
!                       (integer value, for ICE_STREAM==2)

#define GAMMA_SLIDE_SEDI 1.0d0
!                       Sub-melt sliding coefficient, for sediment / ice stream,
!                       in K (for ICE_STREAM==2)

#define P_WEERT_SEDI 3
!                       Weertman exponent p (integer value) for the basal
!                       shear stress, for sediment / ice stream
!                       (for ICE_STREAM==2)

#define Q_WEERT_SEDI 2
!                       Weertman exponent q (integer value) for the basal
!                       pressure, for sediment / ice stream
!                       (for ICE_STREAM==2)

!-------- Geothermal heat flux --------

#define Q_GEO_MOD 2
!                         1 : Constant geothermal heat flux defined
!                             by parameter Q_GEO
!                         2 : Spatially varying geothermal heat flux
!                             read from file

#define Q_GEO 0.0d0
!                       Constant geothermal heat flux (for Q_GEO_MOD==1),
!                       in mW/m2

#define Q_GEO_FILE 'GHF_Greenland_Ver2.0_GridEPSG3413_20km.dat'
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

#define QBM_MARINE 2.0d0
!                        Basal melting rate at the marine ice front,
!                        in m/a water equiv. (for MARINE_ICE_BASAL_MELTING==2,3)

!-------- Basal melting for floating ice (only for MARGIN==3) --------

#define FLOATING_ICE_BASAL_MELTING 4
!                        Basal melting rate for floating ice:
!                        1-3 : GRISLI parameterization.
!                              Basal melting rate at the grounding line
!                              (grounded ice point with floating ice neighbour):
!                              1: computed by the usual energy jump condition
!                                 for grounded ice
!                              2: prescribed by QBM_FLOAT_1
!                              3: weighed average of grounded ice melting
!                                 (computed) and grounding zone melting
!                                 (prescribed by QBM_FLOAT_1)
!                        4   : Parameterization as a function of the thermal
!                              forcing (ocean temperature minus ice shelf basal
!                              temperature)

#define QBM_FLOAT_1  1.0d0
!                       Basal melting rate in the grounding line zone,
!                       in m/a water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==1,2,3)

#define QBM_FLOAT_2  1.0d0
!                       Basal melting rate over the continental shelf,
!                       in m/a water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==1,2,3)

#define QBM_FLOAT_3 10.0d0
!                       Basal melting rate over the abyssal ocean,
!                       in m/a water equiv.

#define Z_ABYSS -2500.0d0
!                       Threshold seabed elevation separating
!                       the continental shelf from the abyssal ocean, in m

#define TEMP_OCEAN -1.8d0
!                       Ambient ocean water temperature, in degC
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define OMEGA_QBM 5.0d0
!                       Sensitivity of basal melting to thermal forcing,
!                       in m/[a*degC^alpha] water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define ALPHA_QBM 1.0d0
!                       Exponent alpha of the thermal forcing
!                       (for FLOATING_ICE_BASAL_MELTING==4)

#define H_W_0 75.0d0
!                       Threshold water column thickness below which
!                       basal melting is reduced (0.0d0 -> no reduction)
!                       (for FLOATING_ICE_BASAL_MELTING==4)

!-------- Data output --------

#define OUT_TIMES 1
!                         Output of times in all files:
!                         1 : In SICOPOLIS years
!                         2 : In astronomical year numbering
!                             [ = signed year CE (AD) ],
!                             that is, SICOPOLIS years + YEAR_ZERO

#define OUTPUT_INIT 0
!                         Output of initial conditions
!                         in time-slice files '.erg' or '.nc'
!                         (for prescribed output time step, OUTPUT==1,3)
!                         and in time-series files '.ser' and '.core':
!                         0 : Initial conditions are not written to
!                             output files
!                         1 : Initial conditions are written to
!                             output files

#define OUTPUT 2
!                         1 : Writing of time-slice data in files
!                             '.erg' or '.nc' with prescribed time step
!                         2 : Writing of time-slice data in files
!                             '.erg' or '.nc' at arbitrarily specified times
!                         3 : Writing of time-slice data (only 2-d fields) in
!                             files '.erg' or '.nc' with prescribed time step
!                             plus
!                             writing of time-slice data
!                             (full set of 2-d and 3-d fields) in files
!                             '.erg' or '.nc' at arbitrarily specified times

#define ERGDAT 1
!                         0 : Only 2-d fields written as time-slice data
!                             (only for OUTPUT==1,2)
!                         1 : Full set of 2-d and 3-d fields written
!                             as time-slice data (only for OUTPUT==1,2)

#define NETCDF 1
!                         1 : Time-slice files unformatted ('.erg')
!                         2 : Time-slice files in netCDF format ('.nc'),
!                             and writing of an extended time-series file
!                             '_ser.nc' in netCDF format

#define DTIME_OUT0 0.0d0
!                             Time step (in a) for writing of
!                             time-slice data (only for OUTPUT==1,3)

#define N_OUTPUT 7
!                             Number of specified times for writing of
!                             time-slice data (only for OUTPUT==2,3,
!                             not more than 20)

#define TIME_OUT0_01 -114000.0d0
#define TIME_OUT0_02  -65000.0d0
#define TIME_OUT0_03  -21000.0d0
#define TIME_OUT0_04   -9000.0d0
#define TIME_OUT0_05   -4000.0d0
#define TIME_OUT0_06    -100.0d0
#define TIME_OUT0_07       0.0d0
#define TIME_OUT0_08 1.11d+11
#define TIME_OUT0_09 1.11d+11
#define TIME_OUT0_10 1.11d+11
#define TIME_OUT0_11 1.11d+11
#define TIME_OUT0_12 1.11d+11
#define TIME_OUT0_13 1.11d+11
#define TIME_OUT0_14 1.11d+11
#define TIME_OUT0_15 1.11d+11
#define TIME_OUT0_16 1.11d+11
#define TIME_OUT0_17 1.11d+11
#define TIME_OUT0_18 1.11d+11
#define TIME_OUT0_19 1.11d+11
#define TIME_OUT0_20 1.11d+11
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

#define VH_MAX 1.0d+04
!                       Lower (-VH_MAX) and upper (+VH_MAX) limits of
!                       horizontal velocities vx_c/t, vy_c/t (in m/a)

#define HD_MIN   0.0d0
#define HD_MAX 500.0d0
!                       Lower and upper limits of the SIA free-surface
!                       diffusivity hdiff (in m^2/s)

#define QBM_MIN 0.0d0
#define QBM_MAX 3.0d0
!                       Lower and upper limits of the basal melting and
!                       drainage rates Q_bm, Q_tld and Q_b_tot
!                       (in m/a water equiv.)

#define AGE_MIN 0.0d0
#define AGE_MAX 1.0d+06
!                       Lower and upper limits of computed ages (in a)

#define MEAN_ACCUM 3.0d-02
!                       Mean accumulation rate over modelled ice sheet
!                       (in mm water equiv./a)
!                       [Only required in case of CALCTHK==2, 5 for
!                       the convergence criterion of the SOR method.
!                       Need not be very precise, a rough estimate is
!                       sufficient.]

#define ALLOW_COST
#define AGE_COST
#define ALLOW_OPENAD
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
