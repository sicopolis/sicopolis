!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Specification file sico_specs_runname.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define RUN_SPECS_HEADER_LAST_CHANGED '2024-10-08'
!                      Date of last change

!-------- Domain --------

#define EISMINT
!                 Simulated domain:
!                   ANT     - Antarctica
!                   GRL     - Greenland
!                   NHEM    - Entire northern hemisphere
!                   LCIS    - Laurentide and Cordilleran ice sheets
!                   SCAND   - Fennoscandian and Eurasian ice sheets
!                   TIBET   - Tibetan ice sheet
!                   ASF     - Austfonna
!                   NPI     - North Patagonian ice field
!                   MOCHO   - Mocho-Choshuenco ice cap
!                   EISMINT - EISMINT (Phase 2 SGE and modifications)
!                   HEINO   - ISMIP HEINO
!                   NMARS   - North polar cap of Mars
!                   SMARS   - South polar cap of Mars
!                   XYZ     - Unspecified domain

#define XYZ_SPECIAL_MODULES 0
!                   Only for unspecified domain XYZ:
!                       0 : Common modules will be used
!                       1 : Special modules 'boundary_m', 'sico_init_m'
!                           and 'sico_vars_m' (in 'src/subroutines/xyz/')
!                           will be used

!-------- Physical parameters --------

#define PARAM_RHO 910.0d0
!       Density of ice (in kg/m3)

#define PARAM_RHO_W 1000.0d0
!       Density of pure water (in kg/m3)

#define PARAM_RHO_SW 1028.0d0
!       Density of sea water (in kg/m3)

#define PARAM_L 3.35d+05
!       Latent heat (in J/kg)

#define PARAM_G 9.81d0
!       Gravity acceleration (in m/s2)

#define PARAM_NUE 1.0d-06
!       Water diffusivity (in kg/(m*s)) 

#define PARAM_BETA 8.7d-04
!       Clausius-Clapeyron gradient (in K/m)

#define PARAM_DELTA_TM_SW 1.85d0
!       Melting point depression of sea water
!       due to its average salinity (in degC)

#define PARAM_OMEGA_MAX 0.01d0
!       Threshold value for the water content

#define PARAM_H_R 2000.0d0
!       Thickness of the modelled lithosphere layer (in m)

#define PARAM_RHO_C_R 2.0d+06
!       Density times specific heat of the lithosphere (in J/(m3*K))

#define PARAM_KAPPA_R 3.0d0
!       Heat conductivity of the lithosphere (in W/(m*K))

#define PARAM_RHO_A 3300.0d0
!       Density of the asthenosphere (in kg/m3)

#define PARAM_R_T 0.0d0
!       Coefficient of the water-content dependence in the rate factor
!       for temperate

#define RF_KAPPA_C_FILE 'RF_KAPPA_C_vialov.dat'
!       Name of the file containing the tabulated values of the
!       temperature-dependent rate factor, heat conductivity and specific heat

!-------- Type of grid, spatial resolution --------

#define GRID 0
!                       0 : Cartesian coordinates in the stereographic plane
!                           without distortion correction
!                       1 : Cartesian coordinates in the stereographic plane
!                           with distortion correction
!                       2 : Geographical coordinates (longitude/latitude)
!                           [not allowed for this application]

#define PLANET_R 6371000.0d0
!                       Radius of the Earth (in m)

#define PLANET_A 6378137.0d0
!                       Semi-major axis of the Earth (in m)

#define PLANET_F_INV 298.257223563d0
!                       Inverse flattening of the Earth

#define STEREO_PROJ_LATD0 71.0d0
!                       Standard parallel (in degN)
!                       (only for GRID==0, 1)

#define STEREO_PROJ_LOND0 0.0d0
!                       Central meridian (in degE)
!                       (only for GRID==0, 1)

#define X0 -750.0d0
!                       x coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define Y0 -750.0d0
!                       y coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define DX 25.0d0
!                       Horizontal grid spacing in km, for GRID==0 or GRID==1
!                           (  5 km    requires IMAX=300 and JMAX=300,
!                              6.25 km requires IMAX=240 and JMAX=240,
!                             10 km    requires IMAX=150 and JMAX=150,
!                             12.5 km  requires IMAX=120 and JMAX=120,
!                             20 km    requires IMAX= 75 and JMAX= 75,
!                             25 km    requires IMAX= 60 and JMAX= 60,
!                             50 km    requires IMAX= 30 and JMAX= 30,
!                             75 km    requires IMAX= 20 and JMAX= 20,
!                            100 km    requires IMAX= 15 and JMAX= 15,
!                            125 km    requires IMAX= 12 and JMAX= 12,
!                            150 km    requires IMAX= 10 and JMAX= 10,
!                            250 km    requires IMAX=  6 and JMAX=  6)

#define IMAX 60
!                       IMAX+1: number of grid points in x-direction
!                               (i=0...IMAX)

#define JMAX 60
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

#define YEAR_ZERO 0.0d0
!                       SICOPOLIS year zero in astronomical year numbering
!                       [ = signed year CE (AD) ]

!!!!! NOTE: All time quantities below refer to the SICOPOLIS calendar. !!!!!

#define TIME_INIT0 0.0d0
!                       Initial time of simulation (in a)

#define TIME_END0 100000.0d0
!                       Final time of simulation (in a)

#define DTIME0 20.0d0
!                       Time step (in a) for computation of velocity
!                       and topography

#define DTIME_TEMP0 20.0d0
!                       Time step (in a) for computation of
!                       temperature, water content and age of the ice

!!! #define YEAR_SEC 31556926.0d0
!                       Conversion from years to seconds;
!                       only required if supposed to be different from
!                       the default value 1 a = 31556925.445 s
!                       (IUPAC-IUGS year for epoch 2000.0)

!-------- Ice sheet dynamics --------

#define DYNAMICS 1
!                         0 : Ice flow velocity set to zero everywhere
!                             (static ice)
!                         1 : SIA for grounded ice,
!                             SSA for floating ice (if existing)
!                         2 : SIA/SStA hybrid for grounded ice,
!                             SSA for floating ice (if existing)
!                         3 : DIVA for grounded ice,
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

#define SSTA_SIA_WEIGH_FCT 1
!                         SStA-SIA weighting factor as a function of the
!                         slip ratio (for DYNAMICS==2 and HYB_MOD==0):
!                         0 : Linear function (continuous transitions)
!                         1 : Cubic function (smooth transitions)
!                         2 : Quintic function (even smoother transitions)

#define HYB_REF_SPEED 30.0d0
!                         Scaling reference speed for hybrid SIA/SStA dynamics
!                         (in m/a, for DYNAMICS==2 and HYB_MOD==1).

!-------- Ice sheet thermodynamics --------

#define CALCMOD -1
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

#define FIN_VISC 1
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

#define ENH_FACT 1.0d0
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

#define ANF_DAT 2
!                         1 : Present initial topography
!                         2 : Ice-free initial topography with
!                             relaxed lithosphere
!                         3 : Initial values from previous
!                             simulation

#define ZS_PRESENT_FILE   'emtp2sge025_topo.nc'
!                             Name of the file containing the present-day
!                             ice-surface topography

#define ZB_PRESENT_FILE   'emtp2sge025_topo.nc'
!                             Name of the file containing the present-day
!                             ice-base topography (only for ANF_DAT==1)

#define ZL_PRESENT_FILE   'emtp2sge025_topo.nc'
!                             Name of the file containing the present-day
!                             lithosphere-surface topography
!                             (only for ANF_DAT==1)

#define ZL0_FILE          'emtp2sge025_topo.nc'
!                             Name of the file containing the topography
!                             of the relaxed lithosphere surface

#define MASK_PRESENT_FILE 'emtp2sge025_topo.nc'
!                             Name of the file containing the present-day
!                             ice-land-ocean mask

#define MASK_REGION_FILE 'none'
!                             Name of the file containing the region mask
!                             ('none' if no file is to be defined)

#define TEMP_INIT 2
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
!                         Prescribed initial temperature (in degC)
!                         (only for ANF_DAT==1 and TEMP_INIT==1)

#define ANFDATNAME 'none'
!                         Initial-value file (only for ANF_DAT==3,
!                         or for ANF_DAT==1 and TEMP_INIT==5)

#define RESTART 1
!                         Restart method (only for ANF_DAT==3):
!                         0 : Legacy restart
!                         1 : Improved restart (smoother transition)

!-------- Lithosphere (bedrock) modelling --------

#define REBOUND 0
!                         0 : No bedrock adjustment
!                         1 : Isostatic bedrock adjustment with
!                             time lag (LLRA model)

#define FRAC_LLRA 1.0d0
!                             Fraction of isostatic compensation in the LLRA
!                             model (REBOUND==1). Range 0 <= FRAC_LLRA <= 1;
!                             0: no bedrock adjustment, 1: full adjustment.

#define TIME_LAG_MOD 1
!                         1 : Constant value for the time lag of the
!                             relaxing asthenosphere (for REBOUND==1)
!                         2 : Spatially varying time lag of the relaxing
!                             asthenosphere read from file (for REBOUND==1)

#define TIME_LAG 3000.0d0
!                         Time lag of the relaxing asthenosphere (in a)
!                         (for TIME_LAG_MOD==1)

#define TIME_LAG_FILE 'none'
!                         Name of the file containing the time lag of the
!                         relaxing asthenosphere (for TIME_LAG_MOD==2)

!-------- Evolution of the ice thickness --------

#define THK_EVOL 1
!                         0 : No evolution of the ice thickness, kept fixed on
!                             the initial thickness
!                         1 : Evolution of the ice thickness
!                         2 : Evolution of the ice thickness, but
!                             the ice topography (zs, zb, zl, H) is nugded
!                             towards a prescribed target with a
!                             time-dependent relaxation time
!                             read from the file TARGET_TOPO_TAU0_FILE.
!                         3 : Evolution of the ice thickness, but
!                             the ice topography (zs, zb, zl, H) is nugded
!                             towards a prescribed target with the
!                             constant relaxation time TARGET_TOPO_TAU0.

#define OCEAN_CONNECTIVITY 1
!                         0 : Ocean connectivity not enforced.
!                         1 : Ocean connectivity enforced.

#define H_ISOL_MAX 1000.0d0
!                         Maximum thickness of isolated ice points (in m)
!                         (if set to 0.0d0, isolated ice points are killed).

#define TARGET_TOPO_TAU0_FILE 'none'
!                         Name of the file containing the time-dependent
!                         relaxation time for
!                         nudging towards target topography
!                         (only for THK_EVOL==2)

#define TARGET_TOPO_TAU0 100.0d0
!                         Relaxation time for
!                         nudging towards target topography
!                         (in a;
!                          only for THK_EVOL==3)

#define TARGET_TOPO_DAT_NAME 'none'
!                         Target-topography file
!                         (only for THK_EVOL==2, 3)

#define MASK_MAXEXTENT_FILE 'none'
!                         Maximum ice extent mask file (only for THK_EVOL>=1)
!                         ('none' if no file is to be defined)

#define CALCTHK 2
!                         Solution of the ice-thickness equation:
!                         1 : Explicit scheme for the diffusive
!                             SIA ice-surface equation
!                         2 : Over-implicit scheme for the diffusive
!                             SIA ice-surface equation
!                             (iterative built-in SOR solver)
!                         4 : Explicit scheme for the general
!                             ice-thickness equation

#define OVI_WEIGHT 1.5d0
!                         Weighting parameter for the over-implicit scheme
!                         (only for CALCTHK==2)

#define OMEGA_SOR 1.0d0
!                         Relaxation parameter for the iterative SOR solver
!                         for systems of linear equations
!                         (0 < OMEGA_SOR < 2, only for CALCTHK==2)

#define ITER_MAX_SOR 1000
!                         Maximum number of iterations for the iterative
!                         SOR solver for systems of linear equations
!                         (only for CALCTHK==2)

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

!-------- Parameters for surface temperature and surface mass balance --------

#define SURFACE_FORCING 1
!                         Surface temperature and mass balance forcings:
!                         1 : Distance-dependent forcings
!                             (original EISMINT set-up by Payne et al. (2000))
!                         2 : Elevation-dependent forcings
!                             (Greve and Blatter 2009, Sect. 5.7.4)

#define TEMP_MIN -1.0d+01   /* in degC     , for SURFACE_FORCING==1 */
#define S_T       0.0d+00   /* in K/km     , for SURFACE_FORCING==1 */
#define B_MAX     3.0d-01   /* in m/a      , for SURFACE_FORCING==1 */
#define S_B       1.0d-02   /* in m/(a*km) , for SURFACE_FORCING==1 */
#define ELD       1.1d+11   /* in km       , for SURFACE_FORCING==1 */

#define TEMP_0   -3.5d0     /* in degC     , for SURFACE_FORCING==2 */
#define GAMMA_T   8.0d0     /* in K/km     , for SURFACE_FORCING==2 */
#define S_0       5.0d-01   /* in m/a      , for SURFACE_FORCING==2 */
#define M_0       5.0d0     /* in m/(a*km) , for SURFACE_FORCING==2 */
#define ELA       1.1d0     /* in km       , for SURFACE_FORCING==2 */

#define MB_ACCOUNT 0
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

!-------- Surface temperature --------

#define TSURFACE 1
!                         1 : delta_ts = DELTA_TS0, steady state
!                         3 : Sinusoidal air-temperature forcing
!                             between delta_ts = 0 C and delta_ts =
!                             -2*SINE_AMPLIT C with period
!                             SINE_PERIOD (in a)
!                         4 : delta_ts from ice-core data
!                             (e.g., GRIP, Vostok)
!                         5 : Arbitrary scenario to be added in
!                             subroutine 'boundary'

#define DELTA_TS0 0.0d0
!                       Constant air-temperature deviation for steady
!                       states (only for TSURFACE==1)

#define SINE_AMPLIT 10.0d0
!                       Amplitude (in C) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define SINE_PERIOD 100000.0d0
!                       Period (in a) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define GRIP_TEMP_FILE 'grip95_temp_nosmooth.dat'
!                       Name of the file containing the ice-core
!                       air-temperature forcing (only for TSURFACE==4)

#define GRIP_TEMP_FACT 1.0d0
!                       Modification factor for ice-core air-temperature
!                       forcing (only for TSURFACE==4)

!-------- Sea level --------

#define SEA_LEVEL 1
!                       1 : Constant sea level z_sl = Z_SL0
!                       3 : Time-dependent sea level read from file

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

#define SLIDE_LAW 0
!                       0 : No-slip
!                       1 : Weertman-Budd sliding law

#define N_SLIDE_REGIONS 1
!                       Number of regions with different sliding laws

#define SLIDE_REGIONS_FILE 'none'
!                       File defining the regions for the sliding laws
!                       (only for N_SLIDE_REGIONS > 1)

#define BASAL_WATER_PRESSURE 0
!                       Basal water pressure:
!                       0 : Zero everywhere
!                       1 : Ocean pressure without cut-off (can become negative)
!                       2 : Ocean pressure with cut-off

#define C_SLIDE 0.0d0
!                       Sliding coefficient, in m/[a*Pa^(p-q)]
!                       (N_SLIDE_REGIONS separate values).

#define C_SLIDE_FILTER_WIDTH 0.0d0
!                       Filtering width (spatial smoothing by Gaussian filter)
!                       for the sliding coefficient, in km.
!                       Set to 0.0d0 for no smoothing.
!                       Values > 0 only make sense
!                       for constant Weertman exponents p and q!

#define GAMMA_SLIDE 0.0d0
!                       Sub-melt sliding coefficient, in K
!                       (N_SLIDE_REGIONS separate values).
!                       Set to 1.11d+11 (or any other very large value)
!                       to allow basal sliding everywhere,
!                       irrespective of the basal temperature.

#define P_WEERT 1
!                       Weertman exponent p (integer) for the basal shear stress
!                       (N_SLIDE_REGIONS separate values)

#define Q_WEERT 0
!                       Weertman exponent q (integer) for the basal pressure
!                       (N_SLIDE_REGIONS separate values)

#define TIME_RAMP_UP_SLIDE 0.0d0
!                       Ramp-up time (in a) for basal sliding:
!                       Sliding starts gradually between the inital time
!                       (TIME_INIT0) and the initial time plus the ramp-up time
!                       (TIME_INIT0+TIME_RAMP_UP_SLIDE).
!                       Set to 0.0d0 for immediate start of basal sliding
!                       (no gradual ramp-up).

#define RED_PRES_LIMIT_FACT 0.35d0
!                       Limiter for the reduced pressure (ice minus water);
!                       ensures that, for SIA dynamics,
!                       the reduced pressure cannot become smaller
!                       than RED_PRES_LIMIT_FACT times the ice pressure

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

!-------- Geothermal heat flux (GHF) --------

#define Q_GEO 42.0d0
!                       Spatially constant GHF (in mW/m2)

#define Q_LITHO 0
!                       0 : No coupled heat-conducting bedrock
!                           (GHF imposed directly at the grounded ice base)
!                       1 : Coupled heat-conducting bedrock
!                           (GHF imposed at the base of the
!                           thermal lithosphere layer of thickness PARAM_H_R)

!-------- Basal melting at the marine ice front --------

#define MARINE_ICE_BASAL_MELTING 1
!                       Basal melting rate at the marine ice front:
!                       1 : Computed by the usual energy jump condition
!                           for grounded ice
!                       2 : Prescribed by QBM_MARINE
!                       3 : Weighed average of grounded ice melting (computed)
!                           and marine ice melting (prescribed by QBM_MARINE)

#define QBM_MARINE 0.0d0
!                       Basal melting rate at the marine ice front,
!                       in m/a water equiv. (for MARINE_ICE_BASAL_MELTING==2,3)

!-------- Basal melting for floating ice (only for MARGIN==3) --------

#define FLOATING_ICE_BASAL_MELTING 1
!                       Basal melting rate for floating ice:
!                       1 : Constant values for the continental shelf
!                           and the abyssal ocean, respectively
!                       4 : Parameterization as a function of the
!                           thermal forcing
!                           (ocean temperature minus ice shelf basal
!                           temperature)

#define QBM_FLOAT_1 0.0d0
!                       Basal melting rate for the continental shelf,
!                       in m/a water equiv.
!                       (for FLOATING_ICE_BASAL_MELTING==1)

#define QBM_FLOAT_3 10.0d0
!                       Basal melting rate for the abyssal ocean,
!                       in m/a water equiv.

#define Z_ABYSS -1.11d+11
!                       Threshold seabed elevation separating
!                       the continental shelf from the abyssal ocean, in m

#define H_W_0 0.0d0
!                       Threshold water column thickness below which
!                       basal melting is reduced, in m
!                       (0.0d0 -> no reduction)

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

#define DTIME_SER0 100.0d0
!                         Time step (in a) for writing of data to
!                         the time-series files (scalar variables)

#define DTIME_OUT0 0.0d0
!                         Time step (in a) for writing of
!                         time-slice data (only for OUTPUT==1,3)

#define N_OUTPUT 3
!                         Number of specified times for writing of
!                         time-slice data (only for OUTPUT==2,3,
!                         not more than 100)

#define TIME_OUT0 [ 5000.0d0, 10000.0d0, 100000.0d0 ]
!                         Times (in a) for writing of time-slice
!                         data (only for OUTPUT==2,3, in increasing
!                         order from #1 to #N_OUTPUT)

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
!                       Lower and upper limits of free-surface
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
#define AGE_MAX 1.0d+06
!                       Lower and upper limits of computed ages (in a)

#define MEAN_ACCUM 5.0d+02
!                       Mean accumulation rate over modelled ice sheet
!                       (in mm water equiv./a)
!                       [Only required in case of CALCTHK==2 for
!                       the convergence criterion of the SOR method.
!                       Need not be very precise, a rough estimate is
!                       sufficient.]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
