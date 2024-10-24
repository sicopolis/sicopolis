!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Specification file sico_specs_runname.h
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------- Basic settings --------

#define RUN_SPECS_HEADER_LAST_CHANGED '2024-10-08'
!                      Date of last change

!-------- Domain --------

#define NMARS
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

#define PARAM_RHO_I 910.0d0
!       Density of ice (in kg/m3)

#define PARAM_RHO_W 1000.0d0
!       Density of pure water (in kg/m3)

#define PARAM_RHO_SW 1028.0d0
!       Density of sea water (in kg/m3)

#define PARAM_L 3.35d+05
!       Latent heat (in J/kg)

#define PARAM_G 3.72d0
!       Gravity acceleration (in m/s2)

#define PARAM_NUE 1.0d-06
!       Water diffusivity (in kg/(m*s)) 

#define PARAM_BETA 3.3d-04
!       Clausius-Clapeyron gradient (in K/m)

#define PARAM_OMEGA_MAX 0.01d0
!       Threshold value for the water content

#define PARAM_RHO_C 2900.0d0
!       Density of crustal material (dust) (in kg/m3)

#define PARAM_KAPPA_C 2.5d0
!       Heat conductivity of crustal material (dust) (in W/(m*K))

#define PARAM_C_C 1000.0d0
!       Specific heat of crustal material (dust) (in J/(kg*K))

#define PARAM_H_R 2000.0d0
!       Thickness of the modelled lithosphere layer (in m)

#define PARAM_RHO_C_R 2.0d+06
!       Density times specific heat of the lithosphere (in J/(m3*K))

#define PARAM_KAPPA_R 3.0d0
!       Heat conductivity of the lithosphere (in W/(m*K))

#define PARAM_RHO_A 3500.0d0
!       Density of the asthenosphere (in kg/m3)

#define PARAM_R_T 181.25d0
!       Coefficient of the water-content dependence in the rate factor
!       for temperate

#define RF_KAPPA_C_FILE 'RF_KAPPA_C_GrWeHu98.nc'
!       Name of the file containing the tabulated values of the
!       temperature-dependent rate factor, heat conductivity and specific heat

!-------- Type of grid, spatial resolution --------

#define GRID 1
!                       0 : Cartesian coordinates in the stereographic plane
!                           without distortion correction
!                       1 : Cartesian coordinates in the stereographic plane
!                           with distortion correction
!                       2 : Geographical coordinates (longitude/latitude)
!                           [not allowed for this application]

#define PLANET_R 3389500.0d0
!                       Radius of Mars (in m)

#define PLANET_A 3396200.0d0
!                       Semi-major axis of Mars (in m)

#define PLANET_F_INV 169.81d0
!                       Inverse flattening of Mars

#define STEREO_PROJ_LATD0 71.0d0
!                       Standard parallel (in degN)
!                       (only for GRID==0, 1)

#define STEREO_PROJ_LOND0 0.0d0
!                       Central meridian (in degE)
!                       (only for GRID==0, 1)

#define X0 -900.0d0
!                       x coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define Y0 -900.0d0
!                       y coordinate (in km) of the origin point (i,j) = (0,0),
!                       for GRID==0 or GRID==1

#define DX 10.0d0
!                       Horizontal grid spacing in km, for GRID==0
!                       or GRID==1
!                       (20 km requires IMAX= 90 and JMAX= 90,
!                        10 km requires IMAX=180 and JMAX=180,
!                         5 km requires IMAX=360 and JMAX=360)

#define IMAX 180
!                       IMAX+1: number of grid points in x-direction
!                               (i=0...IMAX)

#define JMAX 180
!                       JMAX+1: number of grid points in y-direction
!                               (j=0...JMAX)

#define KCMAX 50
!                       KCMAX+1: number of grid points in z-direction
!                                in cold ice (kc=0...KCMAX)

#define KTMAX 10
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

#define YEAR_ZERO 2000.0d0
!                       SICOPOLIS year zero in astronomical year numbering
!                       [ = signed year CE (AD) ]

!!!!! NOTE: All time quantities below refer to the SICOPOLIS calendar. !!!!!

#define TIME_INIT0 -1.0d+07
!                       Initial time of simulation (in a)

#define TIME_END0 0.0d0
!                       Final time of simulation (in a)

#define DTIME0 1.0d+05
!                       Time step (in a) for computation of velocity
!                       and topography

#define DTIME_TEMP0 1.0d+05
!                       Time step (in a) for computation of
!                       temperature, water content and age of the ice

#define DTIME_WSS0 1.0d+05
!                       Time step (in a) for computation of
!                       isostatic steady-state displacement of the lithosphere
!                       (only for REBOUND==2, ELRA model)

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

#define CALCMOD 1
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
!                             [not allowed for this application]

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

#define FRAC_DUST 0.02d0
!                         Average volume fraction of dust in the ice.
!			  Should be within the range 0 <= FRAC_DUST <= 0.56.

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

!                         Note that if the dust content FRAC_DUST is not
!                         equal to zero, these values are multiplied by the
!                         hardening factor exp(-n*b*FRAC_DUST)
!                         [where n is the stress exponent, and b=2].

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

#define ANF_DAT 1
!                         1 : Present initial topography
!                               IMPORTANT: Make sure that the correct
!                               bed-topography file nmarsxx_zl.dat is specified
!                               below.
!                               It must fit the selected dust content
!                               FRAC_DUST, and, in case of REBOUND==0 or 1
!                               (LLRA model), the selected fraction of isostatic
!                               compensation FRAC_LLRA,
!                               or, in case of REBOUND==2 (ELRA model), the 
!                               selected flexural rigidity !!!!!!!!!!!!!!!!!!!!!
!                         2 : Ice-free initial topography with
!                             relaxed lithosphere
!                         3 : Initial values from previous
!                             simulation

#define ZS_PRESENT_FILE   'nmars10_zs.dat'
!                             Name of the file containing the present-day
!                             ice-surface topography

#define ZL_PRESENT_FILE   'nmars10_zl_flexrig3_dust002.dat'
!                             Name of the file containing the present-day
!                             lithosphere-surface topography
!                             (only for ANF_DAT==1)

#define ZL0_FILE          'nmars10_zl0_elra.dat'
!                             Name of the file containing the topography
!                             of the relaxed lithosphere surface

#define MASK_PRESENT_FILE 'nmars10_mask.dat'
!                             Name of the file containing the present-day
!                             ice-land-ocean mask

#define MASK_REGION_FILE 'none'
!                             Name of the file containing the region mask
!                             ('none' if no file is to be defined)

#define TEMP_INIT 3
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

#define TEMP_INIT_VAL -100.0d0
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
!                         1 : Isostatic bedrock adjustment with local
!                             lithosphere and relaxing asthenosphere (LLRA model)
!                         2 : Isostatic bedrock adjustment with elastic
!                             lithosphere and relaxing asthenosphere (ELRA model)

#define FRAC_LLRA 0.65d0
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

#define FLEX_RIG 3.0d+25
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

!-------- Evolution of the ice thickness --------

#define THK_EVOL 0
!                         0 : No evolution of the ice thickness, kept fixed on
!                             the initial thickness
!                         1 : Evolution of the ice thickness

#define OCEAN_CONNECTIVITY 0
!                         0 : Ocean connectivity not enforced.
!                         1 : Ocean connectivity enforced.

#define H_ISOL_MAX 1000.0d0
!                         Maximum thickness of isolated ice points (in m)
!                         (if set to 0.0d0, isolated ice points are killed).

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

#define AGEDIFF 5.0d-11
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

#define TSURFACE 6
!                         1 : delta_ts = DELTA_TS0, steady state
!                         3 : Sinusoidal air-temperature forcing
!                             between delta_ts = 0 C and delta_ts =
!                             +2*SINE_AMPLIT C with period
!                             SINE_PERIOD (in a)
!                         4 : delta_ts from north-polar insolation based
!                             on simplified two-cycle obliquity (main cycle
!                             and first modulation)
!                         5 : delta_ts from north-polar insolation based
!                             on computed obliquity and eccentricity
!                             evolutions
!                         6 : Surface temperature computed by the
!                             local insolation temperature (LIT) scheme

#define DELTA_TS0 0.0d0
!                       Constant air-temperature deviation for steady
!                       states (for TSURFACE==1);
!                       constant offset for the LIT scheme (for TSURFACE==6)

#define SINE_AMPLIT 15.0d0
!                       Amplitude (in C) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define SINE_PERIOD 1.3d+06
!                       Period (in a) for sinusoidal air-temperature
!                       forcing (only for TSURFACE==3)

#define TEMP0_MA_90N -105.0d0
!                       Present mean-annual air temperature (in C) at the
!                       north pole (only for TSURFACE==1, 2, 3)

#define C_MA 2.25d0
!                       Co-latitude factor in the parameterization of the
!                       present mean-annual air temperature
!                       (in C/deg co-lat.) (for TSURFACE==1, 2, 3, 4, 5)

#define GAMMA_MA 0.0d0
!                       Lapse rate in the parameterization of the
!                       present mean-annual air temperature (in C/m)
!                       (for TSURFACE==1, 2, 3, 4, 5)

#define INSOL_MA_90N_FILE 'mars_laskar_const.dat'
!                       Name of the file containing the mean-annual north-polar
!                       insolation (for TSURFACE==5, 6)

#define ALBEDO 0.3d0
!                       Albedo of the residual cap (only for TSURFACE==4, 5, 6)
!                       Standard value for TSURFACE==4, 5:  0.43d0
!                                      for TSURFACE==   6:  0.3d0 (planetary mean)

!-------- Surface accumulation --------

#define ACC_UNIT 2
!                         Unit of input accumulation data:
!                         1 : mm water equiv./a
!                         2 : mm (ice + dust) equiv./a

#define ACC_PRESENT 1
!                         Data for the present accumulation:
!                         1 : Spatially constant value

#define ACC_PRESENT_VAL 1.0d-01
!                         Value of present accumulation,
!                         in mm water or (ice + dust) equiv./a
!                         (only for ACC_PRESENT==1)

#define ACCSURFACE 3
!                         1 : Accumulation coupled linearly to delta_ts,
!                             coupling parameter GAMMA_S
!                         2 : Accumulation coupled exponentially to delta_ts,
!                             coupling parameter GAMMA_S
!                         3 : Accumulation coupled to delta_ts via
!                             Clausius-Clapeyron relation.

#define GAMMA_S 0.0d0
!                       Parameter in the linear or exponential relation
!                       between accumulation and delta_ts
!                       (in 1/C, only for ACCSURFACE==1, 2)

!-------- Surface ablation --------

#define ABLSURFACE 1
!                       1 : Ablation parameterized
!                           by equilibrium-line-distance approach.
!                           Equilibrium-line distance coupled linearly
!                           to delta_ts, coupling parameter GAMMA_ELD
!                       2 : Ablation parameterized
!                           by equilibrium-line-distance approach.
!                           Equilibrium-line distance coupled exponentially
!                           to delta_ts, coupling parameter GAMMA_ELD

#define ELD_0 5.5d+02
!                       Present equilibrium-line distance for the ablation
!                       parameterization (in km, only for ABLSURFACE==1, 2)

#define G_0 2.5d-04
!                       Accumulation-ablation-rate gradient,
!                       in [mm water or (ice + dust) equiv./a]/km
!                       (only for ABLSURFACE==1, 2)

#define GAMMA_ELD 0.0d0
!                       Parameter in the linear or exponential relation
!                       between equilibrium line distance (eld) and delta_ts
!                       (in 1/C, only for ABLSURFACE==1, 2)

#define GAMMA_G 1.0d0
!                       Parameter in the linear relation between
!                       accumulation-ablation-rate gradient (G_0)
!                       and accumulation ratio accum_now/accum_present
!                       (only for ABLSURFACE==1, 2)

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

!-------- Sea level --------

#define SEA_LEVEL 1
!                       1 : Constant sea level z_sl = Z_SL0

#define Z_SL0 -1.11d+11
!                       Constant sea level
!                       (in m, only for SEA_LEVEL==1)

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

#define Q_GEO 35.0d0
!                       Spatially constant GHF (in mW/m2)
!                       (only used if Q_GEO_FILE == 'none', otherwise ignored)

#define Q_GEO_FILE 'none'
!                       Name of the file containing the spatially varying GHF
!                       (set to 'none' if spatially constant GHF
!                       defined by parameter Q_GEO is to be used)

#define Q_LITHO 0
!                       0 : No coupled heat-conducting bedrock
!                           (GHF imposed directly at the grounded ice base)
!                       1 : Coupled heat-conducting bedrock
!                           (GHF imposed at the base of the
!                           thermal lithosphere layer of thickness PARAM_H_R)

!-------- Treatment of Chasma Borealis --------

#define CHASM 1
!                       1 : No special treatment of Chasma Borealis
!                       2 : Special treatment of Chasma Borealis

#define MASK_CHASM_FILE 'nmars10_mask_chasm.dat'
!                       Name of the file containing the chasm mask
!                       (for CHASM==2)

#define TIME_CHASM_INIT0 -1.11d+11
!                       Initial time for active chasm area (in a)

#define TIME_CHASM_END0 1.11d+11
!                       Final time for active chasm area (in a)

#define Q_GEO_CHASM 150.0d0
!                       Geothermal heat flux for the chasm area, in mW/m2
!                       (for CHASM==2, active between TIME_CHASM_INIT0 and
!                       TIME_CHASM_END0)

#define EROSION_CHASM 0.0d0
!                       Wind erosion rate for the chasm area,
!                       in mm water or (ice + dust) equiv./a (for CHASM==2,
!                       active between TIME_CHASM_INIT0 and TIME_CHASM_END0)

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

#define DTIME_SER0 1.0d+05
!                         Time step (in a) for writing of data to
!                         the time-series files (scalar variables)

#define DTIME_OUT0 0.0d0
!                         Time step (in a) for writing of
!                         time-slice data (only for OUTPUT==1,3)

#define N_OUTPUT 1
!                         Number of specified times for writing of
!                         time-slice data (only for OUTPUT==2,3,
!                         not more than 100)

#define TIME_OUT0 [ 0.0d0 ]
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

#define VH_MAX 10.0d0
!                       Lower (-VH_MAX) and upper (+VH_MAX) limits of
!                       horizontal velocities vx_c/t, vy_c/t (in m/a)

#define HD_MIN  0.0d0
#define HD_MAX 20.0d0
!                       Lower and upper limits of the SIA free-surface
!                       diffusivity hdiff (in m^2/s)

#define VISC_MIN 1.0d+10
#define VISC_MAX 1.0d+35
!                       Lower and upper limits of the depth-averaged viscosity
!                       of the SS(t)A, vis_ave_g (in Pa s)

#define QBM_MIN 0.0d0
#define QBM_MAX 0.3d0
!                       Lower and upper limits of the basal melting and
!                       drainage rates Q_bm, Q_tld and Q_b_tot
!                       (in m/a water equiv.)

#define AGE_MIN 0.0d0
#define AGE_MAX 1.0d+09
!                       Lower and upper limits of computed ages (in a)

#define MEAN_ACCUM 1.0d-01
!                       Mean accumulation rate over modelled ice sheet,
!                       in mm water or (ice + dust) equiv./a
!                       [Only required in case of CALCTHK==2 for
!                       the convergence criterion of the SOR method.
!                       Need not be very precise, a rough estimate is
!                       sufficient.]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
