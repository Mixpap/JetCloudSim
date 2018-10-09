#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            9

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  radius_cloud                   0
#define  rho_ism                        1
#define  rho_cloud                      2
#define  rho_jet                        3
#define  lorentz_jet                    4
#define  pressure_gas                   5
#define  pressure_jet                   6
#define  jet_window                     7
#define  gJET                           8

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  UNIT_DENSITY                   1.67e-24
#define  UNIT_LENGTH                    3.0e19
#define  UNIT_VELOCITY                  3.0e10

/* [End] user-defined constants (do not change this line) */
