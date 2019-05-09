#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            13
#define  INTERNAL_BOUNDARY            	YES

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  RESISTIVITY                    NO

/* -- user-defined parameters (labels) -- */

#define  radius_cloud                   0
#define  rho_ism                        1
#define  cloudx                         2
#define  cloudy                         3
#define  rho_jet                        4
#define  lorentz_jet                    5
#define  pressure_gas                   6
#define  pressure_jet_thermal           7
#define  jet_window                     8
#define  r0toR                          9
#define  gJET                           10
#define  gBphi                          11
#define  gcloud                         12

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  UNIT_DENSITY                   1.67e-24
#define  UNIT_LENGTH                    3.0e19
#define  UNIT_VELOCITY                  3.0e10

/* [End] user-defined constants (do not change this line) */
