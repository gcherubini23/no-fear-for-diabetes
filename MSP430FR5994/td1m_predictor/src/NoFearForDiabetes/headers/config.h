#ifndef CONFIG_H
#define CONFIG_H

// ---------------- DIMENSIONS ----------------
#define STATE_SPACE         13      // ['Qsto1','Qsto2','Qgut','Gp','Gt','Gsc','Il','Ip','I1','Id','X','Isc1','Isc2']
#define EXTRA_STATE_SPACE   6       // ['insulin_to_infuse','last_IIR','CHO_to_eat','D','lastQsto','is_eating']
#define INPUT_SPACE         2       // ['CHO', 'IIR']
#define MODEL_INPUT_SPACE   2       // ['CHO_consumed_rate','IIR_dt']
#define NUM_PARAMS          46
#define NUM_TUNING_PARAMS   31
// ------------------- END --------------------

// ---------------- INDICES -------------------
// State x:
#define Q_STO1      0
#define Q_STO2      1
#define Q_GUT       2
#define G_P         3
#define G_T         4
#define G_SC        5
#define I_L         6
#define I_P         7
#define I_1         8
#define I_D         9
#define X           10
#define I_SC1       11
#define I_SC2       12
// State y:
#define INS_TO_INF  0
#define LAST_IIR    1
#define CHO_TO_EAT  2
#define D           3
#define LAST_Q_STO  4
#define IS_EATING   5
// Input u:
#define CHO         0
#define IIR         1
// Input v:
#define CHO_DT      0
#define IIR_DT      1
// Parameters:
#define BW          0
#define EAT_RATE    1
#define BASAL       2
#define U2SS        3
#define GB          4
//  Sensor features
#define TD          5
//  Glucose Kinetics
#define VG          6
#define K1          7
#define K2          8
#define GPB         9
//  Insulin Kinetics
#define VI          10
#define HEB         11
#define CL          12
#define M1          13
#define M2          14
#define M30         15
#define M4          16
#define M5          17
#define IPB         18
#define ILB         19
#define IB          20
//  Rate of Appearance
#define KMAX        21
#define KMIN        22
#define KABS        23
#define KGRI        24
#define F           25
#define BETA        26
#define DELTA       27
//  Endogenous Glucose Production
#define KP1         28
#define KP2         29
#define KP3         30
#define KI          31
#define EGPB        32
//  Utilization
#define FCNS        33
#define GTB         34
#define KM0         35
#define VM0         36
#define VMX         37
#define P2U         38
//  Insulin Infusion
#define KD          39
#define KA1         40
#define KA2         41
#define ISC1SS      42
#define ISC2SS      43
//  Renal Excretion
#define KE1         44
#define KE2         45
// ------------------- END -----------------------

// ---------------- PARAMETERS -------------------
#define TRUE        1
#define FALSE       0
#define EPS         0.00001
#define UPDATE_RATE 20
#define SEC_IN_MIN  1.0f
// Kalman Filter:
#define DT          1
#define CGM_MARD    15
#define MODEL_COV   1000
#define SENSOR_COV  1000
// ------------------- END -----------------------

// -------------- FRAM VARIABLES -----------------
extern float params[NUM_PARAMS];
extern float mA[STATE_SPACE][STATE_SPACE];
extern float mB[STATE_SPACE][MODEL_INPUT_SPACE];
extern float mD[STATE_SPACE];
extern float mQ[STATE_SPACE][STATE_SPACE];
extern float mR;
extern float mH[STATE_SPACE];
extern float x[STATE_SPACE];
extern float P[STATE_SPACE][STATE_SPACE];
extern float y[EXTRA_STATE_SPACE];
extern float u[INPUT_SPACE];
extern long secondsElapsed;
// ------------------- END -----------------------

#endif
