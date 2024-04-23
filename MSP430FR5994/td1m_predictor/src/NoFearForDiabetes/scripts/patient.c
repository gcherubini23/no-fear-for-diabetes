#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "config.h"
#include "driverlib.h"
#include "memory_util.h"
#include "model.h"
#include "util.h"


//#pragma CODE_SECTION(init_patient, ".ramfunc")
void init_patient() {
    float daily_basal = 18;
    float basal = daily_basal / 1440;
    int p_names[NUM_TUNING_PARAMS] = {BW,EAT_RATE,BASAL,GB,TD,VG,K1,K2,VI,HEB,M1,M5,KMAX,KMIN,KABS,F,BETA,DELTA,KP1,KP2,KP3,KI,FCNS,KM0,VMX,P2U,KD,KA1,KA2,KE1,KE2};
    float p_values[NUM_TUNING_PARAMS] = {86.1,5,basal,101.6,10,1.88,0.065,0.079,0.05,0.6,0.19,0.0304,0.0558,0.008,0.0570,0.9,0.82,0.01,2.7,0.0021,0.009,0.0079,1,225.59,0.047,0.0331,0.0164,0.0018,0.0182,0.0005,339};
    update_patient(p_values, p_names, TRUE);
}

