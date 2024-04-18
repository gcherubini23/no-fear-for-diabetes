#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "config.h"
#include "driverlib.h"
#include "memory_util.h"
#include "model.h"
#include "util.h"

// Define what needs to be in FRAM
#pragma PERSISTENT(Q)
float Q[STATE_SPACE][STATE_SPACE] = { 0 };

#pragma PERSISTENT(R)
float R = 0;

#pragma PERSISTENT(K)
float K[STATE_SPACE] = { 0 };

#pragma PERSISTENT(H)
float H[STATE_SPACE] = { 0 };

// ---------- Functions -----------
//#pragma CODE_SECTION(init_ekf, ".ramfunc")
void init_ekf() {

}

//#pragma CODE_SECTION(process_update, ".ramfunc")
void process_update() {

}

//#pragma CODE_SECTION(measurement_update, ".ramfunc")
void measurement_update() {

}

//#pragma CODE_SECTION(predict, ".ramfunc")
void predict() {

}
