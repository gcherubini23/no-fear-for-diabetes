#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "memory_util.h"
#include "config.h"


bool sample_measurement(float* z) {
    *z = 110;
    return TRUE;
//    return FALSE;
}

bool sample_input(float* u) {
    u[CHO] = 0;
    u[IIR] = 0;
//    return TRUE;
    return TRUE;
}
