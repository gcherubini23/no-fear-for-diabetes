#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
//#include <time.h>
#include "memory_util.h"
#include "config.h"


bool sample_measurement(float* z) {
    *z = 110;
    long t;
    t = lr(&secondsElapsed);
    if (t >= UPDATE_RATE * SEC_IN_MIN) {
        lw(&secondsElapsed, 0);
        return TRUE;
    } else {
        return FALSE;
    }
}

bool sample_input(float* u) {
    long t;
    u[CHO] = 0;
    u[IIR] = 0;
    t = lr(&secondsElapsed);
    if (t >= UPDATE_RATE * SEC_IN_MIN) {
            return TRUE;
    } else {
        return FALSE;
    }
}
