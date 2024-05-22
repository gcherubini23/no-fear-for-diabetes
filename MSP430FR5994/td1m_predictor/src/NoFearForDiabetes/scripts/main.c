#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "driverlib.h"
#include "memory_util.h"
#include "model.h"
#include "util.h"
#include "ekf.h"
#include "patient.h"
#include "auxiliary_functions.h"

#pragma PERSISTENT(x)
float x[STATE_SPACE] = { 0 };

#pragma PERSISTENT(P)
float P[STATE_SPACE][STATE_SPACE] = { 0 };

#pragma PERSISTENT(y)
float y[EXTRA_STATE_SPACE] = { 0 };

#pragma PERSISTENT(u)
float u[INPUT_SPACE] = { 0 };

#pragma PERSISTENT(secondsElapsed)
long secondsElapsed = 0;


void init_MCU()
{
    // Stop watchdog timer
    WDTCTL = WDTPW | WDTHOLD;
}

void init_timer() {
    TA0CCTL0 = CCIE;
    TA0CTL = TASSEL__ACLK | MC__UP;
    TA0CCR0 = ACLK_HZ - 1;
}

void init_states() {
    fw(&x[Q_STO1], 0);
    fw(&x[Q_STO2], 0);
    fw(&x[Q_GUT], 0);
    fw(&x[G_P], fr(&params[GPB]));
    fw(&x[G_T], fr(&params[GTB]));
    fw(&x[G_SC], fr(&params[GPB]));
    fw(&x[I_L], fr(&params[ILB]));
    fw(&x[I_P], fr(&params[IPB]));
    fw(&x[I_1], fr(&params[IB]));
    fw(&x[I_D], fr(&params[IB]));
    fw(&x[X], 0);
    fw(&x[I_SC1], fr(&params[ISC1SS]));
    fw(&x[I_SC2], fr(&params[ISC2SS]));

    fw(&P[Q_STO1][Q_STO1], MODEL_COV * 0.01 * 3);
    fw(&P[Q_STO2][Q_STO2], MODEL_COV * 0.01 * 3);
    fw(&P[Q_GUT][Q_GUT], MODEL_COV * 0.01 * 3);
    fw(&P[G_P][G_P], MODEL_COV * 3);
    fw(&P[G_T][G_T], MODEL_COV * 3);
    fw(&P[G_SC][G_SC], MODEL_COV * 3);

    int j;
    for(j = 0; j < EXTRA_STATE_SPACE; j++) {
        fw(&y[j], 0);
    }
}

void ekf_on() {
    float z_k, u_k[INPUT_SPACE];
    bool new_measurement, new_input;

    new_measurement = sample_measurement(&z_k);
    new_input = sample_input(u_k);
    if (new_measurement || new_input) {
        float x_[STATE_SPACE], P_[STATE_SPACE][STATE_SPACE], y_[STATE_SPACE], u_[INPUT_SPACE], horizon;

        readFloatArray(x, x_, STATE_SPACE);
        readFloatArray(y, y_, EXTRA_STATE_SPACE);
        memcpy(P_, P, STATE_SPACE * STATE_SPACE * sizeof(float));
        readFloatArray(u, u_, INPUT_SPACE);
        horizon = lr(&secondsElapsed) / SEC_IN_MIN;
        predict(x_, P_, u_, y_, horizon);

        if (new_measurement) {
            float residual, innovation_cov;
            update_measurement_cov(z_k);
            measurement_update(x_, P_, z_k, &residual, &innovation_cov);
            /* To implement anomaly detector */
        }

        writeFloatArray(x, x_, STATE_SPACE);
        writeFloatArray(y, y_, EXTRA_STATE_SPACE);
        memcpy(P, P_, STATE_SPACE * STATE_SPACE * sizeof(float));
        writeFloatArray(u, u_k, INPUT_SPACE);
    }
}

int main(void) {
    init_MCU();
    init_patient();
    init_ekf();
    init_states();
    init_model();
    init_timer();

    _BIS_SR(LPM3_bits + GIE);

//    while(1) {
//        ekf_on();
//        lw(&secondsElapsed, 60);
//    }
}

#pragma vector=TIMER0_A0_VECTOR
__interrupt void Timer_A0_ISR(void) {
    long time, startTick, endTick, lastExecutionTime;
    time = lr(&secondsElapsed);
    lw(&secondsElapsed, time+1);
    if(time >= UPDATE_RATE) {
//        startTick = TA0R;
        ekf_on();
//        endTick = TA0R;
//        lastExecutionTime = (endTick - startTick) * 10 / ACLK_HZ;
////        lastExecutionTime = (endTick >= startTick) ? (endTick - startTick) : (ACLK_HZ - startTick + endTick);
//        lw(&secondsElapsed, lastExecutionTime);
        lw(&secondsElapsed, 0);
    }
}
