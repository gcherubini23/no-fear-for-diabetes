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

#pragma PERSISTENT(t)
float t = 0;

#pragma PERSISTENT(last_update_t)
float last_update_t = 0;


void init_MCU()
{
    // Stop watchdog timer
    WDTCTL = WDTPW | WDTHOLD;

    P1DIR = 0xFF;
    P1OUT = 0x00;
    P2DIR = 0xFF;
    P2OUT = 0x00;
    P3DIR = 0xFF;
    P3OUT = 0x00;
    P4DIR = 0xFF;
    P4OUT = 0x00;
    P5DIR = 0xFF;
    P5OUT = 0x00;
    P6DIR = 0xFF;
    P6OUT = 0x00;
    P7DIR = 0xFF;
    P7OUT = 0x00;
    P8DIR = 0xFF;
    P8OUT = 0x00;
    PADIR = 0xFF;
    PAOUT = 0x00;
    PBDIR = 0xFF;
    PBOUT = 0x00;
    PCDIR = 0xFF;
    PCOUT = 0x00;
    PDDIR = 0xFF;
    PDOUT = 0x00;

#if 1
    // Configure one FRAM waitstate as required by the device datasheet for MCLK
    // operation beyond 8MHz _before_ configuring the clock system.
    FRCTL0 = FRCTLPW | NWAITS_1;

    // Clock System Setup
    CSCTL0_H = CSKEY_H;                     // Unlock CS registers
    CSCTL1 = DCOFSEL_0;                     // Set DCO to 1MHz
    // Set SMCLK = MCLK = DCO, ACLK = VLOCLK
    CSCTL2 = SELA__VLOCLK | SELS__DCOCLK | SELM__DCOCLK;
    // Per Device Errata set divider to 4 before changing frequency to
    // prevent out of spec operation from overshoot transient
    CSCTL3 = DIVA__4 | DIVS__4 | DIVM__4;   // Set all corresponding clk sources to divide by 4 for errata
    CSCTL1 = DCOFSEL_4 | DCORSEL;           // Set DCO to 16MHz
    // Delay by ~10us to let DCO settle. 60 cycles = 20 cycles buffer + (10us / (1/4MHz))
    __delay_cycles(60);
    CSCTL3 = DIVA__1 | DIVS__1 | DIVM__1;   // Set all dividers to 1 for 16MHz operation
    CSCTL0_H = 0;
#endif

    PMM_unlockLPM5();
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

    fw(&P[Q_STO1][Q_STO1], MODEL_COV * 0.01);
    fw(&P[Q_STO2][Q_STO2], MODEL_COV * 0.01);
    fw(&P[Q_GUT][Q_GUT], MODEL_COV * 0.01);
    fw(&P[G_P][G_P], MODEL_COV);
    fw(&P[G_T][G_T], MODEL_COV);
    fw(&P[G_SC][G_SC], MODEL_COV);

    int j;
    for(j = 0; j < EXTRA_STATE_SPACE; j++) {
        fw(&y[j], 0);
    }
}

void ekf_on() {
    float z_k;
    float u_k[INPUT_SPACE];
    float time = current_time();
    fw(&t, time);
    bool new_measurement = sample_measurement(&z_k);
    bool new_input = sample_input(u_k);
    if (new_measurement || new_input) {
        float x_[STATE_SPACE];
        float P_[STATE_SPACE][STATE_SPACE];
        float y_[STATE_SPACE];
        float u_[INPUT_SPACE];
        readFloatArray(x, x_, STATE_SPACE);
        readFloatArray(y, y_, EXTRA_STATE_SPACE);
//        memcpy(P_, P, sizeof(P));
        readFloatMatrix(P, P_, STATE_SPACE, STATE_SPACE);   // not copying
        readFloatArray(u, u_, INPUT_SPACE);

        float horizon = fr(&t) - fr(&last_update_t);
        predict(x_, P_, u_, y_, horizon);
        if (new_measurement) {
            float residual, innovation_cov;
            measurement_update(x_, P_, z_k, &residual, &innovation_cov);
            /* To implement anomaly detector */
        }
        writeFloatArray(x, x_, STATE_SPACE);
        writeFloatArray(y, y_, EXTRA_STATE_SPACE);
        writeFloatMatrix(P, P_, STATE_SPACE, STATE_SPACE);
        writeFloatArray(u, u_k, INPUT_SPACE);
        fw(&last_update_t, fr(&t));
    }
}

int main(void) {
    init_MCU();
    init_patient();
    init_ekf();
    init_states();

    /* LED for test */
    GPIO_setAsOutputPin(GPIO_PORT_P1, GPIO_PIN0);
    GPIO_setOutputLowOnPin(GPIO_PORT_P1, GPIO_PIN0);
    GPIO_setAsOutputPin(GPIO_PORT_P1, GPIO_PIN1);
    GPIO_setOutputHighOnPin(GPIO_PORT_P1, GPIO_PIN1);

    /* setup interrupt */
    GPIO_setAsInputPinWithPullDownResistor(GPIO_PORT_P1, GPIO_PIN3);
    GPIO_selectInterruptEdge(GPIO_PORT_P1, GPIO_PIN3, GPIO_LOW_TO_HIGH_TRANSITION);
    GPIO_clearInterrupt(GPIO_PORT_P1, GPIO_PIN3);
    GPIO_enableInterrupt(GPIO_PORT_P1, GPIO_PIN3);

    while (true) {
        ekf_on();
    }


    _BIS_SR(LPM4_bits + GIE);
}

#pragma vector=PORT1_VECTOR
__interrupt void Port_1(void)
{
    /* LED for test */
    GPIO_setOutputHighOnPin(GPIO_PORT_P1, GPIO_PIN0);
    GPIO_setOutputLowOnPin(GPIO_PORT_P1, GPIO_PIN1);

    while (true) {
        uint8_t status = GPIO_getInputPinValue(GPIO_PORT_P1, GPIO_PIN3);
        if (status == 0)
            break;
        ekf_on();
    }

    /* LED for test */
    GPIO_setOutputLowOnPin(GPIO_PORT_P1, GPIO_PIN0);
    GPIO_setOutputHighOnPin(GPIO_PORT_P1, GPIO_PIN1);

    GPIO_clearInterrupt(GPIO_PORT_P1, GPIO_PIN3);
}
