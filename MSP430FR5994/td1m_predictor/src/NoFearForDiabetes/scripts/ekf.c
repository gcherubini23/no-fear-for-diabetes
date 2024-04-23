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

// Define what needs to be in FRAM
#pragma PERSISTENT(mQ)
float mQ[STATE_SPACE][STATE_SPACE] = { 0 };

#pragma PERSISTENT(mR)
float mR = 0;

#pragma PERSISTENT(mH)
float mH[STATE_SPACE] = { 0 };

// ---------- Functions -----------
//#pragma CODE_SECTION(init_ekf, ".ramfunc")
void init_ekf() {
    fw(&mQ[G_P][G_P], MODEL_COV);
    fw(&mQ[G_T][G_T], MODEL_COV);
    fw(&mQ[G_SC][G_SC], MODEL_COV);

    fw(&mR, SENSOR_COV);

    float v_g = fr(&params[VG]);
    fw(&mH[G_SC], 1 / v_g);
}

//#pragma CODE_SECTION(update_measurement_cov, ".ramfunc")
void update_measurement_cov(float z) {
    float cov = (CGM_MARD / 100 * z) * (CGM_MARD / 100 * z);
    fw(&mR, cov);
}

//#pragma CODE_SECTION(process_update, ".ramfunc")
void process_update(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], float dt, const float params[NUM_PARAMS]) {
    float xkmin1[STATE_SPACE];
    int i;
    for (i = 0; i < STATE_SPACE; i++) {
        xkmin1[i] = x[i];
    }
    euler_solve(x, u, y, v, dt, params);
    linearize(xkmin1, y, FALSE, params);
    float A_T[STATE_SPACE][STATE_SPACE], A[STATE_SPACE][STATE_SPACE];
    readFloatMatrix(mA, A, STATE_SPACE, STATE_SPACE);
    transpose(A, A_T, STATE_SPACE, STATE_SPACE);
    float temp1[STATE_SPACE][STATE_SPACE];
    matrixMultiply(A, STATE_SPACE, STATE_SPACE, P, STATE_SPACE, STATE_SPACE, temp1);
    matrixMultiply(temp1, STATE_SPACE, STATE_SPACE, A_T, STATE_SPACE, STATE_SPACE, P);
    float Q[STATE_SPACE][STATE_SPACE];
    readFloatMatrix(mQ, Q, STATE_SPACE, STATE_SPACE);
    matrixAddInPlace(P, Q, STATE_SPACE, STATE_SPACE);
}

//#pragma CODE_SECTION(measurement_update, ".ramfunc")
void measurement_update(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], float z, float* residual, float* innovation_cov) {
    float H[1][STATE_SPACE];
    readFloatArray(mH, H[0], STATE_SPACE);
    float R[1][1] = {fr(&mR)};

    float H_T[STATE_SPACE][1];
    transpose(H, H_T, 1, STATE_SPACE);
    float P_H_T[STATE_SPACE][1];
    matrixMultiply(P, STATE_SPACE, STATE_SPACE, H_T, STATE_SPACE, 1, P_H_T);
    float temp1[1][1];
    matrixMultiply(H, 1, STATE_SPACE, P_H_T, STATE_SPACE, 1, temp1);
    matrixAddInPlace(temp1, R, 1, 1);
    *innovation_cov = temp1[0][0];

    float temp2[1][1] = {1 / *innovation_cov};
    float K[STATE_SPACE][1];
    matrixMultiply(P_H_T, STATE_SPACE, 1, temp2, 1, 1, K);
    *residual = z - arrayMultiply(H[0], x, STATE_SPACE);

    float eye[STATE_SPACE][STATE_SPACE] = { 0 };
    int i;
    for (i = 0; i < STATE_SPACE; i++) {
        x[i] += K[i][0] * *residual;
        eye[i][i] = 1;
    }

    float K_H[STATE_SPACE][STATE_SPACE];
    matrixMultiply(K, STATE_SPACE, 1, H, 1, STATE_SPACE, K_H);
    matrixDiffInPlace(eye, K_H, STATE_SPACE, STATE_SPACE);
    float eye_T[STATE_SPACE][STATE_SPACE];
    transpose(eye, eye_T, STATE_SPACE, STATE_SPACE);
    float temp3[STATE_SPACE][STATE_SPACE];
    matrixMultiply(eye, STATE_SPACE, STATE_SPACE, P, STATE_SPACE, STATE_SPACE, temp3);
    matrixMultiply(temp3, STATE_SPACE, STATE_SPACE, eye_T, STATE_SPACE, STATE_SPACE, P);

    float temp4[STATE_SPACE][1];
    matrixMultiply(K, STATE_SPACE, 1, R, 1, 1, temp4);
    float K_T[1][STATE_SPACE], temp5[STATE_SPACE][STATE_SPACE];
    transpose(K, K_T, STATE_SPACE, 1);
    matrixMultiply(temp4, STATE_SPACE, 1, K_T, 1, STATE_SPACE, temp5);

    matrixDiffInPlace(P, temp5, STATE_SPACE, STATE_SPACE);
}

//#pragma CODE_SECTION(predict, ".ramfunc")
void predict(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float horizon) {
    float patient[NUM_PARAMS];
    readFloatArray(params, patient, NUM_PARAMS);
    float t = 0;
    float v[MODEL_INPUT_SPACE] = { 0 };
    float input_u[INPUT_SPACE] = { 0 };
    input_u[CHO] = u[CHO];
    input_u[IIR] = u[IIR];
    while (t < horizon) {
        float step_dt = DT;
        float time_left = horizon - t;
        if (time_left < DT) {
            step_dt = time_left;
        }
        process_update(x, P, input_u, y, v, step_dt, patient);
        input_u[CHO] = 0;
        input_u[IIR] = 0;
        t += step_dt;
    }
}
