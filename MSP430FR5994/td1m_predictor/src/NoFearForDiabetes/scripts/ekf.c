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
    float v_g = fr(&params[VG]);
    fw(&mQ[Q_STO1][Q_STO1], MODEL_COV * 0.01);
    fw(&mQ[Q_STO2][Q_STO2], MODEL_COV * 0.01);
    fw(&mQ[Q_GUT][Q_GUT], MODEL_COV * 0.01);
    fw(&mQ[G_P][G_P], MODEL_COV);
    fw(&mQ[G_T][G_T], MODEL_COV);
    fw(&mQ[G_SC][G_SC], MODEL_COV);

    fw(&mR, SENSOR_COV);

    fw(&mH[G_SC], 1 / v_g);
}

//#pragma CODE_SECTION(update_measurement_cov, ".ramfunc")
void update_measurement_cov(float z) {
    float cov = (CGM_MARD / 100 * z) * (CGM_MARD / 100 * z);
    fw(&mR, cov);
}

//#pragma CODE_SECTION(process_update, ".ramfunc")
void process_update(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], float dt, const float params[NUM_PARAMS]) {
    float xkmin1[STATE_SPACE], A_T[STATE_SPACE][STATE_SPACE], A[STATE_SPACE][STATE_SPACE], temp1[STATE_SPACE][STATE_SPACE], Q[STATE_SPACE][STATE_SPACE];
    int i, j;
    memcpy(xkmin1, x, STATE_SPACE * sizeof(float));
    euler_solve(x, u, y, v, dt, params);
    linearize(xkmin1, y, FALSE, params);
    memcpy(A, mA, STATE_SPACE * STATE_SPACE * sizeof(float));
    memcpy(Q, mQ, STATE_SPACE * STATE_SPACE * sizeof(float));
    SquareMatrixScalarMultiply(A, dt);
    SquareMatrixScalarMultiply(Q, dt);
    for (i = 0; i < STATE_SPACE; i++) {
            A[i][i] += 1;
    }
    SquareMatrixTranspose(A, A_T);
    SquareMatrixMultiply(A, P, temp1);
    SquareMatrixMultiply(temp1, A_T, P);
    SquareMatrixAdd(P, Q);;
}

//#pragma CODE_SECTION(measurement_update, ".ramfunc")
void measurement_update(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], float z, float* residual, float* innovation_cov) {
    float H[STATE_SPACE], P_Ht[STATE_SPACE], KRK[STATE_SPACE][STATE_SPACE], temp[STATE_SPACE][STATE_SPACE], eyet[STATE_SPACE][STATE_SPACE];
    float eye[STATE_SPACE][STATE_SPACE] = {0};
    float* K;
    float cov, weight, observed_state;
    float r = fr(&mR);
    int i, j;
    readFloatArray(mH, H, STATE_SPACE);

    ArrayMatrixMultiply(H, P, P_Ht, FALSE);
    cov = arrayMultiply(P_Ht, H, STATE_SPACE);
    *innovation_cov = cov + r;
    weight = 1 / *innovation_cov;
    arrayScalarMultiply(P_Ht, STATE_SPACE, weight);
    K = P_Ht;
    observed_state = arrayMultiply(H, x, STATE_SPACE);
    *residual = z - observed_state;
    for(i = 0; i < STATE_SPACE; i++) {
        x[i] += K[i] * *residual;
        eye[i][i] = 1;
    }
    for(i = 0; i < STATE_SPACE; i++) {
        for(j = 0; j < STATE_SPACE; j++) {
            eye[i][j] -= K[i] * H[j];
            KRK[i][j] = K[i] * r * K[j];
        }
    }
    SquareMatrixMultiply(eye, P, temp);
    SquareMatrixTranspose(eye, eyet);
    SquareMatrixMultiply(temp, eyet, P);
    SquareMatrixAdd(P, KRK);
}

//#pragma CODE_SECTION(predict, ".ramfunc")
void predict(float x[STATE_SPACE], float P[STATE_SPACE][STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float horizon) {
    float patient[NUM_PARAMS];
    float t = 0;
    float v[MODEL_INPUT_SPACE] = { 0 };
    float input_u[INPUT_SPACE] = { 0 };
    readFloatArray(params, patient, NUM_PARAMS);
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
