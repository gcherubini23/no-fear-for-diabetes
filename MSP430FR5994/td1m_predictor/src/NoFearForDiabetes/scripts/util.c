#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include "driverlib.h"
#include "memory_util.h"
#include "model.h"
#include "config.h"


void arrayScalarMultiply(float* vec, int length, float scalar) {
    int i;
    for (i = 0; i < length; i++) {
        vec[i] *= scalar;
        }
    }

void arrayAddInPlace(float* vec1, const float* vec2, int length, bool make_non_negative) {
    int i;
    for (i = 0; i < length; i++) {
            vec1[i] += vec2[i];
            if ((make_non_negative) && (vec1[i] < 0)) {
                vec1[i] = 0;
                }
            }
}

void matrixAddInPlace(float** A, const float** B, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            A[i][j] += B[i][j];
        }
    }
}

int matrixMultiply(const float** A, int rowsA, int colsA, const float** B, int rowsB, int colsB, float** RES) {
    if (colsA != rowsB) {
        return FALSE;
    }
    int i, j, k;
    for (i = 0; i < rowsA; i++) {
        for (j = 0; j < colsB; j++) {
            RES[i][j] = 0;
            for (k = 0; k < colsA; k++) {
                RES[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return TRUE;
}

void euler_solve(float x[STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], float dt) {
    preprocess(x, u, y, v, dt);
    float dx_dt[STATE_SPACE];
    step(x, y, v, dx_dt);
    arrayScalarMultiply(dx_dt, STATE_SPACE, dt);
    arrayAddInPlace(x, dx_dt, STATE_SPACE, TRUE);
}
