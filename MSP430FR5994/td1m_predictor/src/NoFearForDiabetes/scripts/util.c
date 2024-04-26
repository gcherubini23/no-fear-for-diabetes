#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "driverlib.h"
#include "memory_util.h"
#include "model.h"
#include "config.h"


//#pragma CODE_SECTION(arrayScalarMultiply, ".ramfunc")
void arrayScalarMultiply(float* vec, int length, float scalar) {
    int i;
    for (i = 0; i < length; i++) {
        vec[i] *= scalar;
        }
    }

//#pragma CODE_SECTION(arrayAddInPlace, ".ramfunc")
void arrayAddInPlace(float* vec1, const float* vec2, int length, bool make_non_negative) {
    int i;
    for (i = 0; i < length; i++) {
            vec1[i] += vec2[i];
            if ((make_non_negative) && (vec1[i] < 0)) {
                vec1[i] = 0;
                }
            }
}

//#pragma CODE_SECTION(arrayMultiply, ".ramfunc")
float arrayMultiply(const float* vec1, const float* vec2, int length) {
    int i;
    float sum = 0;
    for (i = 0; i < length; i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}

//#pragma CODE_SECTION(SquareMatrixTranspose, ".ramfunc")
void SquareMatrixTranspose(float A[STATE_SPACE][STATE_SPACE], float At[STATE_SPACE][STATE_SPACE]) {
    int i,j;
    for(i = 0; i < STATE_SPACE; i++) {
        for(j = 0; j < STATE_SPACE; j++) {
            At[i][j] = A[j][i];
        }
    }
}

//#pragma CODE_SECTION(SquareMatrixMultiply, ".ramfunc")
void SquareMatrixMultiply(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE], float RES[STATE_SPACE][STATE_SPACE]) {
    int i, j, k;
    for (i = 0; i < STATE_SPACE; i++) {
        for (j = 0; j < STATE_SPACE; j++) {
            RES[i][j] = 0;
            for (k = 0; k < STATE_SPACE; k++) {
                RES[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

//#pragma CODE_SECTION(SquareMatrixAdd, ".ramfunc")
void SquareMatrixAdd(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE]) {
    int i,j;
    for(i = 0; i < STATE_SPACE; i++) {
        for(j = 0; j < STATE_SPACE; j++) {
            A[i][j] += B[j][i];
        }
    }
}

//#pragma CODE_SECTION(SquareMatrixDiff, ".ramfunc")
void SquareMatrixDiff(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE]) {
    int i,j;
    for(i = 0; i < STATE_SPACE; i++) {
        for(j = 0; j < STATE_SPACE; j++) {
            A[i][j] -= B[j][i];
        }
    }
}

//#pragma CODE_SECTION(ArrayMatrixMultiply, ".ramfunc")
void ArrayMatrixMultiply(float vec[STATE_SPACE], float A[STATE_SPACE][STATE_SPACE], float res[STATE_SPACE], bool left) {
    float val;
    int i, j;
    for (i = 0; i < STATE_SPACE; i++) {
        res[i] = 0;
        for (j = 0; j < STATE_SPACE; j++) {
            if (left) {
                val = A[j][i];
            } else {
                val = A[i][j];
            }
            res[i] += vec[j] * val;
        }
    }
}

//#pragma CODE_SECTION(euler_solve, ".ramfunc")
void euler_solve(float x[STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], float dt, const float params[NUM_PARAMS]) {
    preprocess(x, u, y, v, dt);
    float dx_dt[STATE_SPACE];
    step(x, y, v, dx_dt, params);
    arrayScalarMultiply(dx_dt, STATE_SPACE, dt);
    arrayAddInPlace(x, dx_dt, STATE_SPACE, TRUE);
}

