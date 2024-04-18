#include <msp430.h>
#include <stdio.h>
#include "driverlib.h"

long lr(long *addr)
{
    return __data20_read_long((unsigned long int)addr);
}

void lw(long *addr, long i)
{
    __data20_write_short((unsigned long int)addr, i);
}

int ir(int *addr)
{
    int temp = __data20_read_short((unsigned long int)addr);
    int *i = (int *)&temp;
    return *i;
}

void iw(int *addr, int f)
{
    int *i = (int *)&f;
    __data20_write_short((unsigned long int)addr, *i);
}

float fr(const float *addr)
{
    unsigned long temp = __data20_read_long((unsigned long int)addr);
    float *f = (float *)&temp;
    return *f;
}

void fw(float *addr, float f)
{
    uint32_t *i = (uint32_t *)&f;
    __data20_write_long((unsigned long int)addr, *i);
}

void readFloatArray(const float* source, float* destination, int size) {
    int i;
    for (i = 0; i < size; i++) {
        destination[i] = fr(source + i);
    }
}

void writeFloatArray(float* addr, const float* values, int size) {
    int i;
    for (i = 0; i < size; i++) {
        fw(addr + i, values[i]);
    }
}

void readFloatMatrix(const float** source, float** destination, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            destination[i][j] = fr(source[i] + j);
        }
    }
}

void writeFloatMatrix(float** addr, const float** values, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fw(addr[i] + j, values[i][j]);
        }
    }
}
