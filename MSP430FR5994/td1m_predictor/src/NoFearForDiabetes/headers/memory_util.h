long lr(long *addr);
void lw(long *addr, long i);
int ir(int *addr);
void iw(int *addr, int f);
float fr(const float *addr);
void fw(float *addr, float f);
void readFloatArray(const float* source, float* destination, int size);
void writeFloatArray(float* addr, const float* values, int size);
void readFloatMatrix(const float** source, float** destination, int rows, int cols);
void writeFloatMatrix(float** addr, const float** values, int rows, int cols);

