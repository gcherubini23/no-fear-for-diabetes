void arrayScalarMultiply(float* vec, int length, float scalar);
void arrayAddInPlace(float* vec1, const float* vec2, int length, bool make_non_negative);
float arrayMultiply(const float* vec1, const float* vec2, int length);
void matrixAddInPlace(float** A, const float** B, int rows, int cols);
void matrixDiffInPlace(float** A, const float** B, int rows, int cols);
int matrixMultiply(const float** A, int rowsA, int colsA, const float** B, int rowsB, int colsB, float** RES);
void transpose(const float** A, float** A_T, int rows, int cols);
void euler_solve(float* x, const float* u, float* y, float* v, float dt, const float* params);
