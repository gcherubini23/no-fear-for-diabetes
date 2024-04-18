void arrayScalarMultiply(float* vec, int length, float scalar);
void arrayAddInPlace(float* vec1, const float* vec2, int length, bool make_non_negative);
void matrixAddInPlace(float** A, const float** B, int rows, int cols);
int matrixMultiply(const float** A, int rowsA, int colsA, const float** B, int rowsB, int colsB, float** RES);
void euler_solve(float* x, const float* u, float* y, float* v, float dt);
