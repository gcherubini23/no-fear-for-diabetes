void arrayScalarMultiply(float* vec, int length, float scalar);
void arrayAddInPlace(float* vec1, const float* vec2, int length, bool make_non_negative);
float arrayMultiply(const float* vec1, const float* vec2, int length);
void SquareMatrixTranspose(float A[STATE_SPACE][STATE_SPACE], float At[STATE_SPACE][STATE_SPACE]);
void SquareMatrixMultiply(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE], float RES[STATE_SPACE][STATE_SPACE]);
void SquareMatrixAdd(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE]);
void SquareMatrixDiff(float A[STATE_SPACE][STATE_SPACE], float B[STATE_SPACE][STATE_SPACE]);
void ArrayMatrixMultiply(float vec[STATE_SPACE], float A[STATE_SPACE][STATE_SPACE], float res[STATE_SPACE], bool left);
void euler_solve(float* x, const float* u, float* y, float* v, float dt, const float* params);
