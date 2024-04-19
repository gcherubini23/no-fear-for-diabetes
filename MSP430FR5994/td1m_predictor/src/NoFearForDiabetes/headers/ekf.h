void init_ekf();
void update_measurement_cov(float z);
void process_update(float* x, float** P, const float* u, float* y, float* v, float dt, const float* params);
void measurement_update(float* x, float** P, float z, float* residual, float* innovation_cov);
void predict(float* x, float** P, const float* u, float* y, float horizon);
