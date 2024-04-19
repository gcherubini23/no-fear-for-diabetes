bool update_patient(const float* p_values, const int* p_names, bool initialization);
void preprocess(const float* x, const float* u, float* y, float* v, const float dt);
void gastro_intestinal_tract(const float* x, const float* y, const float* v, float* dx_dt, const float* params);
void glucose_subystem(const float* x, float* dx_dt, const float* params);
void insulin_infusion_subsystem(const float* x, const float* v, float* dx_dt, const float* params);
void step(const float* x, const float* y, const float* v, float* dx_dt, const float* params);
void init_condition(float* x, float* y);
void linearize(const float* x, const float* y, bool patient_updated, const float* params);
