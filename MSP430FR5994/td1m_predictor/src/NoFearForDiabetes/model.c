#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "driverlib.h"
#include "memory_util.h"
#include "model_idx.h"

// Define what needs to be in FRAM
#pragma PERSISTENT(params)
float params[NUM_PARAMS];

// ---------- Functions -----------
#pragma CODE_SECTION(update_patient, ".ramfunc")
void update_patient(const float p_values[NUM_TUNING_PARAMS], const int p_names[NUM_TUNING_PARAMS], bool initialization){
    int i;
    for(i=0; i<NUM_TUNING_PARAMS; i++){
        params[p_names[i]] = p_values[i];
    }
    if (initialization){
            params[CL] = 0.0242 * params[BW];
        }
    params[M30] = params[M1] * params[HEB] / (1 - params[HEB]);
    params[M2] = 3 / 5 * params[CL] / (params[HEB] * params[VI] * params[BW]);
    params[M4] = 2 / 5 * params[CL] / (params[VI] * params[BW]);
    params[IPB] = params[U2SS] / (params[M2] + params[M4] - params[M1] * params[M2] / (params[M1] + params[M30]));
    params[ILB] = params[M2] / (params[M1] + params[M30]) * params[IPB];
    params[IB] = params[IPB] / params[VI];
    params[KGRI] = params[KMAX];
    params[GPB] = params[GB] * params[VG];
    params[EGPB] = params[KP1] - params[KP2] * params[GPB] - params[KP3] * params[IB];
    params[GTB] = 1 / params[K2] * (params[FCNS] - params[EGPB] + params[K1] * params[GPB]);
    params[VM0] = (params[EGPB] - params[FCNS]) * (params[KM0] + params[GTB]) / params[GTB];
    params[ISC1SS] = params[U2SS] / (params[KD] + params[KA1]);
    params[ISC2SS] = params[ISC1SS] * params[KD] / params[KA2];

}

#pragma CODE_SECTION(step, ".ramfunc")
void step(const float x[STATE_SPACE], const float y[EXTRA_STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE]){
    gastro_intestinal_tract(x, y, v, dx_dt);
    glucose_subystem(x, dx_dt);
    insulin_infusion_subsystem(x, v, dx_dt);
}

#pragma CODE_SECTION(preprocess, ".ramfunc")
// Input: [x_k,y_kminus1,u_k,dt] --> Output: [y_k,v_k]
void preprocess(const float x[STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], const float dt){

    float IIR_dt;
    if (y[INS_TO_INF] <= 0) {
        IIR_dt = u[IIR];
    } else {
        if (u[IIR] >= y[LAST_IIR]) {
            IIR_dt = u[IIR];
        } else {
            IIR_dt = y[LAST_IIR];
        }
    }
    float insulin_to_infuse = y[INS_TO_INF] + u[IIR];
    if (insulin_to_infuse < IIR_dt * dt && dt != 0) {
        IIR_dt = insulin_to_infuse / dt;
        insulin_to_infuse = 0;
    } else {
        insulin_to_infuse = insulin_to_infuse - IIR_dt * dt;
        if (insulin_to_infuse < 0.00001) {
            insulin_to_infuse = 0;
        }
    }

    float CHO_consumed_rate;
    if (dt == 0) {
        CHO_consumed_rate = 0;
    } else if ((y[CHO_TO_EAT] / dt >= params[EAT_RATE]) || (u[CHO] / dt >= params[EAT_RATE] && y[CHO_TO_EAT] == 0)) {
        CHO_consumed_rate = params[EAT_RATE];
    } else if ((u[CHO] > 0) && (u[CHO] / dt < params[EAT_RATE]) && (y[CHO_TO_EAT] == 0)) {
        CHO_consumed_rate = u[CHO] / dt;
    } else {
        CHO_consumed_rate = y[CHO_TO_EAT] / dt;
    }
    float CHO_to_eat = u[CHO] + y[CHO_TO_EAT] - CHO_consumed_rate * dt;

    float new_D;
    if (y[IS_EATING]) {
        new_D = CHO_consumed_rate * dt + y[D];
    } else if (CHO_consumed_rate > 0) {
        new_D = CHO_consumed_rate * dt;
    } else {
        new_D = y[D];
    }

    float lastQsto, is_eating;
    if ((CHO_consumed_rate > 0) && !y[IS_EATING]) {
        is_eating = 1.0;
        lastQsto = x[Q_STO1] + x[Q_STO2];
    } else if ((CHO_consumed_rate == 0) && y[IS_EATING]) {
        is_eating = 0.0;
        lastQsto = y[LAST_Q_STO];
    } else {
        is_eating = y[IS_EATING];
        lastQsto = y[LAST_Q_STO];
    }

    v[CHO_DT] = CHO_consumed_rate;
    v[IIR_DT] = IIR_dt;
    y[INS_TO_INF] = insulin_to_infuse;
    y[LAST_IIR] = IIR_dt;
    y[CHO_TO_EAT] = CHO_to_eat;
    y[D] = new_D;
    y[LAST_Q_STO] = lastQsto;
    y[IS_EATING] = is_eating;
}

#pragma CODE_SECTION(gastro_intestinal_tract, ".ramfunc")
void gastro_intestinal_tract(const float x[STATE_SPACE], const float y[EXTRA_STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE]){
    dx_dt[Q_STO1] = - params[KGRI] * x[Q_STO1] + v[CHO_DT] * 1000;

    float Q_sto = x[Q_STO1] + x[Q_STO2];
    float D_bar = y[LAST_Q_STO] + y[D] * 1000;
    float kgut;
    if (D_bar > 0) {
        float aa = 5 / (2 * (1 - params[BETA]) * D_bar);
        float cc = 5 / (2 * params[DELTA] * D_bar);
        kgut = params[KMIN] + (params[KMAX] - params[KMIN]) / 2 * (tanh(aa * (Q_sto - params[BETA] * D_bar)) - tanh(cc * (Q_sto - params[DELTA] * D_bar)) + 2);
    } else {
        kgut = params[KMAX];
    }

    dx_dt[Q_STO2] = params[KGRI] * x[Q_STO1] - kgut * x[Q_STO2];

    dx_dt[Q_GUT] = kgut * x[Q_STO2] - params[KABS] * x[Q_GUT];
}

#pragma CODE_SECTION(glucose_subystem, ".ramfunc")
void glucose_subystem(const float x[STATE_SPACE], float dx_dt[STATE_SPACE]){



}

#pragma CODE_SECTION(insulin_infusion_subsystem, ".ramfunc")
void insulin_infusion_subsystem(const float x[STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE]) {



}



