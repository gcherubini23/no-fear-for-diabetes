#include <msp430.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "driverlib.h"
#include "memory_util.h"
#include "config.h"

// Define what needs to be in FRAM
#pragma PERSISTENT(params)
float params[NUM_PARAMS] = { 0 };

#pragma PERSISTENT(mA)
float mA[STATE_SPACE][STATE_SPACE] = { 0 };

#pragma PERSISTENT(mB)
float mB[STATE_SPACE][MODEL_INPUT_SPACE] = { 0 };

#pragma PERSISTENT(mD)
float mD[STATE_SPACE] = { 0 };


// ---------- Functions -----------
//#pragma CODE_SECTION(update_patient, ".ramfunc")
bool update_patient(const float p_values[NUM_TUNING_PARAMS], const int p_names[NUM_TUNING_PARAMS], bool initialization) {

    float patient[NUM_PARAMS];
    readFloatArray(params, patient, NUM_PARAMS);

    int i;
    for(i=0; i<NUM_TUNING_PARAMS; i++){
        patient[p_names[i]] = p_values[i];
    }
    if (initialization){
        patient[CL] = 0.0242 * patient[BW];
        }
    patient[M30] = patient[M1] * patient[HEB] / (1 - patient[HEB]);
    patient[M2] = 3 / 5 * patient[CL] / (patient[HEB] * patient[VI] * patient[BW]);
    patient[M4] = 2 / 5 * patient[CL] / (patient[VI] * patient[BW]);
    patient[IPB] = patient[U2SS] / (patient[M2] + patient[M4] - patient[M1] * patient[M2] / (patient[M1] + patient[M30]));
    patient[ILB] = patient[M2] / (patient[M1] + patient[M30]) * patient[IPB];
    patient[IB] = patient[IPB] / patient[VI];
    patient[KGRI] = patient[KMAX];
    patient[GPB] = patient[GB] * patient[VG];
    patient[EGPB] = patient[KP1] - patient[KP2] * patient[GPB] - patient[KP3] * patient[IB];
    patient[GTB] = 1 / patient[K2] * (patient[FCNS] - patient[EGPB] + patient[K1] * patient[GPB]);
    patient[VM0] = (patient[EGPB] - patient[FCNS]) * (patient[KM0] + patient[GTB]) / patient[GTB];
    patient[ISC1SS] = patient[U2SS] / (patient[KD] + patient[KA1]);
    patient[ISC2SS] = patient[ISC1SS] * patient[KD] / patient[KA2];

    writeFloatArray(params, patient, NUM_PARAMS);

    return TRUE;
}

//#pragma CODE_SECTION(preprocess, ".ramfunc")
// Input: [x_k,y_kminus1,u_k,dt] --> Output: [y_k,v_k]
void preprocess(const float x[STATE_SPACE], const float u[INPUT_SPACE], float y[EXTRA_STATE_SPACE], float v[MODEL_INPUT_SPACE], const float dt) {

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
        if (insulin_to_infuse < EPS) {
            insulin_to_infuse = 0;
        }
    }

    float CHO_consumed_rate;
    float eat_rate = fr(&params[EAT_RATE]);
    if (dt == 0) {
        CHO_consumed_rate = 0;
    } else if ((y[CHO_TO_EAT] / dt >= eat_rate) || (u[CHO] / dt >= eat_rate && y[CHO_TO_EAT] == 0)) {
        CHO_consumed_rate = eat_rate;
    } else if ((u[CHO] > 0) && (u[CHO] / dt < eat_rate) && (y[CHO_TO_EAT] == 0)) {
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
        is_eating = TRUE;
        lastQsto = x[Q_STO1] + x[Q_STO2];
    } else if ((CHO_consumed_rate == 0) && y[IS_EATING]) {
        is_eating = FALSE;
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

//#pragma CODE_SECTION(gastro_intestinal_tract, ".ramfunc")
void gastro_intestinal_tract(const float x[STATE_SPACE], const float y[EXTRA_STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE], const float params[NUM_PARAMS]) {
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

//#pragma CODE_SECTION(glucose_subystem, ".ramfunc")
void glucose_subystem(const float x[STATE_SPACE], float dx_dt[STATE_SPACE], const float params[NUM_PARAMS]) {
    float Rat = params[F] * params[KABS] * x[Q_GUT] / params[BW];
    float EGPt = params[KP1] - params[KP2] * x[G_P] - params[KP3] * x[I_D];
    if (EGPt < 0) {
        EGPt = 0;
    }
    float Uiit = params[FCNS];
    float Vmt = params[VM0] + params[VMX] * x[X];
    float Kmt = params[KM0];
    float Uidt = Vmt * x[G_T] / (Kmt + x[G_T]);
    float Et = params[KE1] * (x[G_P] - params[KE2]);
    if (Et < 0) {
        Et = 0;
    }

    dx_dt[G_P] = EGPt + Rat - Uiit - Et - params[K1] * x[G_P] + params[K2] * x[G_T];

    dx_dt[G_T] = -Uidt + params[K1] * x[G_P] - params[K2] * x[G_T];

    dx_dt[G_SC] = -1 / params[TD] * x[G_SC] + 1 / params[TD] * x[G_P];
}

//#pragma CODE_SECTION(insulin_infusion_subsystem, ".ramfunc")
void insulin_infusion_subsystem(const float x[STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE], const float params[NUM_PARAMS]) {
    float insulin = v[IIR_DT] * 6000 / params[BW];

    dx_dt[I_L] = -(params[M1] + params[M30]) * x[I_L] + params[M2] * x[I_P];

    dx_dt[I_SC1] = insulin - (params[KD] + params[KA1]) * x[I_SC1];

    dx_dt[I_SC2] = params[KD] * x[I_SC1] - params[KA2] * x[I_SC2];

    float Rit =  params[KA1] * x[I_SC1] + params[KA2] * x[I_SC2];

    dx_dt[I_P] = -(params[M2] + params[M4]) * x[I_P] + params[M1] * x[I_L] + Rit;
    if (abs(dx_dt[I_P]) <= EPS) {
        dx_dt[I_P] = 0;
    }

    float It = x[I_P] / params[VI];

    dx_dt[X] = params[P2U] * (-x[X] + It - params[IB]);

    dx_dt[I_1] = params[KI] * (It - x[I_1]);

    dx_dt[I_D] = params[KI] * (x[I_1] - x[I_D]);
}

//#pragma CODE_SECTION(step, ".ramfunc")
void step(const float x[STATE_SPACE], const float y[EXTRA_STATE_SPACE], const float v[MODEL_INPUT_SPACE], float dx_dt[STATE_SPACE], const float params[NUM_PARAMS]) {
    gastro_intestinal_tract(x, y, v, dx_dt, params);
    glucose_subystem(x, dx_dt, params);
    insulin_infusion_subsystem(x, v, dx_dt, params);
}

//#pragma CODE_SECTION(init_condition, ".ramfunc")
void init_condition(float x[STATE_SPACE], float y[EXTRA_STATE_SPACE]) {
    x[Q_STO1] = 0;
    x[Q_STO2] = 0;
    x[Q_GUT] = 0;
    x[G_P] = fr(&params[GPB]);
    x[G_T] = fr(&params[GTB]);
    x[G_SC] = fr(&params[GPB]);
    x[I_L] = fr(&params[ILB]);
    x[I_P] = fr(&params[IPB]);
    x[I_1] = fr(&params[IB]);
    x[I_D] = fr(&params[IB]);
    x[X] = 0;
    x[I_SC1] = fr(&params[ISC1SS]);
    x[I_SC2] = fr(&params[ISC2SS]);

    y[INS_TO_INF] = 0;
    y[LAST_IIR] = 0;
    y[CHO_TO_EAT] = 0;
    y[D] = 0;
    y[LAST_Q_STO] = 0;
    y[IS_EATING] = 0;
}

//#pragma CODE_SECTION(linearize, ".ramfunc")
void linearize(const float x[STATE_SPACE], const float y[EXTRA_STATE_SPACE], bool patient_updated, const float patient[NUM_PARAMS]) {
    if (patient_updated) {
        fw(&mA[Q_STO1][Q_STO1], -patient[KGRI]);
        fw(&mA[Q_GUT][Q_GUT], -patient[KABS]);
        fw(&mA[G_P][Q_GUT], patient[F] * patient[KABS] / patient[BW]);
        fw(&mA[G_P][G_T], patient[K2]);
        fw(&mA[G_T][G_P], patient[K1]);
        fw(&mA[G_SC][G_P], 1 / patient[TD]);
        fw(&mA[G_SC][G_SC], -1 / patient[TD]);
        fw(&mA[I_L][I_L], -patient[M1] - patient[M30]);
        fw(&mA[I_L][I_P], patient[M2]);
        fw(&mA[I_P][I_L], patient[M1]);
        fw(&mA[I_P][I_P], -patient[M2] - patient[M4]);
        fw(&mA[I_P][I_SC1], patient[KA1]);
        fw(&mA[I_P][I_SC2], patient[KA2]);
        fw(&mA[I_1][I_P], patient[KI] / patient[VI]);
        fw(&mA[I_1][I_1], -patient[KI]);
        fw(&mA[I_D][I_1], patient[KI]);
        fw(&mA[I_D][I_D], -patient[KI]);
        fw(&mA[X][I_P], patient[P2U] / patient[VI]);
        fw(&mA[X][X], -patient[P2U]);
        fw(&mA[I_SC1][I_SC1], -patient[KD] - patient[KA1]);
        fw(&mA[I_SC2][I_SC1], patient[KD]);
        fw(&mA[I_SC2][I_SC2], -patient[KA2]);

        fw(&mD[X], -patient[IB] * patient[P2U]);

        fw(&mB[Q_STO1][CHO_DT], 1000);
        fw(&mB[I_SC1][IIR_DT], 6000 / patient[BW]);
    }

    float dEGPt_dGp, dEGPt_dId, bias_Gp_EGPt;
    if (patient[KP1] - patient[KP2] * x[G_P] - patient[KP3] * x[I_D] > 0) {
        dEGPt_dGp = -patient[KP2];
        dEGPt_dId = -patient[KP3];
        bias_Gp_EGPt = patient[KP1];
    } else {
        dEGPt_dGp = 0;
        dEGPt_dId = 0;
        bias_Gp_EGPt = 0;
    }

    float dEt_dGp, bias_Gp_Et;
    if (x[G_P] - patient[KE2] > 0) {
        dEt_dGp = patient[KE1];
        bias_Gp_Et = -patient[KE1] * patient[KE2];
    } else {
        dEt_dGp = 0;
        bias_Gp_Et = 0;
    }

    float dUidt_dGt = (patient[VM0] + patient[VMX] * x[X]) * patient[KM0] / ((patient[KM0] + x[G_T]) * (patient[KM0] + x[G_T]));
    float dUidt_dX = patient[VMX] * x[G_T] / (patient[KM0] + x[G_T]);

    float dQsto2_dQsto1, dQsto2_dQsto2, dQgut_dQsto1, dQgut_dQsto2;
    float Q_sto = x[Q_STO1] + x[Q_STO2];
    float D_bar = y[LAST_Q_STO] + y[D] * 1000;
    if (D_bar > 0) {
        float aa = 5 / (2 * (1 - patient[BETA]) * D_bar);
        float cc = 5 / (2 * patient[DELTA] * D_bar);
        float T1 = tanh(aa * (Q_sto - patient[BETA] * D_bar));
        float T2 = tanh(cc * (Q_sto - patient[DELTA] * D_bar));
        float T3 = x[Q_STO2] * (patient[KMAX] - patient[KMIN]) / 2 * (aa * (1 - T1 * T1) + cc * (1 - T2 * T2));
        dQsto2_dQsto1 = patient[KGRI] - T3;
        dQsto2_dQsto2 = -patient[KMIN] - (patient[KMAX] - patient[KMIN]) / 2 * (T1 - T2 + 2) - T3;
        dQgut_dQsto1 = T3;
        dQgut_dQsto2 = -dQsto2_dQsto2;
    } else {
        dQsto2_dQsto1 = patient[KGRI];
        dQsto2_dQsto2 = -patient[KMAX];
        dQgut_dQsto1 = 0;
        dQgut_dQsto2 = -dQsto2_dQsto2;
    }

    fw(&mA[Q_STO2][Q_STO1], dQsto2_dQsto1);
    fw(&mA[Q_STO2][Q_STO2], dQsto2_dQsto2);
    fw(&mA[Q_GUT][Q_STO1], dQgut_dQsto1);
    fw(&mA[Q_GUT][Q_STO2], dQgut_dQsto2);
    fw(&mA[G_P][G_P], -patient[K1] + dEGPt_dGp - dEt_dGp);
    fw(&mA[G_P][G_P], -patient[K1] + dEGPt_dGp - dEt_dGp);
    fw(&mA[G_P][I_D], dEGPt_dId);
    fw(&mA[G_T][G_T], -patient[K2] - dUidt_dGt);
    fw(&mA[G_T][X], -dUidt_dX);

    fw(&mD[G_P], bias_Gp_EGPt - patient[FCNS] - bias_Gp_Et);
}
