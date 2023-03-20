#ifndef __FOC_H
#define __FOC_H

#include <fixedptc.h>

/* Constants */
#define FIXEDPT_ONE_SQRTTHREE   fixedpt_rconst(0.577350269)     // 1/sqrt(3)

#define FIXEDPT_TS              fixedpt_fromint(1024)           // Virtual PWM period

/* Datatype */
typedef struct Clarke_Data_st
{
    fixedpt alpha;
    fixedpt beta;
}Clarke_Data;

typedef struct Park_Data_st
{
    fixedpt q;
    fixedpt d;
}Park_Data;

typedef struct SVM_Vector_st
{
    fixedpt t_u;
    fixedpt t_v;
    fixedpt t_w;
    uint32_t sector;
}SVM_Vector;

/* Public Functions */
void ClarkeTransform( Clarke_Data* res, fixedpt i_u, fixedpt i_v, fixedpt i_w );
void ParkTransform( Park_Data* res, Clarke_Data* data, fixedpt phi_angle );
void InvParkTransform( Clarke_Data* res, Park_Data* data, fixedpt phi_angle );
void SVModulation( SVM_Vector* res, Clarke_Data* data, fixedpt V_Motor, const fixedpt PWM_period );

#endif