#include <foc.h>

/* Implementation */
void ClarkeTransform( Clarke_Data* res, fixedpt i_u, fixedpt i_v, fixedpt i_w )
{
    res->alpha = i_u;
    res->beta = fixedpt_mul( FIXEDPT_ONE_SQRTTHREE, fixedpt_add( i_u, fixedpt_mul( FIXEDPT_TWO, i_v ) ) );  // beta = (i_u + 2i_v) / sqrt(3)
}

void ParkTransform( Park_Data* res, Clarke_Data* data, fixedpt phi_angle )
{
    fixedpt cos_phi = fixedpt_cos( phi_angle );
    fixedpt sin_phi = fixedpt_sin( phi_angle );

    res->d = fixedpt_add( fixedpt_mul( data->alpha, cos_phi ), fixedpt_mul( data->beta, sin_phi ) ); // i_d = i_alpha*cos(phi) + i_beta*(sin_phi)
    res->q = fixedpt_sub( fixedpt_mul( data->beta, cos_phi ), fixedpt_mul( data->alpha, sin_phi ) );  // i_q = -i_alpha*sin(phi) + i_beta*(cos_phi)
}

void InvParkTransform( Clarke_Data* res, Park_Data* data, fixedpt phi_angle )
{
    fixedpt cos_phi = fixedpt_cos( phi_angle );
    fixedpt sin_phi = fixedpt_sin( phi_angle );

    res->alpha = fixedpt_sub( fixedpt_mul( data->d, cos_phi ), fixedpt_mul( data->q, sin_phi ) ); // i_alpha = i_d*cos(phi) - i_q*(sin_phi)
    res->beta = fixedpt_sub( fixedpt_mul( data->q, cos_phi ), fixedpt_mul( data->d, sin_phi ) );  // i_beta = -i_d*sin(phi) + i_q*(cos_phi)
}

/* Segment border */
#define FIXEDPT_V1  0
#define FIXEDPT_V2  fixedpt_rconst(3.14159265358979323846 / 3)          // pi / 3
#define FIXEDPT_V3  fixedpt_rconst(2 * 3.14159265358979323846 / 3)      // 2pi / 3
#define FIXEDPT_V4  fixedpt_rconst(3.14159265358979323846)              // 3pi / 3 = pi
#define FIXEDPT_V5  fixedpt_rconst(-2 * 3.14159265358979323846 / 3)     // -2pi / 3
#define FIXEDPT_V6  fixedpt_rconst(-3.14159265358979323846 / 3)         // -pi / 3

#define FIXEDPT_SQRT_THREE fixedpt_rconst(1.73205080757)                // sqrt(3)

void SVModulation( SVM_Vector* res, Clarke_Data* data, fixedpt V_Motor, const fixedpt PWM_period )
{
    fixedpt v_amp, v_angle, v_angle_normalized;
    fixedpt t1, t2, t0, t0_by_2;         // PWM time in each vector
    fixedpt ma;                 // Temp

    /* Cartesian to Polar */
    v_amp = fixedpt_sqrt( fixedpt_add( fixedpt_mul( data->alpha, data->alpha ), fixedpt_mul( data->beta, data->beta ) ) );  // v_amp = sqrt( alpha^2 + beta^2 )
    v_angle = fixedpt_atan2( data->alpha, data->beta );

    /* Normalize angle */
    if( v_angle < FIXEDPT_V5 )
    {
        res->sector = 4;
        v_angle_normalized = fixedpt_add( v_angle, FIXEDPT_V4 );
    }
    else if( v_angle < FIXEDPT_V6 )
    {
        res->sector = 5;
        v_angle_normalized = fixedpt_add( v_angle, FIXEDPT_V3 );
    }
    else if( v_angle < FIXEDPT_V1 )
    {
        res->sector = 6;
        v_angle_normalized = fixedpt_add( v_angle, FIXEDPT_V2 );
    }
    else if( v_angle < FIXEDPT_V2 )
    {
        res->sector = 1;
        v_angle_normalized = v_angle;
    }
    else if( v_angle < FIXEDPT_V3 )
    {
        res->sector = 2;
        v_angle_normalized = fixedpt_sub( v_angle, FIXEDPT_V2 );
    }
    else
    {
        res->sector = 3;
        v_angle_normalized = fixedpt_sub( v_angle, FIXEDPT_V3 );
    }

    /* Calculate time for each vector */
    ma = fixedpt_div( fixedpt_mul( FIXEDPT_SQRT_THREE, fixedpt_mul( PWM_period, v_amp ) ), V_Motor );   // ma = (sqrt(3) * PWM_period * v_amp) / V_Motor
    t1 = fixedpt_mul( ma, fixedpt_sin( fixedpt_sub( FIXEDPT_V2, v_angle_normalized ) ) );               // t1 = ma * sin( pi/3 - v_angle_normalized )
    t2 = fixedpt_mul( ma, fixedpt_sin( v_angle_normalized ) );                                          // t2 = ma * sin( v_angle_normalized )
    t0 = fixedpt_sub( PWM_period, fixedpt_add( t1, t2 ) );                                              // t0 = PWM_period - (t1 + t2)

    /* PWM mapping */
    t0_by_2 = fixedpt_div( t0, fixedpt_fromint(2) ) ;           // t0 / 2
    switch(res->sector)
    {
        case 1:
        res->t_u = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_v = fixedpt_add( t2, t0_by_2 );
        res->t_w = t0_by_2;
        break;

        case 2:
        res->t_v = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_u = fixedpt_add( t2, t0_by_2 );
        res->t_w = t0_by_2;
        break;

        case 3:
        res->t_v = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_w = fixedpt_add( t2, t0_by_2 );
        res->t_u = t0_by_2;
        break;

        case 4:
        res->t_w = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_v = fixedpt_add( t2, t0_by_2 );
        res->t_u = t0_by_2;
        break;

        case 5:
        res->t_w = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_u = fixedpt_add( t2, t0_by_2 );
        res->t_v = t0_by_2;
        break;

        default:
        res->t_u = fixedpt_add( t1, fixedpt_add( t2, t0_by_2 ) );
        res->t_w = fixedpt_add( t2, t0_by_2 );
        res->t_v = t0_by_2;
        break;
    }
}
