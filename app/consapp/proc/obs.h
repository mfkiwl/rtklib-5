/*
 Author: Dr. Yudan Yi
 Date: 09/20/2021
*/
#ifndef _OBS_H_
#define _OBS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/* data struct used in the engine */
/*
*                       0     1     2     3     4
*           --------------------------------------
*            GPS       L1    L2    L5     -     -
*            GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
*            Galileo   E1    E5b   E5a   E6   E5ab
*            QZSS      L1    L2    L5    L6     -
*            SBAS      L1     -    L5     -     -
*            BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
*            NavIC     L5     S     -     -     -
*/

#ifndef MAX_FRQ 
#define MAX_FRQ 5
#endif

typedef struct
{
	uint8_t sys; /* satellite sys ID G,R,E,C,J,S */
	uint8_t prn; /* satellite prn number/slot number */
	uint8_t sat;
	uint8_t code[MAX_FRQ]; /* measurement code */
	double P[MAX_FRQ]; /* code measurement in meter */
	double L[MAX_FRQ]; /* phase measurement in cycle */
	double D[MAX_FRQ]; /* doppler measurement in HZ */
	uint8_t S[MAX_FRQ]; /* signal SNR/CN0 */
	uint8_t LLI[MAX_FRQ]; /* lost-lock indicator */
	uint32_t lock[MAX_FRQ];
    double wave[MAX_FRQ];
}sat_obs_t;

typedef struct {
    double  rs[6];
    double  dts[2];
    double  var;
    int svh;
    double azel[2];    /*azimuth,elevation*/
    double e[3];       /*partial deviation*/
    double tgd;        /* tgd*/
    double r;          /* vector */
    double tro;        /* tropospheric */
}sat_vec_t;

#ifndef MAX_SAT 
#define MAX_SAT (65)
#endif

typedef struct
{
    uint8_t n;
    sat_obs_t obs[MAX_SAT];
    sat_vec_t vec[MAX_SAT];
    uint16_t wk;
    double time;
    double pos[3];
}epoch_t;


#ifdef __cplusplus
}
#endif

#endif