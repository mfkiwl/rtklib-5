/*
 Author: Dr. Yudan Yi
 Date: 09/20/2021
*/
#ifndef _ENGINE_H_
#define _ENGINE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#include "obs.h"

int add_rcv_data(int staid, epoch_t* epoch);
int engine_proc();

#ifdef __cplusplus
}
#endif

#endif