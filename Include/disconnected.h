
/*******************************************************************************
*
* File observables.h
*
* Functions for measuring observables
*
*******************************************************************************/

#ifndef DISCONNECTED_H
#define DISCONNECTED_H

#include "data_storage.h"
#include "spinor_field.h"
#include "suN.h"
#include <stdio.h>

void measure_bilinear_loops_spinorfield(spinor_field* prop,spinor_field* source,int src_id, int n_mom, storage_switch swc, data_storage_array **ret);
void measure_bilinear_loops_4spinorfield(spinor_field* prop,spinor_field* source,int src_id,int tau,int col,int eo,storage_switch swc, data_storage_array **ret);

#ifndef GAUGE_SPN // FIXFORSPN
void measure_loops(double* m, int nhits,int conf_num, double precision,int source_type,int n_mom,storage_switch swc, data_storage_array **ret);
#endif

#endif
