/*******************************************************************************
*
* Check of the glueballs operators for active against passive rotations.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "update.h"
#include "observables.h"
#include "random.h"
#include "logger.h"
#include "communications.h"
#include "representation.h"
#include "glueballs.h"
#include "setup.h"
#include "hr_complex.h"
#include "check_utils_3_functions.c"

static void all_g_op(double complex *pa)
{
    suNg_field *_u = u_gauge_wrk();
    start_gf_sendrecv(_u);
    int i;

    for (i = 0; i < total_n_glue_op; i++)
        pa[i] = 0.;

    for (i = 0; i < n_active_slices; i++)
        eval_all_glueball_ops(active_slices_list[i], pa);

    for (i = 0; i < total_n_glue_op; i++)
        pa[i] /= n_active_slices * NP_T;

    global_sum((double *)(pa), 2 * total_n_glue_op);
}

int main(int argc, char *argv[])
{

    int return_value = 0;
    int idx_wrk;
    double complex *op, *rop;
    int j, ret;
    int **space_rotations;
    int **inverse_space_rotations;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    initialize_spatial_active_slices(NULL);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    op = malloc(total_n_glue_op * sizeof(complex double));
    rop = malloc(total_n_glue_op * sizeof(complex double));


    lprintf("MAIN", 0, "Measuring all the glueballs operators on the original configuration\n");
    all_g_op(op);

    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Requesting, one workspace gauge field\n");
    idx_wrk = reserve_wrk_space();
    lprintf("MAIN", 0, "done.\n\n");
    space_rotations = direct_spatial_rotations();
    inverse_space_rotations = inverse_spatial_rotations();

    lprintf("MAIN", 0, "Resetting and initializing (rotating) workspace gauge field once for each of the 48 cubic rotations\n\n");

    for (j = 0; j < 48; j++)
    {
        lprintf("MAIN", 0, "Rotation %d  %d->%d   %d->%d   %d->%d   %d->%d\n", j, 0, space_rotations[j][0], 1, space_rotations[j][1], 2, space_rotations[j][2], 3, space_rotations[j][3]);

        assign_spatial_rotated_wrkspace(inverse_space_rotations[j], idx_wrk);

        lprintf("MAIN", 0, "Measuring all the glueballs operators on the rotated configuration\n");
       
        all_g_op(rop);

        ret = fullcheck(j, rop, op);

        if (ret == 0)
            lprintf("MAIN", 0, "active - passive rotation  %d: Pass\n", j);
        else
        {
            lprintf("MAIN", 0, "active - passive rotation  %d: Fail\n", j);
            return_value += ret;
        }
        lprintf("MAIN", 0, "done.\n\n");
    }
    free(op);
    free(rop);
    finalize_process();
    return return_value;
}
