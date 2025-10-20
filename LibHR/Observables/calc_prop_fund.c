//
//  calc_prop_fund.c
//  hr
//
//  Created by Paul Xiao on 2021/4/7.
//  Copyright © 2021 Claudio Pica. All rights reserved.
//

/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen, Jarno Rantaharju            *
*                                                                           *
*                                                                           *
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include "io.h"
#include "random.h"
#include "communications.h"
#include "ranlux.h"

#define PI 3.141592653589793238462643383279502884197

//Helps QMR solver find more accurate solutions
#undef GAUSSIAN_NOISE
//#define GAUSSIAN_NOISE


static double hmass_pre;

static void D_pre_fund(spinor_field_fund *out, spinor_field_fund *in){
#ifdef WITH_CLOVER
  Cphi_eopre(hmass_pre,out,in);
#else
  Dphi_eopre_fund(hmass_pre,out,in);
#endif
}

static void H_pre_fund(spinor_field_fund *out, spinor_field_fund *in){
  g5Dphi_eopre_fund(hmass_pre,out,in);
}

static void H2_pre_fund(spinor_field_fund *out, spinor_field_fund *in){
  g5Dphi_eopre_sq_fund(hmass_pre,out,in);
}



static int init=0;
static int init_odd=0;
static int init_eig=0;
static int neigs=0;

static mshift_par QMR_par;
static double *shift;
static double *mass;
#ifdef GAUSSIAN_NOISE
static spinor_field_fund *QMR_noise;
static spinor_field_fund *QMR_resdn;
#endif

/* PX: for Mix Rep. */
static spinor_field_fund *resd_fund;
static spinor_field_fund *tmp_fund;
static spinor_field_fund *tmp_odd_fund;

// EVA parameters
double *eva_val;
static spinor_field_fund *eva_vec;
static spinor_field_fund *tmp_sf;

enum {_g5QMR=0, _MINRES, _CG, _CG_4F};

/* Initialises the propagator, nm is the number of masses for multimass solver,
   m is the array of masses, and acc is the inverter accuracy
*/

void init_propagator_eo_fund(int nm, double *m, double acc){
  int i;
#ifdef GAUSSIAN_NOISE
  int cgiter=0;
  double norm;
#endif

  if(init==0) {
    shift=(double*)malloc(sizeof(double)*(nm));
    mass=(double*)malloc(sizeof(double)*(nm));
    hmass_pre=m[0]; /* we can put any number here!!! */
    for(i=0;i<nm;++i){
      mass[i]=m[i];
      shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+m[i])*(4.+m[i]);
    }
    QMR_par.n = nm;
    QMR_par.shift = shift;
    QMR_par.err2 = .5*acc;
    QMR_par.max_iter = 0;
 
    resd_fund=alloc_spinor_field_f_fund(QMR_par.n,&glat_even);
    tmp_fund=alloc_spinor_field_f_fund(1,&glat_even);

// not yet for mix rep.
#ifdef GAUSSIAN_NOISE
    QMR_noise=alloc_spinor_field_f(nm+1,&glat_even);
    QMR_resdn=QMR_noise+1;
    /* noisy background */
    gaussian_spinor_field(QMR_noise);
    norm=sqrt(spinor_field_sqnorm_f(QMR_noise));
    spinor_field_mul_f(QMR_noise,1./norm,QMR_noise);
#endif
    init=1;
  }
#ifdef GAUSSIAN_NOISE
  /* invert noise */
  for(i=0;i<QMR_par.n;++i) spinor_field_zero_f(&QMR_resdn[i]);
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, QMR_noise, QMR_resdn);
  lprintf("Z2SEMWALL",10,"QMR_eo MVM = %d\n",cgiter);
#endif
}


void free_propagator_eo_fund() {
  error(init==0,1,"calc_prop.c","propagator not initialized!");

  free_spinor_field_f_fund(tmp_fund);
  free_spinor_field_f_fund(resd_fund);

  free(shift);
  free(mass);

// not yet ready for mix rep.
#ifdef GAUSSIAN_NOISE
  free_spinor_field_f_fund(QMR_noise);
#endif
  if (init_odd){
    free_spinor_field_f_fund(tmp_odd_fund);
    init_odd = 0;
  }
  if (init_eig){
    free(eva_val);
    free_spinor_field_f_fund(eva_vec);
    init_eig = 0;
  }
  init=0;
}

/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/

// To calculate the fundamental fermion propagator


static void calc_propagator_core_fund(spinor_field_fund *psi, spinor_field_fund *eta, int solver) {

  start_sf_sendrecv_fund(eta);
  complete_sf_sendrecv_fund(eta);

  spinor_field_fund qprop_mask_eta;
  spinor_field_fund qprop_mask_psi;
    
  int cgiter=0;
  if (init_odd==0){
    tmp_odd_fund = alloc_spinor_field_f_fund(1,&glat_odd);
    init_odd=1;
  }
  error(init==0,1,"calc_prop_fund.c","calc_propagator_core method not initialized!");

  /* Construct source
     eta_even' = eta_even - D_eo D_oo^-1 eta_odd
   */
  qprop_mask_eta=*eta;
  qprop_mask_eta.type=&glat_odd;
  qprop_mask_eta.ptr=eta->ptr+glat_odd.master_shift;
  spinor_field_mul_f_fund(tmp_odd_fund,(1./( 4.+ mass[0] )),&qprop_mask_eta);
  Dphi_fund_(tmp_fund,tmp_odd_fund);
  qprop_mask_eta.type=&glat_even;
  qprop_mask_eta.ptr=eta->ptr;
  spinor_field_sub_f_fund(tmp_fund,&qprop_mask_eta,tmp_fund);
#ifdef GAUSSIAN_NOISE
  spinor_field_add_assign_f_fund(tmp_fund,QMR_noise_fund);
#endif
  //spinor_field_sub_f(resd,&qprop_mask_eta,tmp);

//if the solution vector is empty use zero guess
if( spinor_field_sqnorm_f_fund(psi) < 1e-28 ){
  spinor_field_zero_f_fund(resd_fund);
} else {
      psi[0].type = &glat_even;
          spinor_field_mul_f_fund(resd_fund,1/(4.+mass[0]),psi);
      psi[0].type = &glattice;
}

#ifdef GAUSSIAN_NOISE
  cgiter+=g5QMR_mshift_fund(&QMR_par, &D_pre, tmp, resd);
#else
  if(solver == _CG){
    qprop_mask_eta.type=&glat_even;
    spinor_field_copy_f_fund( &qprop_mask_eta, tmp_fund );
    spinor_field_g5_assign_f_fund( &qprop_mask_eta );
    H_pre_fund(tmp_fund, &qprop_mask_eta);
    cgiter+=cg_mshift_fund(&QMR_par, &H2_pre_fund, tmp_fund, resd_fund);
  } else if(solver == _MINRES){
    spinor_field_g5_f_fund(tmp_fund,tmp_fund);
    cgiter+=MINRES_mshift_fund(&QMR_par, &H_pre_fund, tmp_fund, resd_fund);
  } else {
      cgiter+=g5QMR_mshift_fund(&QMR_par, &D_pre_fund, tmp_fund, resd_fund);
  }

  //spinor_field_g5_f(tmp,tmp);
  //cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp, resd);

  //cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);

  //cg stuff
  /*qprop_mask_eta.type=&glat_even;
  spinor_field_copy_f( &qprop_mask_eta, tmp );
  spinor_field_g5_assign_f( &qprop_mask_eta );
  H_pre(tmp, &qprop_mask_eta);
  cgiter+=cg_mshift(&QMR_par, &H2_pre, tmp, resd);*/
#endif

#ifdef GAUSSIAN_NOISE
  spinor_field_sub_assign_f_fund(resd,QMR_resdn);
#endif
  /* compute solution
     psi_even = D_ee*resd_e
     psi_odd = D_oo^-1*eta_odd-D_oe resd_e
  */

  qprop_mask_psi=*psi;
  qprop_mask_psi.type=&glat_even;
  spinor_field_mul_f_fund(&qprop_mask_psi,(4.+mass[0]),resd_fund);

  qprop_mask_psi.type=&glat_odd;
  qprop_mask_psi.ptr=psi->ptr+glat_odd.master_shift;
  Dphi_fund_(&qprop_mask_psi,resd_fund);
  
  spinor_field_sub_f_fund(&qprop_mask_psi,tmp_odd_fund,&qprop_mask_psi);

  ++cgiter; /* One whole call*/
  lprintf("CALC_PROP_CORE_fund",10,"QMR_eo MVM = %d\n",cgiter);


   start_sf_sendrecv_fund(psi);
   complete_sf_sendrecv_fund(psi);

}

void calc_propagator_fund(spinor_field_fund *psi, spinor_field_fund* eta, int ndilute){
    int beta,i,n_masses;
    double *m;
    m = mass;
    n_masses=QMR_par.n;
    QMR_par.n=1;
    for (beta=0;beta<ndilute;++beta){
        for (i=0;i<n_masses;++i){
            lprintf("CALC_PROPAGATOR_fund",10,"n masses=%d, mass = %g\n",n_masses, mass[0]);
            hmass_pre = mass[0];
#ifdef WITH_CLOVER
            calc_propagator_clover(&psi[beta*n_masses+i],&eta[beta]);
#else
            calc_propagator_core_fund(&psi[beta*n_masses+i],&eta[beta],_g5QMR);
            //lprintf("CALC_PROPAGATOR_fund",0,"beta*n_masses+i =  %d\n",beta*n_masses+i);
#endif
            mass++;
        }
        mass = m;
    }
    QMR_par.n = n_masses;
    hmass_pre = mass[0];
}
