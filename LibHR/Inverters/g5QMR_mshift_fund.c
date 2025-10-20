//
//  g5QMR_mshift_fund.c
//  hr
//
//  Created by Paul Xiao on 2021/4/7.
//  Copyright © 2021 Claudio Pica. All rights reserved.
//

/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *
* All rights reserved.                                                      *
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "memory.h"
#include "update.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

// Uncomment to use BiCGstab instead of MINRES (does not work with multishift!)
//#define USE_BICGSTAB

/* DA ELIMINARE */
#include "observables.h"
#undef NDEBUG


#define _print_par(c)                \
printf("[%d] " #c " = %e\n",cgiter,c)

/* Calculate the Givens rotation (c,s) that maps:
* |  c   s | | a | = | 1/r |
* | -s^+ c | | b |   | 0 |
*/
#define _Givens_rot(r,c,s,a,b)                        \
{ if (fabs(a)<1.e-15) { (r)=1./(b); (c)=0.; (s)=1.; }            \
  else { (s)=(b)/(a); (c)=1./sqrt(1.+(s)*(s)); (s)*=(c); (r)=1./((c)*(a)+(s)*(b)); } \
}



/*
* performs the multi-shifted QMR inversion for g5-hermitean matrices:
* out[i] = (M-(par->shift[i]))^-1 in
* returns the number of cg iterations done.
*/


static int g5QMR_mshift_core_fund(short *valid, mshift_par *par, spinor_operator_fund M, spinor_field_fund *in, spinor_field_fund *out){
  
  spinor_field_fund **q1,**q2;
  spinor_field_fund *p1, *p2, *Mp;
  spinor_field_fund *sptmp, *memall;
  
  double alpha, beta, delta, rho, innorm2;
  double *r, *s1, *s2, *c1, *c2;
  double maxm;
  
  int i;
  int cgiter;
  unsigned int notconverged;
  
  unsigned short *flags;
  
  /* fare qualche check sugli input */
  /* par->n deve essere almeno 1! */
  assert(par->n>0);
  _TWO_SPINORS_MATCHING(in,&out[0]);
  _ARRAY_SPINOR_MATCHING(out,par->n)
  
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
  * objects of the same type are allocated together
  */
  memall = alloc_spinor_field_f_fund(2*(par->n)+3,in->type);
  q1 = (spinor_field_fund**)malloc(sizeof(spinor_field_fund*)*par->n);
  q2 = (spinor_field_fund**)malloc(sizeof(spinor_field_fund*)*par->n);
  for(i=0; i<par->n; i++) {
    q1[i] = memall+i;
    q2[i] = memall+par->n+i;
  }
  p1 = memall+2*par->n;
  p2 = p1+1;
  Mp = p2+1;
  
  r = (double *)malloc(sizeof(double)*5*(par->n));
  s1 = r+(par->n);
  s2 = s1+(par->n);
  c1 = s2+(par->n);
  c2 = c1+(par->n);
  
  flags=(unsigned short *)malloc(sizeof(unsigned short)*(par->n));
  
  /* init recursion */
  cgiter = 0;
  notconverged=par->n;
  
  spinor_field_copy_f_fund(p2, in); /* trial solution = 0 */
  innorm2=spinor_field_sqnorm_f_fund(in);
  if(par->n==1) { /* in this case is not a multishift and we use as starting vector out[0] */
    M(Mp,&out[0]);
    spinor_field_mul_add_assign_f_fund(Mp,-par->shift[0],&out[0]);
    spinor_field_sub_f_fund(p2,p2,Mp);
  }
  rho=sqrt(spinor_field_sqnorm_f_fund(p2));
  lprintf("INVERTER",20,"g5QMR: rho init: %1.8e\n",rho*rho/innorm2);
  
  spinor_field_mul_f_fund(p2,1./rho,p2);
  spinor_field_zero_f_fund(p1);
  for (i=0; i<(par->n); ++i) {
    r[i]=rho;
    c2[i]=c1[i]=1.;
    s1[i]=s2[i]=0.;
    if (par->n!=1) spinor_field_zero_f_fund(&out[i]); /* if no multishift we start with the trial solution */
    spinor_field_zero_f_fund(q1[i]);
    spinor_field_zero_f_fund(q2[i]);
    flags[i]=1;
  }
  delta = spinor_field_g5_prod_re_f_fund(p2,p2);
  beta = delta;
  
  /* cg recursion */
  do {
    ++cgiter;
    
    M(Mp,p2);
    spinor_field_mul_add_assign_f_fund(Mp,-par->shift[0],p2);
    
    /* compute alpha */
    alpha = spinor_field_g5_prod_re_f_fund(p2,Mp)/delta;
    
    /* update p1, p2 */
    spinor_field_mul_add_assign_f_fund(Mp,-beta,p1);
    spinor_field_mul_add_assign_f_fund(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;
    
    /* update rho */
    rho=sqrt(spinor_field_sqnorm_f_fund(p2));
    
    maxm=1.e-10; /* to check if the error is going down */
    
    for (i=0; i<(par->n); ++i) { /* update solutions */
      if (flags[i]) {
        double a, t, e, d, m;
        a=(alpha-par->shift[i]+par->shift[0]);
        t=s1[i]*beta;
        e=c2[i]*c1[i]*beta+s2[i]*a;
        m=c2[i]*a-s2[i]*c1[i]*beta;
        s1[i]=s2[i];
        c1[i]=c2[i];
        _Givens_rot(d,c2[i],s2[i],m,rho);
        
        maxm=(maxm>fabs(m))?maxm:fabs(m);
        
        /* update q */
        spinor_field_lc_f_fund(Mp,-t*d,q1[i],-e*d,q2[i]);
        spinor_field_mul_add_assign_f_fund(Mp,d,p1);
        sptmp=q1[i];
        q1[i]=q2[i];
        q2[i]=Mp;
        Mp=sptmp;/* swap q1[i]<-q2[i]<-Mp and Mp point to q1 */
        
        /* update solution */
        spinor_field_mul_add_assign_f_fund(&out[i],c2[i]*r[i],q2[i]);
        
        /* update residuum */
        r[i]*=-s2[i];
        
        if(fabs(r[i]*r[i])*(double)(1+cgiter)<par->err2*innorm2){
          flags[i]=0;
          --notconverged;
          lprintf("INVERTER",50,"g5QMR: vect %d converged at iter: %d\n",i,cgiter);
        } else {
          lprintf("INVERTER",50,"g5QMR: vect %d iter: %d res: %1.8e s2: %1.8e m=%1.8e\n",i,cgiter,fabs(r[i]*r[i])*(double)(1+cgiter)/innorm2,s2[i],m);
        }
        
      }
    }
    
    if (fabs(rho)> 1.e-13) {
      double olddelta;
      
      /*normalize p2 */
      spinor_field_mul_f_fund(p2,1./rho,p2);
      
      /* update delta and beta */
      olddelta=delta;
      delta = spinor_field_g5_prod_re_f_fund(p2,p2);
      
      lprintf("INVERTER",40,"g5QMR: delta=%1.8e rho=%1.8e [iter=%d]\n",delta,rho,cgiter);
      if(fabs(delta)<1.e-6 || maxm<1.e-3) { /* the method has failed ! */
        lprintf("INVERTER",40,"g5QMR: method failed! delta=%e [iter=%d]\n",delta,cgiter);
        notconverged=0; /* exit loop */
      }
      
      /* update beta */
      beta=rho*delta/olddelta;
      
    } else { /* system has converged */
      notconverged=0;
      lprintf("INVERTER",40,"g5QMR: rho < 1.e-13 ! (system converged)\n");
    }
    
    
    #ifndef NDEBUG
    if(cgiter%100==0){
      lprintf("INVERTER",40,"g5QMR: [%d] res[0]=%e notconverged=%d\n",cgiter,fabs(r[0]*r[0])*(double)(1+cgiter)/innorm2,notconverged);
    }
    #endif
    
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);
  
  /* test results */
  for(i=0;i<par->n;++i){
    double norm;
    M(Mp,&out[i]);
    ++cgiter;
    if(par->shift[i]!=0.) {
      spinor_field_mul_add_assign_f_fund(Mp,-par->shift[i],&out[i]);
    }
    spinor_field_sub_f_fund(Mp,Mp,in);
    norm=spinor_field_sqnorm_f_fund(Mp)/innorm2;
    valid[i]=1;
    if (fabs(norm)>par->err2){
      valid[i]=0;
      lprintf("INVERTER",30,"g5QMR failed on vect %d: err2 = %1.8e > %1.8e\n",i,norm,par->err2);
    } else {
      lprintf("INVERTER",20,"g5QMR inversion: err2 = %1.8e < %1.8e\n",norm,par->err2);
    }
  }
  
  /* free memory */
  free_spinor_field_f_fund(memall);
  free(q1);
  free(q2);
  free(r);
  
  free(flags);
  
  /* return number of cg iter */
  if(par->n==1) ++cgiter;
  return cgiter;
}

static double sh_fund;
static spinor_operator_fund g5Herm_fund;

static void Herm_fund(spinor_field_fund *out, spinor_field_fund *in){
  g5Herm_fund(out,in);
  if(sh_fund!=0.) {
    spinor_field_mul_add_assign_f_fund(out,-sh_fund,in);
  }
  spinor_field_g5_f_fund(out,out);
}


int g5QMR_mshift_fund(mshift_par *par, spinor_operator_fund M, spinor_field_fund *in, spinor_field_fund *out)
{
  int cgiter;
  int n;
  mshift_par orig;
  short *valid;
  int loccg;
  int msiter;

  orig=*par; /* save par */
  valid=malloc(sizeof(short)*orig.n);
  cgiter=g5QMR_mshift_core_fund(valid,par,M,in,out);
  msiter=cgiter;

/*#ifdef USE_BICGSTAB
  if(par->n == 1)
  {
      if(valid[0] == 0)
      {
          cgiter += BiCGstab(par,M,in,out);
      }
  }
  else
  {
#endif*/

  /* if some vector has not converged try non-multishift */
  par->n=1;
  for(n=0;n<orig.n;++n){
    if (valid[n]==0) {
      #ifndef NDEBUG
      lprintf("INVERTER",20,"Doing non multishift on vector %d\n",n);
      #endif
      par->shift=orig.shift+n;
      par->max_iter=0;
      loccg=g5QMR_mshift_core_fund(valid+n,par,M,in,out+n);
      cgiter+=loccg;

      if (valid[n]==0) {
        /* revert to MINRES */
        MINRES_par Mpar;

        lprintf("INVERTER",20,"g5QMR_mshift: using MINRES on vect %d\n",n);

        Mpar.err2=par->err2;
        Mpar.max_iter=0;

        g5Herm_fund=M;
        sh_fund=par->shift[0];

        spinor_field_g5_f_fund(in,in); /* multiply input by g5 for MINRES */
        loccg=MINRES_fund(&Mpar,&Herm_fund,in,&out[n],(n==0)?&out[n]:&out[n-1]);
        spinor_field_g5_f_fund(in,in); /* restore input vector */

        cgiter+=loccg;
      }
    }
  }
#ifdef USE_BICGSTAB
  }
#endif

  /* this is for debug purposes */
  *par=orig; /* restore par */
  free(valid);
  lprintf("INVERTER",10,"g5QMR_mshift: cgiter (mshift,tot) = %d ; %d\n",msiter,cgiter);

  return cgiter;
}

