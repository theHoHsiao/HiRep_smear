/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File Dphi_fund.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given
* double-precision spinor field for fundamental representation
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "error.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "spinor_field.h"
#include "geometry.h"
#include "communications.h"
#include "memory.h"
#include "clover_tools.h"

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */


/*
 * Init of Dphi_fund
 */

int init_dirac_fund=1;
spinor_field_fund *gtmp_fund=NULL;
spinor_field_fund *etmp_fund=NULL;
spinor_field_fund *otmp_fund=NULL;

static void free_mem() {
    if (gtmp_fund!=NULL) { free_spinor_field_f_fund(gtmp_fund); etmp_fund=NULL; }
    if (etmp_fund!=NULL) { free_spinor_field_f_fund(etmp_fund); etmp_fund=NULL; }
    if (otmp_fund!=NULL) { free_spinor_field_f_fund(otmp_fund); otmp_fund=NULL; }
    init_dirac_fund=1;
}

void init_Dirac_fund() {
    if (init_dirac_fund) {
        gtmp_fund=alloc_spinor_field_f_fund(1,&glattice);
        etmp_fund=alloc_spinor_field_f_fund(1,&glat_even);
        otmp_fund=alloc_spinor_field_f_fund(1,&glat_odd);
        atexit(&free_mem);
        init_dirac_fund=0;
    }
}


/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_fund_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM_fund() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	MVMcounter=0; /* reset counter */

	return res;
}


/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r,u,s)\
    _suNg_multiply(vtmp,(u),(s));\
    _vector_mulc_g((r),eitheta[0],vtmp)

#define _suNf_theta_T_inverse_multiply(r,u,s)\
    _suNg_inverse_multiply(vtmp,(u),(s));\
    _vector_mulc_star_g((r),eitheta[0],vtmp)

#else

#define _suNf_theta_T_multiply(r,u,s) _suNg_multiply((r),(u),(s))
#define _suNf_theta_T_inverse_multiply(r,u,s) _suNg_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r,u,s)\
_suNg_multiply(vtmp,(u),(s));\
_vector_mulc_g((r),eitheta[1],vtmp)

#define _suNf_theta_X_inverse_multiply(r,u,s)\
_suNg_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_g((r),eitheta[1],vtmp)

#else

#define _suNf_theta_X_multiply(r,u,s) _suNg_multiply((r),(u),(s))
#define _suNf_theta_X_inverse_multiply(r,u,s) _suNg_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r,u,s)\
_suNg_multiply(vtmp,(u),(s));\
_vector_mulc_g((r),eitheta[2],vtmp)

#define _suNf_theta_Y_inverse_multiply(r,u,s)\
_suNg_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_g((r),eitheta[2],vtmp)

#else

#define _suNf_theta_Y_multiply(r,u,s) _suNg_multiply((r),(u),(s))
#define _suNf_theta_Y_inverse_multiply(r,u,s) _suNg_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r,u,s)\
_suNg_multiply(vtmp,(u),(s));\
_vector_mulc_g((r),eitheta[3],vtmp)

#define _suNf_theta_Z_inverse_multiply(r,u,s)\
_suNg_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_g((r),eitheta[3],vtmp)

#else

#define _suNf_theta_Z_multiply(r,u,s) _suNg_multiply((r),(u),(s))
#define _suNf_theta_Z_inverse_multiply(r,u,s) _suNg_inverse_multiply((r),(u),(s))

#endif




/*
 * This function defines the massless Dirac operator
 * It can act on spinors defined on the whole lattice
 * or on spinors with definite parity
 */

void Dphi_fund_(spinor_field_fund *out, spinor_field_fund *in)
{

   error((in==NULL)||(out==NULL),1,"Dphi_fund_ [Dphi_fund.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"Dphi_fund_ [Dphi_fund.c]",
         "Input and output fields must be different");

#ifndef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd,1,"Dphi_ [Dphi.c]", "Spinors don't match! (1)");
   error(out->type==&glat_odd && in->type!=&glat_even,1,"Dphi_ [Dphi.c]", "Spinors don't match! (2)");
   error(out->type==&glattice && in->type!=&glattice,1,"Dphi_ [Dphi.c]", "Spinors don't match! (3)");
#endif


//   lprintf("Dphi_in", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(in),*(in+1),*(in+2),*(in+3),*(in+4),*(in+5),*(in+6),*(in+7));
   ++MVMcounter; /* count matrix calls */
   if(out->type==&glattice) ++MVMcounter;

/************************ loop over all lattice sites *************************/
   /* start communication of input spinor field */
   start_sf_sendrecv_fund(in);

  _PIECE_FOR(out->type,ixp) {
     if(ixp==out->type->inner_master_pieces) {
       _OMP_PRAGMA( master )
       /* wait for spinor to be transfered */
       complete_sf_sendrecv_fund(in);
       _OMP_PRAGMA( barrier )
     }
     _SITE_FOR(out->type,ixp,ix) {

       int iy;
       suNg *up,*um;/*JW suNgfull?*/
       suNg_vector psi,chi;/*JW*/
       suNg_spinor *r=0,*sp,*sm;/*JW*/
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
       suNg_vector vtmp;/*JW*/
#endif

       r=_FIELD_AT(out,ix);

       /******************************* direction +0 *********************************/

//       lprintf("Dphi_in", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(in),*(in+1),*(in+2),*(in+3),*(in+4),*(in+5),*(in+6),*(in+7));
       iy=iup(ix,0);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_fund(ix,0);/*JW*/
//       up=pu_gauge_f_fund(ix,0);/*JW*/
//       lprintf("up",10,"%1.6f\n",(*up).c[0].re);

       _vector_add_g(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[0]=chi;
       (*r).c[2]=chi;

       _vector_add_g(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_T_multiply(chi,(*up),psi);

       (*r).c[1]=chi;
       (*r).c[3]=chi;

       /******************************* direction -0 *********************************/

       iy=idn(ix,0);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_fund(iy,0);/*JW*/

       _vector_sub_g(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_sub_assign_g((*r).c[2],chi);

       _vector_sub_g(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_sub_assign_g((*r).c[3],chi);

       /******************************* direction +1 *********************************/

       iy=iup(ix,1);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_fund(ix,1);/*JW*/

       _vector_i_add_g(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_i_sub_assign_g((*r).c[3],chi);

       _vector_i_add_g(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_X_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_i_sub_assign_g((*r).c[2],chi);

       /******************************* direction -1 *********************************/

       iy=idn(ix,1);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_fund(iy,1);/*JW*/

       _vector_i_sub_g(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_i_add_assign_g((*r).c[3],chi);

       _vector_i_sub_g(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_i_add_assign_g((*r).c[2],chi);

       /******************************* direction +2 *********************************/

       iy=iup(ix,2);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_fund(ix,2);

       _vector_add_g(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_add_assign_g((*r).c[3],chi);

       _vector_sub_g(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_Y_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_sub_assign_g((*r).c[2],chi);

       /******************************* direction -2 *********************************/

       iy=idn(ix,2);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_fund(iy,2);

       _vector_sub_g(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_sub_assign_g((*r).c[3],chi);

       _vector_add_g(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_add_assign_g((*r).c[2],chi);

       /******************************* direction +3 *********************************/

       iy=iup(ix,3);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_fund(ix,3);

       _vector_i_add_g(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_i_sub_assign_g((*r).c[2],chi);

       _vector_i_sub_g(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_i_add_assign_g((*r).c[3],chi);

       /******************************* direction -3 *********************************/

       iy=idn(ix,3);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f_fund(iy,3);

       _vector_i_sub_g(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[0],chi);
       _vector_i_add_assign_g((*r).c[2],chi);

       _vector_i_add_g(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_g((*r).c[1],chi);
       _vector_i_sub_assign_g((*r).c[3],chi);

       /******************************** end of loop *********************************/

       _spinor_mul_g(*r,-0.5,*r);

     } /* SITE_FOR */
   } /* PIECE FOR */
}



/*
 * this function takes 2 spinors defined on the whole lattice
 */
void Dphi_fund(double m0, spinor_field_fund *out, spinor_field_fund *in)
{
   double rho;
#ifdef ROTATED_SF
   int ix,iy,iz,index;
   suNg_spinor *r, *sp;/*JW*/
   double SFrho;
   suNg_spinor tmp;/*JW*/
#endif /* ROTATED_SF */

   error((in==NULL)||(out==NULL),1,"Dphi_fund [Dphi_fund.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"Dphi_fund [Dphi_fund.c]",
         "Input and output fields must be different");

   apply_BCs_on_spinor_field_fund(in);

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

   Dphi_fund_(out, in);

   rho=4.+m0;
   spinor_field_mul_add_assign_f_fund(out,rho,in);/*JW where do I find this function?*/

#ifdef ROTATED_SF
   SFrho=3.*_update_par.SF_ds+_update_par.SF_zf-4.;

	if(COORD[0] == 0) {
		for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
			index=ipt(1,ix,iy,iz);
			r=_FIELD_AT(out,index);
            sp=_FIELD_AT(in,index);
			_spinor_mul_add_assign_g(*r,SFrho,*sp);

			_spinor_pminus_g(tmp,*sp);
			_spinor_g5_assign_g(tmp);
			if(_update_par.SF_sign==1) {
				_spinor_i_add_assign_g(*r,tmp);
			} else {
				_spinor_i_sub_assign_g(*r,tmp);
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
			index=ipt(T-1,ix,iy,iz);
			r=_FIELD_AT(out,index);
            sp=_FIELD_AT(in,index);
			_spinor_mul_add_assign_g(*r,SFrho,*sp);

			_spinor_pplus_g(tmp,*sp);
			_spinor_g5_assign_g(tmp);
			if(_update_par.SF_sign==1) {
				_spinor_i_add_assign_g(*r,tmp);
			} else {
				_spinor_i_sub_assign_g(*r,tmp);
			}
		}
	}
#endif /* ROTATED_SF */

  apply_BCs_on_spinor_field_fund(out);
}

void g5Dphi_fund(double m0, spinor_field_fund *out, spinor_field_fund *in)
{
   double rho;
#ifdef ROTATED_SF
   int ix,iy,iz,index;
   suNg_spinor *r, *sp;
   double SFrho;
   suNg_spinor tmp;
#endif /* ROTATED_SF */

   error((in==NULL)||(out==NULL),1,"g5Dphi_fund [Dphi_fund.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"g5Dphi_fund [Dphi_fund.c]",
         "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"g5Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

   apply_BCs_on_spinor_field_fund(in);

   Dphi_fund_(out, in);

   rho=4.+m0;
   spinor_field_mul_add_assign_f_fund(out,rho,in);/*JW Where do I find this function?*/

#ifdef ROTATED_SF
   SFrho=3.*_update_par.SF_ds+_update_par.SF_zf-4.;

   if(COORD[0] == 0) {
     for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
       index=ipt(1,ix,iy,iz);
       r=_FIELD_AT(out,index);
       sp=_FIELD_AT(in,index);
       _spinor_mul_add_assign_g(*r,SFrho,*sp);

       _spinor_pminus_g(tmp,*sp);
       _spinor_g5_assign_g(tmp);
       if(_update_par.SF_sign==1) {
           _spinor_i_add_assign_g(*r,tmp);
       } else {
	 _spinor_i_sub_assign_g(*r,tmp);
       }
     }
   }
   if(COORD[0] == NP_T-1) {
     for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
       index=ipt(T-1,ix,iy,iz);
       r=_FIELD_AT(out,index);
       sp=_FIELD_AT(in,index);
       _spinor_mul_add_assign_g(*r,SFrho,*sp);

       _spinor_pplus_g(tmp,*sp);
       _spinor_g5_assign_g(tmp);
       if(_update_par.SF_sign==1) {
	 _spinor_i_add_assign_g(*r,tmp);
       } else {
	 _spinor_i_sub_assign_g(*r,tmp);
       }
     }
   }
#endif /* ROTATED_SF */

   spinor_field_g5_assign_f_fund(out);/*JW Where do I find this function?*/

   apply_BCs_on_spinor_field_fund(out);
}



/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the even lattice
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre_fund(double m0, spinor_field_fund *out, spinor_field_fund *in)
{
//  lprintf("Dphi pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(in),*(in+1),*(in+2),*(in+3),*(in+4),*(in+5),*(in+6),*(in+7));
  double rho;

  error((in==NULL)||(out==NULL),1,"Dphi_eopre_fund [Dphi_fund.c]",
	"Attempt to access unallocated memory space");

  error(in==out,1,"Dphi_eopre_fund [Dphi_fund.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"Dphi_eopre_fund " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field_fund(in);

  /* alloc memory for temporary spinor field */
  if (init_dirac_fund) { init_Dirac_fund(); init_dirac_fund=0; }

  Dphi_fund_(otmp_fund, in);
  apply_BCs_on_spinor_field_fund(otmp_fund);
  Dphi_fund_(out, otmp_fund);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_fund(out,rho,in);/*JW Where do I find this function?*/
  spinor_field_minus_f_fund(out,out);/*JW Where do I find this function?*/
  apply_BCs_on_spinor_field_fund(out);
//  lprintf("Dphi pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(out),*(out+1),*(out+2),*(out+3),*(out+4),*(out+5),*(out+6),*(out+7));
}


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */

void Dphi_oepre_fund(double m0, spinor_field_fund *out, spinor_field_fund *in)
{
  double rho;

  error((in==NULL)||(out==NULL),1,"Dphi_oepre_fund [Dphi_fund.c]",
	"Attempt to access unallocated memory space");

  error(in==out,1,"Dphi_oepre_fund [Dphi_fund.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_odd || in->type!=&glat_odd,1,"Dphi_oepre_fund " __FILE__, "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field_fund(in);
  /* alloc memory for temporary spinor field */
  if (init_dirac_fund) { init_Dirac_fund(); init_dirac_fund=0; }

  Dphi_fund_(etmp_fund, in);
  apply_BCs_on_spinor_field_fund(etmp_fund);
  Dphi_fund_(out, etmp_fund);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_fund(out,rho,in);/*JW Where do I find this function?*/
  spinor_field_minus_f_fund(out,out);/*JW Where do I find this function?*/

  apply_BCs_on_spinor_field_fund(out);
}



void g5Dphi_eopre_fund(double m0, spinor_field_fund *out, spinor_field_fund *in)
{
//  lprintf("Dphi pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(in),*(in+1),*(in+2),*(in+3),*(in+4),*(in+5),*(in+6),*(in+7));
  double rho;

  error((in==NULL)||(out==NULL),1,"g5Dphi_eopre_fund [Dphi_fund.c]",
	"Attempt to access unallocated memory space");

  error(in==out,1,"Dphi_eopre_fund [Dphi_fund.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"g5Dphi_eopre_fund " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field_fund(in);

//  lprintf("PF",10,"%1.6f\n",(*in).c[0].c[0]re);
  /* alloc memory for temporary spinor field */
  if (init_dirac_fund) { init_Dirac_fund(); init_dirac_fund=0; }

//  lprintf("Dphi pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(in),*(in+1),*(in+2),*(in+3),*(in+4),*(in+5),*(in+6),*(in+7));
  Dphi_fund_(otmp_fund, in);
  apply_BCs_on_spinor_field_fund(otmp_fund);
  Dphi_fund_(out, otmp_fund);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_fund(out,rho,in);/*JW Where do I find this function?*/
  spinor_field_minus_f_fund(out,out);/*JW Where do I find this function?*/
  spinor_field_g5_assign_f_fund(out);/*JW Where do I find this function?*/

  apply_BCs_on_spinor_field_fund(out);
//  lprintf("Dphi pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(out),*(out+1),*(out+2),*(out+3),*(out+4),*(out+5),*(out+6),*(out+7));
}

/* g5Dphi_eopre ^2 */
void g5Dphi_eopre_sq_fund(double m0, spinor_field_fund *out, spinor_field_fund *in) {
  /* alloc memory for temporary spinor field */
  if (init_dirac_fund) { init_Dirac_fund(); init_dirac_fund=0; }

  g5Dphi_eopre_fund(m0, etmp_fund, in);
  g5Dphi_eopre_fund(m0, out, etmp_fund);

}

/* g5Dhi ^2 */
void g5Dphi_sq_fund(double m0, spinor_field_fund *out, spinor_field_fund *in) {
  /* alloc memory for temporary spinor field */
  if (init_dirac_fund) { init_Dirac_fund(); init_dirac_fund=0; }

#ifdef ROTATED_SF
  /*the switch of the SF_sign is needed to take care of the antihermiticity of the boundary term of the dirac operator*/
  g5Dphi_fund(m0, gtmp_fund, in);
  _update_par.SF_sign = -_update_par.SF_sign;
  g5Dphi_fund(m0, out, gtmp_fund);
  _update_par.SF_sign = -_update_par.SF_sign;
#else
  g5Dphi_fund(m0, gtmp_fund, in);
  g5Dphi_fund(m0, out, gtmp_fund);
#endif

}

//Twisted mass operator for even odd preconditioned case - which is not implemented for multireps
/************************************************* - which is not implemented for multireps
 * Dirac operators with clover term:             *
 * Cphi = Dphi + clover                          *
 * Cphi_eopre = D_ee - D_eo D_oo^-1 D_oe         *
 * Cphi_diag = D_oo or D_ee                      *
 * Cphi_diag_inv = D_oo^-1 or D_ee^-1            *
 *************************************************/

