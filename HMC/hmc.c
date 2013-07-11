/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "hmc_utils.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"

#include "cinfo.c"


/* Mesons parameters */
typedef struct _input_mesons {
  char make[256];
  double precision;
  int nhits;
  
  /* for the reading function */
  input_record_t read[4];
  
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"make mesons", "mes:make = %s", STRING_T, (varname).make},\
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"number of noisy sources per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


input_mesons mes_var = init_input_mesons(mes_var);

/* Polyakov-loop parameters */
typedef struct _input_polyakov {
  char make[256];
  
  /* for the reading function */
  input_record_t read[2];
  
} input_polyakov;

#define init_input_polyakov(varname) \
{ \
  .read={\
    {"make polyakov loops", "poly:make = %s", STRING_T, (varname).make},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


input_polyakov poly_var = init_input_polyakov(poly_var);


/* Lowest-eigenvalue parameters */
typedef struct _input_eigval {
  char make[256];
  int nevt; /* search space dimension */
  int nev; /* number of accurate eigenvalues */
  int kmax; /* max degree of polynomial */
  int maxiter; /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
  
  /* for the reading function */
  input_record_t read[8];
  
} input_eigval;

#define init_input_eigval(varname) \
{ \
  .read={\
    {"make lowest eigenvalues", "eva:make = %s", STRING_T, (varname).make},\
    {"search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt},\
    {"number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev},\
    {"max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax},\
    {"max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter},\
    {"absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1},\
    {"relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

input_eigval eigval_var = init_input_eigval(eigval_var);



/* flow control variable */
hmc_flow flow=init_hmc_flow(flow);


char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";
static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: [-i <input file>] [-o <output file>] [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}


double h2evamass=0.;
static void H2eva(spinor_field *out, spinor_field *in){
  g5Dphi_sq(h2evamass, out, in);
}



int main(int argc,char *argv[]) {
  int i, acc, rc;
  char sbuf[128];
 
  read_cmdline(argc,argv);
  
  /* setup process communications */
  setup_process(&argc,&argv);
  
  /* read global variables file */
  read_input(glb_var.read,input_filename);
  
  setup_replicas();
  
  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    freopen(error_filename,"w",stderr);
  }
  
  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS);
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);
  
  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */
  
  /* setup random numbers */
  if(glb_var.rlxd_state[0]!='\0') {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file %s\n",glb_var.rlxd_state);
    read_ranlxd_state(glb_var.rlxd_state);
  } else {
    lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed+MPI_PID);
    rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  }
  
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* read input for measures */
  read_input(glb_var.read,input_filename);
  read_input(mes_var.read,input_filename);
  read_input(poly_var.read,input_filename);
  read_input(eigval_var.read,input_filename);
  
  /* Init Monte Carlo */
  init_mc(&flow, input_filename);
  lprintf("MAIN",0,"MVM during HMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

  
  if(strcmp(mes_var.make,"true")==0) {
    lprintf("OBSERVABLES",0,"Inverter precision for mesons = %e\n",mes_var.precision);
    lprintf("OBSERVABLES",0,"Number of noisy sources for mesons per cnfg = %d\n",mes_var.nhits);
  }
  
  if(strcmp(eigval_var.make,"true")==0) {
    lprintf("OBSERVABLES",0,"EVA Search space dimension  (eva:nevt) = %d\n",eigval_var.nevt);
    lprintf("OBSERVABLES",0,"EVA Number of accurate eigenvalues (eva:nev) = %d\n",eigval_var.nev);
    lprintf("OBSERVABLES",0,"EVA Max degree of polynomial (eva:kmax) = %d\n",eigval_var.kmax);
    lprintf("OBSERVABLES",0,"EVA Max number of subiterations (eva:maxiter) = %d\n",eigval_var.maxiter);
    lprintf("OBSERVABLES",0,"EVA Absolute precision  (eva:omega1) = %e\n",eigval_var.omega1);
    lprintf("OBSERVABLES",0,"EVA Relative precision (eva:omega2) = %e\n",eigval_var.omega2);
  }

  h2evamass=flow.hmc_v->hmc_p.mass;
  
  double *eva_vals=NULL;
  spinor_field *eva_vecs=NULL;
  spinor_field *eva_ws=NULL;
  if(strcmp(eigval_var.make,"true")==0) {
    eva_vals=malloc(sizeof(double)*eigval_var.nevt);
    eva_vecs=alloc_spinor_field_f(eigval_var.nevt+2,&glattice);
    eva_ws=eva_vecs+eigval_var.nevt;
  }
  
  rc=acc=0;
  for(i=flow.start;i<flow.end;++i) {
    int rr;
    double perc;
    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    
    gettimeofday(&start,0);
    
    rr=update_hmc();
    
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Trajectory #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else if(rr!=0) {
      acc++;
    }
    rc++;
    perc=(acc==0)?0.:(float)(100*acc)/(float)(rc);
    
    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM (f;d) = %ld ; %ld\n",i,acc,rc,perc,getMVM_flt(),getMVM());
    
    if((i%flow.save_freq)==0) {
      save_conf(&flow, i);
    }
    
    if((i%flow.meas_freq)==0) {
      /* plaquette */
      lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
      
      /* Mesons */
      if(strcmp(mes_var.make,"true")==0) {
        z2semwall_mesons(i,mes_var.nhits,1,&(flow.hmc_v->hmc_p.mass),mes_var.precision);
      }
      
      /* Polyakov loops */
      if(strcmp(poly_var.make,"true")==0) {
        polyakov();
      }
      
      /* Lowest eigenvalues */
      if(strcmp(eigval_var.make,"true")==0) {
        double max;
        max_H(&H2eva, &glattice, &max);
        max*=1.1;
        
        int status;
        int ie=eva(eigval_var.nev, eigval_var.nevt, 0, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1, eigval_var.omega2, &H2eva, eva_ws, eva_vecs, eva_vals, &status);
        while (ie!=0) { /* if failed restart EVA */
          lprintf("MAIN",0,"Restarting EVA!\n");
          ie=eva(eigval_var.nev, eigval_var.nevt, 2, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1, eigval_var.omega2, &H2eva, eva_ws, eva_vecs, eva_vals, &status);
        }
        
        for (int n=0;n<eigval_var.nev;++n) {
          lprintf("LOWEIG",0,"Eig %d = %1.15e\n",n,eva_vals[n]);
        }
      }
    }
  }
  /* save final configuration */
  if(((--i)%flow.save_freq)!=0) {
    save_conf(&flow, i);
  }
  
  /* Only save state if we have a file to save to */
  if(glb_var.rlxd_state[0]!='\0') {
    lprintf("MAIN",0,"Saving rlxd state to file %s\n",glb_var.rlxd_state);
    write_ranlxd_state(glb_var.rlxd_state);
  }
  
  /* finalize Monte Carlo */
  end_mc();
  
  /* close communications */
  finalize_process();
  
  return 0;
  
}