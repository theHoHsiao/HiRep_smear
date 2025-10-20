/***************************************************************************\
* Copyright (c) 2013, Rudy Arthur, Ari Hietanen                             *
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File spectrum.h
* 
* Functions for measuring spectrum
*
*******************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H
#include "suN.h"
#include "inverters.h"
#include <stdio.h>

void measure_spectrum_semwall(int nm, double* m, int nhits,int conf_num, double precision);
void measure_spectrum_discon_semwall(int nm, double* m, int nhits,int conf_num, double precision);
void measure_spectrum_discon_gfwall(int nm, double* m, int conf_num, double precision);
void measure_spectrum_discon_volume(int nm, double* m, int conf_num, double precision, int dil);
void measure_spectrum_gfwall(int nm, double* m, int conf_num, double precision);
void measure_spectrum_pt(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision);
void measure_spectrum_semwall_ext(int nm, double* m, int nhits,int conf_num, double precision);
void measure_spectrum_pt_ext(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision);
void measure_spectrum_semwall_fixedbc(int dt, int nm, double* m, int nhits,int conf_num, double precision);
void measure_spectrum_pt_fixedbc(int tau, int dt, int nm, double* m, int n_mom,int nhits,int conf_num, double precision);
void measure_spectrum_gfwall_fixedbc(int dt, int nm, double* m, int conf_num, double precision);
void measure_formfactor_pt(int ti, int tf, int nm, double* m, int n_mom, int conf_num, double precision);
void measure_formfactor_fixed(int ti, int tf, int dt, int nm, double* m, int n_mom, int conf_num, double precision);
void measure_conserved_formfactor_fixed(int ti, int tf, int dt, int nm, double* m, int n_mom, int conf_num, double precision);

void measure_diquark_semwall_background(int nm, double* m, int nhits,int conf_num, double precision,double Q, int n);
void measure_baryons(double* m,int conf_num, double precision);


/* For measuring spectrum with a four fermion interaction */
void measure_spectrum_ff_semwall(int nm, double* m, int nhits,int conf_num, double precision);
void measure_spectrum_discon_ff_semwall(int nm, double* m, int nhits, int degree_hopping, int nhits_hopping,int conf_num, double precision);
void measure_spectrum_ff_pt(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision);
void measure_spectrum_semwall_ff_ext(int nm, double* m, int nhits,int conf_num, double precision);

/* smearing Measurements*/
void measure_smearing_ss(int t, int x, int y, int z, int nm, double* m, int n_mom, int nhits, int conf_num, double precision, double epsilon_source, int Nsmear_source, double epsilon_sink, int Nsmear_sink, double APE_epsilon, int APE_N);


/* For measuring meson with fundamental fermions in Anti. Rep. */
void measure_spectrum_pt_fund(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision);
void measure_smearing_ss_fund(int t, int x, int y, int z, int nm, double* m, int n_mom, int nhits, int conf_num, double precision, double epsilon_source, int Nsmear_source, double epsilon_sink, int Nsmear_sink, double APE_epsilon, int APE_N);

/* Chimera baryon measurements*/
void measure_chimera_smear(int t, int x, int y, int z, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int n_mom, int nhits, int conf_num, double precision, double e_source_fund, double e_sink_fund, double e_source_anti, double e_sink_anti, int Nsmear_source, int Nsmear_sink, double APE_epsilon, int APE_N, int mes_fund, int mes_anti);

#endif
