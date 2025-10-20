//
//  chimera_measurements.c
//  hr
//
//  Created by Ho Hsiao on 2021/4/14.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "spectrum.h"
#include "gaugefix.h"
#include "meson_observables.h"



void measure_chimera_smear(int t, int x, int y, int z, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int n_mom, int nhits, int conf_num, double precision, double e_source_fund, double e_sink_fund, double e_source_anti, double e_sink_anti, int Nsmear_source, int Nsmear_sink, double APE_epsilon, int APE_N, int mes_fund, int mes_anti){
    
    char label[128];
    
    if (APE_N != 0){
    represent_gauge_field_APE();
    }
    
    // calculate fundamental fermion propagator
    spinor_field_fund* source_fund = alloc_spinor_field_f_fund(4,&glattice);
    spinor_field_fund* prop_fund = alloc_spinor_field_f_fund(4*nm_fund*NG, &glattice);
    init_propagator_eo_fund(nm_fund, m_fund, precision);
    
    // calculate antisymmetry fermion propagator
    spinor_field* source_anti = alloc_spinor_field_f(4,&glattice);
    spinor_field* prop_anti = alloc_spinor_field_f(4*nm_anti*NF, &glattice);
    init_propagator_eo(nm_anti, m_anti, precision);
    
    
    int k;
    for (k=0;k<NG;++k){
        
        if (APE_N != 0){
            create_smeared_source_with_APE_fund(source_fund, t, x, y, z, k, e_source_fund, Nsmear_source);
        }
        else{
            create_smeared_source_fund(source_fund, t, x, y, z, k, e_source_fund, Nsmear_source);
        }
        
        calc_propagator_fund(prop_fund + 4*k*nm_fund, source_fund , 4);//4 for spin components
        if (mes_fund){measure_mesons_fund(meson_correlators, prop_fund+4*k*nm_fund, source_fund, nm_fund, t);}
    }
    
    //calc_propagator_fund(prop_fund, source_fund , 4*NG);//4 for spin components
    
    if (mes_fund){
    lprintf("CHIM",0,"Fundamentl meson \n");
    sprintf(label, "source_N%d_sink_N%d_fund",Nsmear_source, 0);
    print_mesons(meson_correlators,1.,conf_num,nm_fund,m_fund,GLB_T,n_mom,label);
    }
    
    for (k=0;k<NF;++k){
        
        if (APE_N != 0){
            create_smeared_source_with_APE(source_anti, t, x, y, z, k, e_source_anti, Nsmear_source);
        }
        else{
            create_smeared_source(source_anti, t, x, y, z, k, e_source_anti, Nsmear_source);
        }
        
        calc_propagator(prop_anti + 4*k*nm_anti, source_anti , 4);//4 for spin components
        
        if (mes_anti){measure_mesons(meson_correlators,prop_anti+4*k*nm_anti, source_anti, nm_anti, t);}
    }
    if (mes_anti){
    lprintf("CHIMERA",0,"Antisymmetry meson \n");
    sprintf(label, "source_N%d_sink_N%d_anti",Nsmear_source, 0);
    print_mesons(meson_correlators,1.,conf_num,nm_anti,m_anti,GLB_T,n_mom,label);
    }
   
    sprintf(label, "source_N%d_sink_N%d", Nsmear_source, 0);
    
    contract_chimera_O(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    contract_chimera_OP(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    contract_chimera_OV(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    //contract_chimera_OV_Ps(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    //contract_chimera_OV_noP(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    //contract_chimera_OV_JW(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
    
    int N_diff = 10; // measure sink smearing every 10 steps
    for (int N = 0; N<Nsmear_sink; N += N_diff){
        lprintf("SMEAR",0,"sink smearing with APE steps: ");
        for(int steps=0;steps<N_diff;steps++){
            lprintf("SMEAR",0,"%d...", steps+N+1);
            
            if (APE_N != 0){
                smeared_propagator_with_APE_fund(prop_fund, nm_fund, e_sink_fund);
                smeared_propagator_with_APE(prop_anti, nm_anti, e_sink_anti);
            }
            
            else{
                smeared_propagator_fund(prop_fund, nm_fund, e_sink_fund);
                smeared_propagator(prop_anti, nm_anti, e_sink_anti);
            }
            
        }
        if (mes_fund){
        for (int c=0;c<NG;++c){
             if (n_mom>1){
                 lprintf("MAIN",0,"no momenta yet\n");
                 //measure_point_mesons_momenta(meson_correlators,prop+4*c, source, nm, t, n_mom);
             }
             else{
                measure_mesons_fund(meson_correlators,prop_fund+4*c, source_fund, nm_fund, t);
             }
        }
        sprintf(label, "source_N%d_sink_N%d_fund",Nsmear_source, N+N_diff);
        print_mesons(meson_correlators,1.,conf_num,nm_fund,m_fund,GLB_T,n_mom,label);
        }
        
        if (mes_anti){
        for (int c=0;c<NF;++c){
             if (n_mom>1){
                 lprintf("MAIN",0,"no momenta yet\n");

                 //measure_point_mesons_momenta(meson_correlators,prop+4*c, source, nm, t, n_mom);
             }
             else{
                measure_mesons(meson_correlators,prop_anti+4*c, source_anti, nm_anti, t);
             }
        }
        sprintf(label, "source_N%d_sink_N%d_anti",Nsmear_source, N+N_diff);
        print_mesons(meson_correlators,1.,conf_num,nm_anti,m_anti,GLB_T,n_mom,label);
        }
        
        sprintf(label, "source_N%d_sink_N%d", Nsmear_source, N+N_diff);
        
        contract_chimera_O(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
        contract_chimera_OP(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
        contract_chimera_OV(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
        //contract_chimera_OV_Ps(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);
        //contract_chimera_OV_noP(prop_fund, prop_anti, nm_fund, nm_anti, m_fund, m_anti,t, conf_num, label);

    }
    
    free_propagator_eo();
    free_propagator_eo_fund();
    free_spinor_field_f(source_anti);
    free_spinor_field_f_fund(source_fund);
    free_spinor_field_f(prop_anti);
    free_spinor_field_f_fund(prop_fund);
    
}

