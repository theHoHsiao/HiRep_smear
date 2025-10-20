//
//  chimera_baryon.c
//  hr
//
//  Created by Ho Hsiao on 2021/4/15 .
//
// To calculate chimera baryons in Sp(2N) gauge theory

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
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"
#include <string.h>
#include "meson_observables.h"
#include "complex.h"
#define PI 3.141592653589793238462643383279502884197

static void op_propagator(suNf_propagator* out, suNf_propagator* in, int i){
  switch (i){
  case 5: _g5_propagator( (*out),(*in) ); break;
  case 0: _g0_propagator( (*out),(*in) ); break;
  case 1: _g1_propagator( (*out),(*in) ); break;
  case 2: _g2_propagator( (*out),(*in) ); break;
  case 3: _g3_propagator( (*out),(*in) ); break;
  default: break;
  }
}

static void op_propagator_fund(suNg_propagator* out, suNg_propagator* in, int i){
  switch (i){
  case 5: _g5_propagator_fund( (*out),(*in) ); break;
  case 0: _g0_propagator_fund( (*out),(*in) ); break;
  case 1: _g1_propagator_fund( (*out),(*in) ); break;
  case 2: _g2_propagator_fund( (*out),(*in) ); break;
  case 3: _g3_propagator_fund( (*out),(*in) ); break;
  default: break;
  }
}


void contract_chimera_O(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    char channel_name[][5] = {"OA", "OC", "OAC", "OCA"};
    char label_plus[256];
    int ix, tc, idx;
    complex dqC, dqA, dqCA, dqAC, tmp;

    suNg Omega, Otr;
    suNf_propagator Santi, Santi_tmp, Seven, Sodd;
    suNg_propagator SQC, SQA, SQCA, SQAC, SQ2, Stmp_fund;
    

    complex corr[GLB_T][4][2]; // 4 for 2 gamma structure and their cross channels. 2 for even odd party projection
    
    _symplectic(Omega);
    _suNg_transpose(Otr,Omega);
    
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));
    
    // construct eas antisymmetrising matrix
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        memset(corr, 0, sizeof(corr));
        _propagator_zero(Santi);
        _propagator_zero_fund(SQC);
        _propagator_zero_fund(SQA);
        _propagator_zero_fund(SQCA);
        _propagator_zero_fund(SQAC);
        _propagator_zero_fund(SQ2);
        _propagator_zero_fund(Stmp_fund);
        
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                //lprintf("CHIMERA",0,"(a,beta,nm_anti,i_nas) = (%d,%d,%d,%d) \n",a,beta,nm_anti,i_nas);
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                //lprintf("CHIMERA",0,"(a,beta,nm_fund,i_nf) = (%d,%d,%d,%d) \n",a,beta,nm_fund,i_nf);
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ2, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ2);
            _propagator_transpose_fund(SQ2, Stmp_fund);
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ2);
            _propagator_transpose_fund(SQ2, Stmp_fund);
            
            // Gamma^1 = Cg5 = g0g2g5 ---> SQA
            
            _propagator_g5g0g2_fund(SQA, SQ2); // S g5g0g2
            _propagator_transpose_fund(Stmp_fund, SQA); // S^T g5^T g0^2 g2^T
            _propagator_g5g0g2_fund(SQA, Stmp_fund);
            
            // Gamma^1 = g0Cg5 = g2g5 ---> SQC
            
            _g5_propagator_fund(SQC, SQ2); // SQ1 = g5 S
            _g2_propagator_fund(Stmp_fund, SQC); // Stmp= g2g5 S
            _propagator_transpose_fund(SQC, Stmp_fund); // S^T g5^T g2^T
            _g5_propagator_fund(Stmp_fund, SQC); //  g5 S^T g5^T g2^T
            _g2_propagator_fund(SQC, Stmp_fund); //  g2 g5 S^T g5^T g2^T = (g2g5 S (g2g5)^T)^T
            
            // Gamma^1 = g0g2g5_g2g5 (cross channel) ---> SQCA
            
            _g5_propagator_fund(Stmp_fund, SQ2); // SQ1 = g5 S
            _g2_propagator_fund(SQCA,Stmp_fund); // Stmp= g2g5 S
            _propagator_transpose_fund(Stmp_fund, SQCA); // S^T g5^T g2^T
            _g5g0g2_propagator_fund(SQCA, Stmp_fund); // g5g0g2 S^T g5^T g2^T
            
            // Gamma^1 = g2g5_g0g2g5 (cross channel) ---> SQAC
            
            _g5g0g2_propagator_fund(Stmp_fund, SQ2); // g5g0g2 S
            _propagator_transpose_fund(SQAC, Stmp_fund); // S^T g5^T g0^T g2^T
            _g5_propagator_fund(Stmp_fund, SQAC); //  g5 S^T g5^T g0^T g2^T
            _g2_propagator_fund(SQAC, Stmp_fund); //  g2 g5 S^T g5^T g0^T g2^T
            
            // Gamma^2 = id & parity projection
            _g0_propagator(Santi_tmp, Santi);
            
            _propagator_add(Seven, Santi, Santi_tmp); // Stmp2 = (1+g0)S
            _propagator_mul_assign(Seven, 0.5);
            
            _propagator_sub(Sodd, Santi, Santi_tmp); // Stmp2 = (1-g0)S
            _propagator_mul_assign(Sodd, 0.5);
            
            
            for(int ip = 0; ip < NF; ++ip)
            for(int jp = 0; jp < NF; ++jp)
            {
                _complex_0(dqA);
                _complex_0(dqC);
                _complex_0(dqAC);
                _complex_0(dqCA);
                _complex_0(tmp);

                for(int l = 0; l < NG; ++l)
                for(int lp = 0; lp < NG; ++lp)
                for(int k = 0; k < NG; ++k)
                for(int kp = 0; kp < NG; ++kp)
                {
                    for(int alpha = 0; alpha < 4; ++alpha)
                    for(int gamma = 0; gamma < 4; gamma++)
                    {
                        tmp.re = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                        tmp.im = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                        
                        _complex_mul_assign(dqC, tmp, SQC.c[k].c[gamma].c[alpha].c[kp]);
                        _complex_mul_assign(dqA, tmp, SQA.c[k].c[gamma].c[alpha].c[kp]);
                        
                        _complex_mul_assign(dqAC, tmp, SQAC.c[k].c[gamma].c[alpha].c[kp]);
                        _complex_mul_assign(dqCA, tmp, SQCA.c[k].c[gamma].c[alpha].c[kp]);
                    }
                }
            
                for(int i = 0; i < 4; i++){
                    // A
                    _complex_mul_assign(corr[tc][0][0], Seven.c[ip].c[i].c[i].c[jp], dqA);
                    _complex_mul_assign(corr[tc][0][1], Sodd.c[ip].c[i].c[i].c[jp], dqA);
                    
                    // C
                    _complex_mul_assign(corr[tc][1][0], Seven.c[ip].c[i].c[i].c[jp], dqC);
                    _complex_mul_assign(corr[tc][1][1], Sodd.c[ip].c[i].c[i].c[jp], dqC);
                    
                    // AC
                    _complex_mul_assign(corr[tc][2][0], Seven.c[ip].c[i].c[i].c[jp], dqAC);
                    _complex_mul_assign(corr[tc][2][1], Sodd.c[ip].c[i].c[i].c[jp], dqAC);
                    
                    // CA
                    _complex_mul_assign(corr[tc][3][0], Seven.c[ip].c[i].c[i].c[jp], dqCA);
                    _complex_mul_assign(corr[tc][3][1], Sodd.c[ip].c[i].c[i].c[jp], dqCA);
                }
            }
        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP
    
    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n");
    
    global_sum((double*)corr, 2*2*4*GLB_T);
    
    // Print correlator
    
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    
    for (int i=0; i<4; i++){
        lprintf("CHIMERA",0,"conf #%d %s Chimera_%s_even_re= ", conf_num, label_plus,channel_name[i]);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][i][0].re);
        }
        
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_%s_even_im= ", conf_num, label_plus,channel_name[i]);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][i][0].im);
        }
        
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_%s_odd_re= ", conf_num, label_plus,channel_name[i]);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][i][1].re);
        }
        
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_%s_odd_im= ", conf_num, label_plus,channel_name[i]);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][i][1].im);
        }
        lprintf("CHIMERA",0,"\n");
        
    }
        
    //memset(corr, 0, sizeof(corr));
    }// END MASS LOOP
}

void contract_chimera_OP(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    int ix, tc, idx;
    char label_plus[256];
    complex dqB, tmp;
    suNg Omega, Otr;
    suNf_propagator Santi, Santi_tmp, Santi_tmp2, Seven, Sodd;
    suNg_propagator SQB, SQ2, Stmp_fund;
    complex corr[GLB_T][2]; //  2 for even odd party projection
    
    _symplectic(Omega);_suNg_transpose(Otr,Omega);
    
    // construct eas antisymmetrising matrix
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    // Mass LOOP
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        memset(corr, 0, sizeof(corr));
        _propagator_zero(Santi);_propagator_zero(Santi_tmp);_propagator_zero(Santi_tmp2);
        _propagator_zero(Seven);_propagator_zero(Sodd);
        _propagator_zero_fund(SQB);_propagator_zero_fund(SQ2);_propagator_zero_fund(Stmp_fund);
        
    
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ2, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ2);
            _propagator_transpose_fund(SQ2, Stmp_fund);
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ2);
            _propagator_transpose_fund(SQ2, Stmp_fund);
            
            // Gamma^1 = Cg5 = g0g2g5 ---> SQB
            
            _propagator_g0g2_fund(SQB, SQ2);
            _propagator_transpose_fund(Stmp_fund, SQB);
            _propagator_g0g2_fund(SQB, Stmp_fund); // SQ1 = ( C S (C)^T )^T
            
            
            // Gamma^2 = g5 & parity projection
            _g5_propagator(Seven, Santi);
            _g0_propagator(Santi_tmp, Seven);
            _propagator_add(Santi_tmp2, Seven, Santi_tmp); // Stmp2 = (1+g0) g5 S
            _g5_propagator(Seven, Santi_tmp2); // Seven = g5 (1+g0) g5 S
            _propagator_mul_assign(Seven, 0.5);
            
            
            _g5_propagator(Sodd, Santi);
            _g0_propagator(Santi_tmp, Sodd);
            _propagator_sub(Santi_tmp2, Sodd, Santi_tmp); // Stmp2 = (1-g0) g5 S
            _g5_propagator(Sodd, Santi_tmp2); // Sodd = g5 (1-g0) g5 S
            _propagator_mul_assign(Sodd, 0.5);
            
            for(int ip = 0; ip < NF; ++ip)
            for(int jp = 0; jp < NF; ++jp)
            {
                _complex_0(dqB);
                _complex_0(tmp);

                for(int l = 0; l < NG; ++l)
                for(int lp = 0; lp < NG; ++lp)
                for(int k = 0; k < NG; ++k)
                for(int kp = 0; kp < NG; ++kp)
                {
                    for(int alpha = 0; alpha < 4; ++alpha)
                    for(int gamma = 0; gamma < 4; gamma++)
                    {
                        tmp.re = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                        tmp.im = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                        
                        _complex_mul_assign(dqB, tmp, SQB.c[k].c[gamma].c[alpha].c[kp]);
                        
                    }
                }
            
                for(int i = 0; i < 4; i++){
                    _complex_mul_assign(corr[tc][0], Seven.c[ip].c[i].c[i].c[jp], dqB);
                    _complex_mul_assign(corr[tc][1], Sodd.c[ip].c[i].c[i].c[jp], dqB);

                }
            }
        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP

    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n");
    
    global_sum((double*)corr, 2*2*GLB_T);
        
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    
    // Print correlator
    
    lprintf("CHIMERA",0,"conf #%d %s Chimera_OB_even_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr[t][0].re);
    }
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OB_even_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr[t][0].im);
    }
    lprintf("CHIMERA",0,"\n");
    
    lprintf("CHIMERA",0,"conf #%d %s Chimera_OB_odd_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr[t][1].re);
    }
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OB_odd_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr[t][1].im);
    }
    lprintf("CHIMERA",0,"\n");
        
    } // END Mass LOOP
    
}

void contract_chimera_OV_noP(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    int ix, tc, idx;
    char label_plus[256];
    complex dq, tmp;
    suNg Omega, Otr;
    suNf_propagator Santi;
    suNg_propagator SQ1, SQ2, SQ, Stmp_fund;
    complex corr[GLB_T][3][3]; // Spin 1/2 state:  3 for g_i structure. 2 for even odd party projection
    
    _symplectic(Omega);
    _suNg_transpose(Otr,Omega);
    
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));

    // construct eas antisymmetrising matrix
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        memset(corr, 0, sizeof(corr));
        _propagator_zero(Santi);
        _propagator_zero_fund(SQ);_propagator_zero_fund(SQ1);
        _propagator_zero_fund(SQ2);_propagator_zero_fund(Stmp_fund);
        
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ); // O^T S
            _propagator_transpose_fund(SQ, Stmp_fund); // S^T O
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ); // O S^T O
            _propagator_transpose_fund(SQ, Stmp_fund); // O^T S O^T
        
            // C^T S^T C
            _propagator_g0g2_fund(Stmp_fund, SQ); //  S C
            _propagator_transpose_fund(SQ1, Stmp_fund); //  C^T S^T
            _propagator_g0g2_fund(Stmp_fund, SQ1); //  Stmp = C^T S^T C
            
            for (int mu=0;mu<3;mu++) {
                for (int nu=0;nu<3;nu++){
                    
                    // SQ2*SQ1 = g_mu S g_nu C^T S^T C
                    op_propagator_fund(&SQ1, &Stmp_fund, nu+1); // SQ1 = g_nu C^T S^T C
                    op_propagator_fund(&SQ2, &SQ, mu+1);        // SQ2 = g_mu S
                    
                    
                    // parity projection
                    /*
                    _g0_propagator(Santi_tmp, Santi);
                    _propagator_add(Seven12, Santi_tmp, Santi); // Stmp2 = (1+g0) Sas
                    _propagator_mul_assign(Seven12, 0.5);
                    */
                    
                    // Contraction
                    for(int ip = 0; ip < NF; ++ip)
                    for(int jp = 0; jp < NF; ++jp)
                    {
                        _complex_0(dq);
                        _complex_0(tmp);

                        for(int l = 0; l < NG; ++l)
                        for(int lp = 0; lp < NG; ++lp)
                        for(int k = 0; k < NG; ++k)
                        for(int kp = 0; kp < NG; ++kp)
                        {
                            for(int alpha = 0; alpha < 4; ++alpha)
                            for(int gamma = 0; gamma < 4; gamma++)
                            {
                                tmp.re = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                                tmp.im = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                                
                                _complex_mul_assign(dq, tmp, SQ1.c[k].c[gamma].c[alpha].c[kp]);
                                
                            }
                        }
                    
                        for(int i = 0; i < 4; i++){

                            // Trace spin
                            _complex_mul_assign(corr[tc][mu][nu], Santi.c[ip].c[i].c[i].c[jp], dq);
                        }
                    } //END contraction
                    
                } //END lnz LOOP
            } //END mu & nu LOOP
        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP

    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n NO projections: \n");
    
    global_sum((double*)corr, 2*3*3*GLB_T);
    
    // Print correlator
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    
    for (int mu=0;mu<3;mu++) for (int nu=0;nu<3;nu++) {
        
        lprintf("CHIMERA",0,"conf #%d %s Chimera_OV_g%d_g%d_re= ", conf_num, label_plus, mu+1, nu+1);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][mu][nu].re);
        } //
        
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV_g%d_g%d_im= ", conf_num, label_plus, mu+1, nu+1);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t][mu][nu].im);
        } //
        lprintf("CHIMERA",0,"\n");
    }
    } // END Mass LOOP

}

//#define PNORM 1.00//0.333333333333333333

void contract_chimera_OV(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    int ix, tc, idx;
    char label_plus[256];
    complex dq, tmp;

    suNg Omega, Otr;
    suNf_propagator Santi, Santi_tmp, Santi_tmp2, Seven12, Sodd12, Seven32, Sodd32;
    suNg_propagator SQ1, SQ2, SQ, Stmp_fund;
    complex corr_12[GLB_T][3][3][2]; // Spin 1/2 state:  3 for g_i structure. 2 for even odd party projection
    complex corr_32[GLB_T][3][3][2]; // Spin 3/2 state:  3 for g_i structure. 2 for even odd party projection
    
    _symplectic(Omega);_suNg_transpose(Otr,Omega);
    
    // construct eas antisymmetrising matrix
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        
        memset(corr_12, 0, sizeof(corr_12));memset(corr_32, 0, sizeof(corr_32));
        _propagator_zero(Santi);_propagator_zero(Santi_tmp);_propagator_zero(Santi_tmp2);
        _propagator_zero(Seven12);_propagator_zero(Sodd12);_propagator_zero(Seven32);_propagator_zero(Sodd32);
        _propagator_zero_fund(SQ);_propagator_zero_fund(SQ1);
        _propagator_zero_fund(SQ2);_propagator_zero_fund(Stmp_fund);
        
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ); // O^T S
            _propagator_transpose_fund(SQ, Stmp_fund); // S^T O
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ); // O S^T O
            _propagator_transpose_fund(SQ, Stmp_fund); // O^T S O^T
        
            // C^T S^T C
            _propagator_g0g2_fund(Stmp_fund, SQ); //  S C
            _propagator_transpose_fund(SQ1, Stmp_fund); //  C^T S^T
            _propagator_g0g2_fund(Stmp_fund, SQ1); //  Stmp = C^T S^T C
            
            for (int mu=0;mu<3;mu++) for (int nu=0;nu<3;nu++){
                for (int lnz=0;lnz<3;lnz++){
                    
                    // SQ2*SQ1 = g_lnz S g_nu C^T S^T C
                    op_propagator_fund(&SQ1, &Stmp_fund, nu+1); // SQ1 = g_nu C^T S^T C
                    op_propagator_fund(&SQ2, &SQ, lnz+1);        // SQ2 = g_lnz S
                    
                    // Spin 1/2 projection
                    op_propagator(&Santi_tmp, &Santi, lnz+1); // Sas = g_lnz S
                    op_propagator(&Santi_tmp2, &Santi_tmp, mu+1); // Sas = g_mu g_lnz S
                    _propagator_mul_assign(Santi_tmp2, 0.3333333333333333);
                    
                    // -->  g_mu g_lnz S_as [g_lnz C^T S^T C  g_nu S]
                    
                    // parity projection
                    
                    _g0_propagator(Santi_tmp, Santi_tmp2);
                    
                    _propagator_add(Seven12, Santi_tmp2, Santi_tmp); // Stmp2 = (1+g0) Sas
                    _propagator_mul_assign(Seven12, 0.5);
                    
                    _propagator_sub(Sodd12, Santi_tmp2, Santi_tmp); // Stmp2 = (1-g0) Sas
                    _propagator_mul_assign(Sodd12, 0.5);
                    
                    // Spin 3/2 projection
                    if (mu == lnz){
                        _propagator_sub(Santi_tmp, Santi, Santi_tmp2); // Sas = S - 1/3 g_mu g_nu S
                    }
                    else{
                        _id_propagator(Santi_tmp, Santi_tmp2); // not so smart here yet
                        _propagator_mul_assign(Santi_tmp, -1); // Sas =  - 1/3 g_mu g_nu S
                    }
                    
                    // parity projection
                    
                    _g0_propagator(Santi_tmp2, Santi_tmp);
                    
                    _propagator_add(Seven32, Santi_tmp, Santi_tmp2); // Stmp2 = (1+g0) Sas
                    _propagator_mul_assign(Seven32, 0.5);
                    
                    _propagator_sub(Sodd32, Santi_tmp, Santi_tmp2); // Stmp2 = (1-g0) Sas
                    _propagator_mul_assign(Sodd32, 0.5);
                    
                    
                    // Contraction
                    for(int ip = 0; ip < NF; ++ip)
                    for(int jp = 0; jp < NF; ++jp)
                    {
                        _complex_0(dq);
                        _complex_0(tmp);

                        for(int l = 0; l < NG; ++l)
                        for(int lp = 0; lp < NG; ++lp)
                        for(int k = 0; k < NG; ++k)
                        for(int kp = 0; kp < NG; ++kp)
                        {
                            for(int alpha = 0; alpha < 4; ++alpha)
                            for(int gamma = 0; gamma < 4; gamma++)
                            {
                                tmp.re = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                                tmp.im = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                                
                                _complex_mul_assign(dq, tmp, SQ1.c[k].c[gamma].c[alpha].c[kp]);
                                
                            }
                        }
                    
                        for(int i = 0; i < 4; i++){

                            // Trace spin
                            _complex_mul_assign(corr_12[tc][mu][nu][0], Seven12.c[ip].c[i].c[i].c[jp], dq);
                            _complex_mul_assign(corr_12[tc][mu][nu][1], Sodd12.c[ip].c[i].c[i].c[jp], dq);
                            _complex_mul_assign(corr_32[tc][mu][nu][0], Seven32.c[ip].c[i].c[i].c[jp], dq);
                            _complex_mul_assign(corr_32[tc][mu][nu][1], Sodd32.c[ip].c[i].c[i].c[jp], dq);
                            
                        }
                        
                    } //END contraction
                } //END lnz LOOP
            } //END mu & nu LOOP
        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP

    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n With Spin and parity projection:\n");
    
    global_sum((double*)corr_12, 2*GLB_T*3*3*2);
    global_sum((double*)corr_32, 2*GLB_T*3*3*2);
    
    // Print correlator
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    
    lprintf("CHIMERA",0,"conf #%d %s Chimera_OV12_even_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0][0].re+corr_12[t][1][1][0].re+corr_12[t][2][2][0].re);
    } //
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV12_even_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0][0].im+corr_12[t][1][1][0].im+corr_12[t][2][2][0].im);
    }//
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV12_odd_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0][1].re+corr_12[t][1][1][1].re+corr_12[t][2][2][1].re);
    }
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV12_odd_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0][1].im+corr_12[t][1][1][1].im+corr_12[t][2][2][1].im);
    }
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_even_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0][0].re+corr_32[t][1][1][0].re+corr_32[t][2][2][0].re);
    }
        
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_even_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0][0].im+corr_32[t][1][1][0].im+corr_32[t][2][2][0].im);
    }
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_odd_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0][1].re+corr_32[t][1][1][1].re+corr_32[t][2][2][1].re);
    }
        
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_odd_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0][1].im+corr_32[t][1][1][1].im+corr_32[t][2][2][1].im);
    }
    
    lprintf("CHIMERA",0,"\n");
        
    } // END Mass LOOP
    
    
}

void contract_chimera_OV_Ps(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    int ix, tc, idx;
    char label_plus[256];
    complex dq, tmp;
    suNg Omega, Otr;
    suNf_propagator Santi, Santi_tmp, Santi_tmp2;
    suNg_propagator SQ1, SQ2, SQ, Stmp_fund;

    complex corr_12[GLB_T][3][3]; // Spin 1/2 state:  3 for g_i structure.
    complex corr_32[GLB_T][3][3]; // Spin 3/2 state:  3 for g_i structure.
    
    _symplectic(Omega);
    _suNg_transpose(Otr,Omega);
    
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));
    
    
    // construct eas antisymmetrising matrix
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        memset(corr_12, 0, sizeof(corr_12));
        memset(corr_32, 0, sizeof(corr_32));
        _propagator_zero(Santi);_propagator_zero(Santi_tmp);_propagator_zero(Santi_tmp2);
        _propagator_zero_fund(SQ);_propagator_zero_fund(SQ1);
        _propagator_zero_fund(SQ2);_propagator_zero_fund(Stmp_fund);
        
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ); // O^T S
            _propagator_transpose_fund(SQ, Stmp_fund); // S^T O
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ); // O S^T O
            _propagator_transpose_fund(SQ, Stmp_fund); // O^T S O^T
        
            // C^T S^T C
            _propagator_g0g2_fund(Stmp_fund, SQ); //  S C
            _propagator_transpose_fund(SQ1, Stmp_fund); //  C^T S^T
            _propagator_g0g2_fund(Stmp_fund, SQ1); //  Stmp = C^T S^T C
            
            for (int mu=0;mu<3;mu++) for (int nu=0;nu<3;nu++){
                for (int lnz=0;lnz<3;lnz++){
                    
                    // SQ2*SQ1 = g_lnz S g_nu C^T S^T C
                    op_propagator_fund(&SQ1, &Stmp_fund, nu+1); // SQ1 = g_nu C^T S^T C
                    op_propagator_fund(&SQ2, &SQ, lnz+1);        // SQ2 = g_lnz S
                    
                    // Spin 1/2 projection
                    op_propagator(&Santi_tmp, &Santi, lnz+1); // Sas = g_lnz S
                    op_propagator(&Santi_tmp2, &Santi_tmp, mu+1); // Sas = g_mu g_lnz S
                    _propagator_mul_assign(Santi_tmp2, 0.3333333333333333);
                    
                    // -->  g_mu g_lnz S_as [g_lnz C^T S^T C  g_nu S]
                    
                    
                    // Spin 3/2 projection
                    if (mu == lnz){
                        _propagator_sub(Santi_tmp, Santi, Santi_tmp2); // Sas = S - 1/3 g_mu g_nu S
                    }
                    else{
                        _id_propagator(Santi_tmp, Santi_tmp2); // not so smart here yet
                        _propagator_mul_assign(Santi_tmp, -1); // Sas =  - 1/3 g_mu g_nu S
                    }
                    
                    
                    // Contraction
                    for(int ip = 0; ip < NF; ++ip)
                    for(int jp = 0; jp < NF; ++jp)
                    {
                        _complex_0(dq);
                        _complex_0(tmp);

                        for(int l = 0; l < NG; ++l)
                        for(int lp = 0; lp < NG; ++lp)
                        for(int k = 0; k < NG; ++k)
                        for(int kp = 0; kp < NG; ++kp)
                        {
                            for(int alpha = 0; alpha < 4; ++alpha)
                            for(int gamma = 0; gamma < 4; gamma++)
                            {
                                tmp.re = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                                tmp.im = eas[ip][l][kp]*SQ2.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                                
                                _complex_mul_assign(dq, tmp, SQ1.c[k].c[gamma].c[alpha].c[kp]);
                                
                            }
                        }
                    
                        for(int i = 0; i < 4; i++){

                            // Trace spin
                            _complex_mul_assign(corr_12[tc][mu][nu], Santi_tmp2.c[ip].c[i].c[i].c[jp], dq);
                            _complex_mul_assign(corr_32[tc][mu][nu], Santi_tmp.c[ip].c[i].c[i].c[jp], dq);
                            
                        }
                        
                    } //END contraction
                } //END lnz LOOP
            } //END mu & nu LOOP
        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP

    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n With Spin projection only:\n");
    
    global_sum((double*)corr_12, 2*GLB_T*3*3);
    global_sum((double*)corr_32, 2*GLB_T*3*3);
    
    // Print correlator
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    
    lprintf("CHIMERA",0,"conf #%d %s Chimera_OV12_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0].re+corr_12[t][1][1].re+corr_12[t][2][2].re);
    } //
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV12_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_12[t][0][0].im+corr_12[t][1][1].im+corr_12[t][2][2].im);
    }//
    
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_re= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0].re+corr_32[t][1][1].re+corr_32[t][2][2].re);
    }
        
    lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV32_im= ", conf_num, label_plus);
    for(int t = 0; t < GLB_T; t++){
        lprintf("CHIMERA",0,"%e ",corr_32[t][0][0].im+corr_32[t][1][1].im+corr_32[t][2][2].im);
    }
    lprintf("CHIMERA",0,"\n");lprintf("CHIMERA",0,"\n");
    
    for (int mu=0;mu<3;mu++) for (int nu=0;nu<3;nu++) {
        
        lprintf("CHIMERA",0,"conf #%d %s Chimera_OV12_mu%d_nu%d_re= ", conf_num, label_plus, mu+1, nu+1);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr_12[t][mu][nu].re);
        } //
        
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV12_mu%d_nu%d_im= ", conf_num, label_plus, mu+1, nu+1);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr_12[t][mu][nu].im);
        } //
        lprintf("CHIMERA",0,"\n");
    }

    } // END Mass LOOP
    
}

void contract_chimera_OV_JW(spinor_field_fund* psi_fund, spinor_field* psi_anti, int nm_fund, int nm_anti, double* m_fund, double* m_anti, int tau, int conf_num, char* label){
    
    int ix, tc, idx;
    char label_plus[256];
    complex dq, tmp;
    suNg Omega, Otr;
    suNf_propagator Santi;
    suNg_propagator SQ1, SQ, Stmp_fund;
    complex corr[GLB_T]; // Spin 1/2 state:  3 for g_i structure. 2 for even odd party projection
    
    _symplectic(Omega);
    _suNg_transpose(Otr,Omega);
    
    double eas[NF][NG][NG];
    memset(eas, 0, sizeof(eas));

    // construct eas antisymmetrising matrix
    int A, B, C;
    C = 0;
    for(A = 1; A < NG; A++)
    for(B = 0; B < A; B++){
        if(A == B+NG/2 ){// diagonal of the submatrices
            if(B==0) continue;
            double entry = 1.0/sqrt(2*B*(B+1));
            for(int d = 0;d<B;++d) // diagonal
            {
                eas[C][d][d+NG/2] =  entry;
                eas[C][d+NG/2][d] = -entry;
            }
            
            eas[C][B][B+NG/2] = -B*entry;
            eas[C][B+NG/2][B] = B*entry;
            }
            
        else{
            eas[C][A][B] = 0.7071067811865476;
            eas[C][B][A] = -0.7071067811865476;
        }
        C ++;
    }
    
    for (int i_nf=0; i_nf < nm_fund; i_nf++)for (int i_nas=0; i_nas < nm_anti; i_nas++){
        
        memset(corr, 0, sizeof(corr));
        _propagator_zero(Santi);
        _propagator_zero_fund(SQ);_propagator_zero_fund(SQ1);_propagator_zero_fund(Stmp_fund);
        
    for(int t = 0; t < T; t++)
    {
        tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
        
        
        for(int x = 0; x < X; x++)
        for(int y = 0; y < Y; y++)
        for(int z = 0; z < Z; z++)
        {
            ix = ipt(t, x, y, z);
            
            for(int a = 0; a < NF; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = 4*a*nm_anti+beta*nm_anti + i_nas;
                _propagator_assign(Santi, *_FIELD_AT(&psi_anti[idx],ix), a, beta);
            }
            
            for(int a = 0; a < NG; ++a)
            for(int beta = 0; beta < 4; beta++)
            {
                idx = a*4*nm_fund +beta*nm_fund + i_nf;
                _propagator_assign_fund(SQ, *_FIELD_AT(&psi_fund[idx],ix), a, beta);
            }
            
            _suNf_prop_multiply_fund(Stmp_fund, Otr, SQ); // O^T S
            _propagator_transpose_fund(SQ, Stmp_fund); // S^T O
            _suNf_prop_multiply_fund(Stmp_fund, Omega, SQ); // O S^T O
            _propagator_transpose_fund(SQ, Stmp_fund); // O^T S O^T
        
            _g5g3_propagator_fund(Stmp_fund, SQ);
            _propagator_g0_fund(SQ1, Stmp_fund);
                    
                    // Contraction
                    for(int ip = 0; ip < NF; ++ip)
                    for(int jp = 0; jp < NF; ++jp)
                    {
                        _complex_0(dq);
                        _complex_0(tmp);

                        for(int l = 0; l < NG; ++l)
                        for(int lp = 0; lp < NG; ++lp)
                        for(int k = 0; k < NG; ++k)
                        for(int kp = 0; kp < NG; ++kp)
                        {
                            for(int alpha = 0; alpha < 4; ++alpha)
                            for(int gamma = 0; gamma < 4; gamma++)
                            {
                                tmp.re = eas[ip][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k];
                                tmp.im = eas[ip][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k];
                                
                                _complex_mul_assign(dq, tmp, SQ1.c[kp].c[alpha].c[gamma].c[k]);
                                
                            }
                        }
                    
                        for(int i = 0; i < 4; i++){

                            // Trace spin
                            _complex_mul_assign(corr[tc], Santi.c[ip].c[i].c[i].c[jp], dq);
                        }
                    } //END contraction

        }// END SPATIAL LOOP
    } // END TEMPORAL LOOP

    lprintf("CHIMERA",0,"CHIMERA BARYON CONTRACTION DONE!\n NO projections: \n");
    
    global_sum((double*)corr, 2*GLB_T);
    
    // Print correlator
    sprintf(label_plus, "mf=%2.6f mas=%2.6f %s", m_fund[i_nf], m_anti[i_nas], label);
    

        lprintf("CHIMERA",0,"conf #%d %s Chimera_OV_JW_re= ", conf_num, label_plus);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t].re);
        } //
        
    
        lprintf("CHIMERA",0,"\nconf #%d %s Chimera_OV_JW_im= ", conf_num, label_plus);
        for(int t = 0; t < GLB_T; t++){
            lprintf("CHIMERA",0,"%e ",corr[t].im);
        } //
        
        lprintf("CHIMERA",0,"\n");
        
    } // END Mass LOOP

}
