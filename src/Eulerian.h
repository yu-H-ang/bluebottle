#ifndef _EULERIAN_H
#define _EULERIAN_H

#include "bluebottle.h"
//======================================================================
// Eulerian

void Eulerian_read_input(void);
void Eulerian_show_config(void);
int Eulerian_init(void);
void Eulerian_clean(void);
void cuda_Eulerian_malloc(void);
void cuda_Eulerian_free(void);
void cuda_Eulerian_push(void);
void cuda_Eulerian_pull(void);

void cuda_compute_bubble_diameter(void);
void cuda_compute_particle_velz(void);
void cuda_compute_coupling_forcing(void);
void cuda_compute_mdot(void);
//======================================================================
// number density equation

#define UNIFORM 0
#define RANDOM 271828

extern real bubble_density;
extern real bubble_radius;
extern int bubble_init_cond;
extern real bubble_init_cond_uniform_m;
extern int bubble_init_cond_random_min;
extern int bubble_init_cond_random_max;
typedef struct numdenBC_struct {
  int nW;
  int nE;
  int nS;
  int nN;
  int nB;
  real nBD;
  int nT;
  real nTD;
} numdenBC_struct;
extern numdenBC_struct numdenBC;
extern real *numden;
extern real totalnumden;
extern real *w_b;
extern real **_numden;
extern real **_nextnumden;
extern real **_w_b;
extern real **_f_z_coupling_numden;

void cuda_numberdensity_BC(void);
void cuda_numberdensity_march(void);
void cuda_numberdensity_compute_totalnumden(void);
//======================================================================
// concentration equation

extern real concen_diff;
extern real concen_diss;
extern int concen_init_cond;
extern real concen_init_cond_uniform_m;
extern int concen_init_cond_random_max;
extern int concen_init_cond_random_min;
extern numdenBC_struct concenBC;
extern real *concen;
extern real **_concen;
extern real **_nextconcen;
extern real **_velmag;
extern real **_mdot;

void cuda_concentration_BC(void);
void cuda_concentration_march(void);
//======================================================================
// bubble mass equation

extern real *bubmas;
extern real **_bubmas;
extern real **_nextbubmas;
extern real **_bubdia;
extern real **_bubdiafz;
extern int bubmas_init_cond;
extern numdenBC_struct bubmasBC;
extern real bubmas_init_cond_uniform_m;

extern real pressure_atm;
extern real rho_atm;
extern real grav_acc;

void cuda_bubblemass_BC(void);
void cuda_bubblemass_march(void);
//======================================================================
#endif
