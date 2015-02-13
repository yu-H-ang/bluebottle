#ifndef _EULERIAN_H
#define _EULERIAN_H

#include "bluebottle.h"

//======================================================================
// number density equation

#define UNIFORM 0

extern real bubble_density;
extern real bubble_radius;
extern int bubble_init_cond;
extern real bubble_init_cond_uniform_m;
typedef struct numdenBC_struct {
  int nW;
  int nE;
  int nS;
  int nN;
  int nB;
  int nT;
} numdenBC_struct;
extern numdenBC_struct numdenBC;
extern real *numden;
extern real *u_p;
extern real *v_p;
extern real *w_p;
extern real **_numden;
extern real **_nextnumden;
extern real **_u_p;
extern real **_v_p;
extern real **_w_p;

void numberdensity_read_input(void);
int numberdensity_init(void);
void numberdensity_clean(void);
void cuda_numberdensity_malloc(void);
void cuda_numberdensity_push(void);
void cuda_numberdensity_pull(void);
void cuda_numberdensity_free(void);
void cuda_numberdensity_BC(void);
void cuda_numberdensity_march(void);
void cuda_numberdensity_particle_velz(void);

//======================================================================
// concentration equation
extern real concen_diff;
extern real concen_diss;
extern int concen_init_cond;
extern real concen_init_cond_uniform_m;
extern numdenBC_struct concenBC;
extern real *concen;
extern real **_concen;
extern real **_nextconcen;

void concentration_read_input(void);
int concentration_init(void);
void concentration_clean(void);
//======================================================================
#endif
