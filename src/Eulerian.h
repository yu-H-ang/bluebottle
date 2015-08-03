#ifndef _EULERIAN_H
#define _EULERIAN_H

#include "bluebottle.h"
//======================================================================
// Eulerian

void Eulerian_read_input(void);
void Eulerian_show_config(void);
int Eulerian_init(void);
int Eulerian_init_parameters(void);
void Eulerian_clean(void);
void cuda_Eulerian_malloc(void);
void cuda_Eulerian_free(void);
void cuda_Eulerian_push(void);
void cuda_Eulerian_pull(void);

void cuda_add_terminal_velz(void);
void cuda_compute_mdot(void);
void cuda_coupling_force_preparation(void);
void cuda_compute_coupling_force(void);

void cgns_flow_field_Eulerian(void);
void cgns_flow_field_Eulerian2(void);
void cgns_grid_Eulerian(void);

void in_restart_Eulerian(void);
void out_restart_Eulerian(void);

extern int turb_switch;
#define ON 1
#define OFF -1
int domain_init_Eulerian(void);
//======================================================================
// number density equation

#define UNIFORM 0
#define RANDOM 271828
#define GAUSSIAN 0
#define HYPERBOLICTAN 1

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
	real nBDt;
	int nT;
	real nTD;
} numdenBC_struct;
extern numdenBC_struct numdenBC;
extern real *numden;
extern real totalnumden;
extern real *w_b;
extern real *ter_cell;
extern real **_numden;
extern real **_nextnumden;
extern real **_w_b;
extern real **_ter_cell;
extern real **_f_z_coupling_numden;
extern real *numGau;
extern real **_numGau;
extern real **_numGau_temp;
extern real *masGau;
extern real **_masGau;
extern real **_masGau_temp;
typedef struct BubGen_struct {
	int  BGis;
	int  BGie;
	int  BGjs;
	int  BGje;
	int  BGks;
	int  BGke;
	int  BubGen_type;
	real BubGen_dia;
	real BubGen_bubblenumber;
	real BubGen_amplitude;
	real BubGen_sigmaX;
	real BubGen_sigmaY;
	real BubGen_sigmaZ;
	real BubGen_x0;
	real BubGen_y0;
	real BubGen_z0;
	real BubGen_Lx1;
	real BubGen_Lx2;
	real BubGen_Ly1;
	real BubGen_Ly2;
	real BubGen_epsx;
	real BubGen_epsy;
} BubGen_struct;
extern BubGen_struct BubbleGenerator;
extern real *BGndot;
extern real **_BGndot;
extern real *BGmdot;
extern real **_BGmdot;
extern real TerminalVelocityLimit;

void cuda_numberdensity_BC(void);
void cuda_num_mas_BC_compute(void);
void cuda_numberdensity_march(void);
void cuda_numberdensity_compute_totalnumden(void);
//======================================================================
// concentration equation

extern real concen_diff;
extern real concen_atm;
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
extern real *concen_sat;
extern real **_concen_sat;

void cuda_concentration_BC(void);
void cuda_concentration_march(void);
//======================================================================
// bubble mass equation

#define GRAVITY_X 1
#define GRAVITY_Y 2
#define GRAVITY_Z 3

extern real *bubmas;
extern real **_bubmas;
extern real **_nextbubmas;
extern real *bubdia;
extern real **_bubdia;
extern real **_bubdiafz;
extern int bubmas_init_cond;
extern numdenBC_struct bubmasBC;
extern real bubmas_init_cond_uniform_m;
extern int gravity_direction;
extern real *bubden;
extern real *bubden_face;
extern real **_bubden;
extern real **_bubden_face;

extern real pressure_atm;
extern real rho_atm;
extern real grav_acc;

void cuda_bubblemass_BC(void);
void cuda_bubblemass_march(void);
//======================================================================
#endif
