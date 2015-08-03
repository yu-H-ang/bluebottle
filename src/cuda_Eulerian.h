#ifndef _CUDA_EULERIAN_H
#define _CUDA_EULERIAN_H

extern "C"
{
#include "bluebottle.h"
#include "Eulerian.h"
}

__global__ void kernel_numberdensity_particle_velz(dom_struct *dom,
                                                   real *velpz,
                                                   real *velz,
                                                   real *ter);
 
__global__ void kernel_march_numberdensity(real dt,
                                           dom_struct *dom,
                                           real *numden,
                                           real *numden1,
                                           real *ux,
                                           real *uy,
                                           real *uz,
                                           real *numdot);
 
__global__ void kernel_numberdensity_update_numden(dom_struct *dom,
                                                   real *numden,
                                                   real *numden1);

__global__ void kernel_fz_coupling_numden_generate(real *z_nd,
                                                   real *nd,
                                                   real *wb,
                                                   dom_struct *dom);
 
__global__ void kernel_march_concentration(real dt,
                                           dom_struct *dom,
                                           real *conc,
                                           real *conc1,
                                           real *ux,
                                           real *uy,
                                           real *uz,
                                           real *mdot,
                                           real *numd,
                                           real ccdiff);
  
__global__ void kernel_concentration_update_concen(dom_struct *dom,
                                                   real *conc,
                                                   real *conc1);
 
__global__ void kernel_numden_inner_copy(dom_struct *dom,
                                         real *numd,
                                         real *numd_tmp);

__global__ void kernel_compute_mdot(dom_struct *dom,
                                    real *numd,
                                    real *mas,
                                    real *conc,
                                    real *dia,
                                    real *ter,
                                    real *mdot,
                                    real *ccsat,
                                    real ccdiff,
                                    real nu);

__global__ void BC_p_B_D(real *p,
                         dom_struct *dom,
                         real bc);

__global__ void BC_p_T_D(real *p,
                         dom_struct *dom,
                         real bc);

__global__ void kernel_compute_bubble_diameter(dom_struct *dom,
                                               real *mas,
                                               real *numd,
                                               real *dia,
                                               real *rho_b,
                                               real rho_fluid,
                                               real pre_a,
                                               real rho_a,
                                               real gravacc);

__global__ void kernel_forcing_add_z_field_bubble(real scale,
                                                  real *ndfz,
                                                  real *diafz,
                                                  real *fz,
                                                  dom_struct *dom,
                                                  real *bdenf,
                                                  real rho_fluid);


__global__ void kernel_inner_scalarfield_update_x(dom_struct *dom,
                                                  real *left,
                                                  real *right);
 
__global__ void kernel_march_bubblemass(real dt,
                                        dom_struct *dom,
                                        real *bubm,
                                        real *bubm1,
                                        real *ux,
                                        real *uy,
                                        real *uz,
                                        real *numd,
                                        real *mdot,
                                        real *bubgenmdot);
 
 __global__ void kernel_compute_bubble_diameterfz(dom_struct *dom,
                                                 real *diafz,
                                                 real *dia,
                                                 real *wb);

__global__ void kernel_compute_terminal_velocity(dom_struct *dom,
                                                 real *ter,
                                                 real *bdia,
                                                 real *bubd,
                                                 real rho_fluid,
                                                 real constant,
                                                 real limit);

__global__ void BC_p_B_D_square(real *p,
                                dom_struct *dom,
                                real bc);

__global__ void BC_p_B_D_Gaussian(real *p,
                                dom_struct *dom,
                                real *bc);

__global__ void kernel_compute_bubble_generator(dom_struct *dom,
                                                real *left,
                                                real *right,
                                                real sc);

#endif
