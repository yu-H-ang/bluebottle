#ifndef _CUDA_EULERIAN_H
#define _CUDA_EULERIAN_H

extern "C"
{
#include "bluebottle.h"
#include "Eulerian.h"
}

__global__ void kernel_numberdensity_particle_velz(real val,
                                                   real* velpz,
                                                   real *velz,
                                                   real *dia,
                                                   dom_struct *dom);
 
__global__ void kernel_march_numberdensity(real dt,
                                           dom_struct *dom,
                                           real *numden,
                                           real *numden1,
                                           real *ux,
                                           real *uy,
                                           real *uz);
 
__global__ void kernel_numberdensity_update_numden(dom_struct *dom,
                                                   real *numden,
                                                   real *numden1);

__global__ void kernel_fz_coupling_numden_generate(real *z_nd,
                                                   real *nd,
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
                                    real *conc,
                                    real *dia,
                                    real *mdot,
                                    real *velmag,
                                    real *ux,
                                    real *uy, 
                                    real *uz,
                                    real ccdiss,
                                    real ccdiff,
                                    real nu);

__global__ void BC_p_B_D(real *p,
                         dom_struct *dom,
                         real bc);

__global__ void BC_p_T_D(real *p,
                         dom_struct *dom,
                         real bc);

__global__ void kernel_compute_bubble_diameter(dom_struct *dom,
                                               real *vol,
                                               real *dia);

__global__ void kernel_forcing_add_z_field_bubble(real scale,
                                                  real *nd,
                                                  real *dia,
                                                  real *fz,
                                                  dom_struct *dom);


__global__ void kernel_inner_scalarfield_update_x(dom_struct *dom,
                                                  real *left,
                                                  real *right);
 
__global__ void kernel_march_bubblevolume(real dt,
                                          dom_struct *dom,
                                          real *bubv,
                                          real *bubv1,
                                          real *ux,
                                          real *uy,
                                          real *uz,
                                          real *mdot,
                                          real rhog);
 
 __global__ void kernel_compute_bubble_diameterfz(dom_struct *dom,
                                                  real *dia,
                                                  real *diafz);
 
#endif
