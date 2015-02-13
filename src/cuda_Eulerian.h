#ifndef _CUDA_EULERIAN_H
#define _CUDA_EULERIAN_H

extern "C"
{
#include "bluebottle.h"
#include "Eulerian.h"
}

__global__ void kernel_numberdensity_particle_velz(real val, real* wp, real *w, dom_struct *dom);
__global__ void kernel_numberdensity_march(real dt, dom_struct *dom, real *numden, real *numden1, real *ux, real *uy, real *uz);
__global__ void kernel_numberdensity_update_numden(dom_struct *dom, real *numden, real *numden1);


#endif
