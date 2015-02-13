#include "cuda_Eulerian.h"
#include <math.h>

__global__ void kernel_numberdensity_particle_velz(real val, real* wp, real *w, dom_struct *dom)
{
	int ti = blockIdx.x * blockDim.x + threadIdx.x;
	int tj = blockIdx.y * blockDim.y + threadIdx.y;
	
	for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
		if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
			wp[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] + val;
		}
	}
}

__global__ void kernel_numberdensity_march(real dt, dom_struct *dom, real *numden, real *numden1, real *ux, real *uy, real *uz)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
	int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;
	
	if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
		for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
			
			int C   = i       + tj      *dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			int Cx0 = (i - 1) + tj      *dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			int Cx1 = (i + 1) + tj      *dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			int Cy0 = i       + (tj - 1)*dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			int Cy1 = i       + (tj + 1)*dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			int Cz0 = i       + tj      *dom->Gcc._s1b + (tk - 1)*dom->Gcc._s2b;
			int Cz1 = i       + tj      *dom->Gcc._s1b + (tk + 1)*dom->Gcc._s2b;
			int fx0 = i       + tj      *dom->Gfx._s1b + tk      *dom->Gfx._s2b;
			int fx1 = (i + 1) + tj      *dom->Gfx._s1b + tk      *dom->Gfx._s2b;
			int fy0 = i       + tj      *dom->Gfy._s1b + tk      *dom->Gfy._s2b;
			int fy1 = i       + (tj + 1)*dom->Gfy._s1b + tk      *dom->Gfy._s2b;
			int fz0 = i       + tj      *dom->Gfz._s1b + tk      *dom->Gfz._s2b;
			int fz1 = i       + tj      *dom->Gfz._s1b + (tk + 1)*dom->Gfz._s2b;
			
			real n_x1 = 0.5 * (1 - copysign(1.0, ux[fx1])) * numden[Cx1] + 0.5 * (1 + copysign(1.0, ux[fx1])) * numden[C];
			real n_y1 = 0.5 * (1 - copysign(1.0, uy[fy1])) * numden[Cy1] + 0.5 * (1 + copysign(1.0, uy[fy1])) * numden[C];
			real n_z1 = 0.5 * (1 - copysign(1.0, uz[fz1])) * numden[Cz1] + 0.5 * (1 + copysign(1.0, uz[fz1])) * numden[C];
			
			real n_x0 = 0.5 * (1 - copysign(1.0, ux[fx0])) * numden[C] + 0.5 * (1 + copysign(1.0, ux[fx0])) * numden[Cx0];
			real n_y0 = 0.5 * (1 - copysign(1.0, uy[fy0])) * numden[C] + 0.5 * (1 + copysign(1.0, uy[fy0])) * numden[Cy0];
			real n_z0 = 0.5 * (1 - copysign(1.0, uz[fz0])) * numden[C] + 0.5 * (1 + copysign(1.0, uz[fz0])) * numden[Cz0];
			
			real nx = dt / dom->dx * (ux[fx1] * n_x1 - ux[fx0] * n_x0);
			real ny = dt / dom->dy * (uy[fy1] * n_y1 - uy[fy0] * n_y0);
			real nz = dt / dom->dz * (uz[fz1] * n_z1 - uz[fz0] * n_z0);
			numden1[C] = numden[C] - nx - ny - nz;
		}
	}
}

__global__ void kernel_numberdensity_update_numden(dom_struct *dom, real *numden, real *numden1)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
	int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;
	
	if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
		for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
			int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
			numden[C] = numden1[C];
		}
	}
}
