#include "cuda_Eulerian.h"
#include <math.h>
#include "bluebottle.h"

__global__ void kernel_numberdensity_particle_velz(real val, real* velzp, real *velz, real *dia, dom_struct *dom)
{
	int ti = blockIdx.x * blockDim.x + threadIdx.x;
	int tj = blockIdx.y * blockDim.y + threadIdx.y;
	
	for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
		if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
			int C = ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b;
			velzp[C] = velz[C] + val * dia[C] * dia[C];
			
			//if(ti==32 && tj==32 && k<3) printf("compute vel: k==%d, diafz==%f, velz==%f\n", k, dia[C], velzp[C]);
		}
	}
}

__global__ void kernel_march_numberdensity(real dt,
                                           dom_struct *dom,
                                           real *numden,
                                           real *numden1,
                                           real *ux,
                                           real *uy,
                                           real *uz)
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
			
			real n_x1 = 0.5 * (1.0 - copysign(1.0, ux[fx1])) * numden[Cx1] + 0.5 * (1.0 + copysign(1.0, ux[fx1])) * numden[C];
			real n_y1 = 0.5 * (1.0 - copysign(1.0, uy[fy1])) * numden[Cy1] + 0.5 * (1.0 + copysign(1.0, uy[fy1])) * numden[C];
			real n_z1 = 0.5 * (1.0 - copysign(1.0, uz[fz1])) * numden[Cz1] + 0.5 * (1.0 + copysign(1.0, uz[fz1])) * numden[C];
			
			real n_x0 = 0.5 * (1.0 - copysign(1.0, ux[fx0])) * numden[C] + 0.5 * (1.0 + copysign(1.0, ux[fx0])) * numden[Cx0];
			real n_y0 = 0.5 * (1.0 - copysign(1.0, uy[fy0])) * numden[C] + 0.5 * (1.0 + copysign(1.0, uy[fy0])) * numden[Cy0];
			real n_z0 = 0.5 * (1.0 - copysign(1.0, uz[fz0])) * numden[C] + 0.5 * (1.0 + copysign(1.0, uz[fz0])) * numden[Cz0];
			
			real nx = dt / dom->dx * (ux[fx1] * n_x1 - ux[fx0] * n_x0);
			real ny = dt / dom->dy * (uy[fy1] * n_y1 - uy[fy0] * n_y0);
			real nz = dt / dom->dz * (uz[fz1] * n_z1 - uz[fz0] * n_z0);
			numden1[C] = numden[C] - nx - ny - nz;
		}
	}
}
/*
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
*/
/*
__global__ void kernel_fx_coupling_numden_generate(real *x_nd, real *nd, dom_struct *dom)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	int tk = blockIdx.y * blockDim.y + threadIdx.y;
	for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
		if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
			int C  = i       + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b;
			int C0 = (i - 1) + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b;
			x_nd[C] = (nd[C] + nd[C0])/2.0;
		}
	}
}

__global__ void kernel_fy_coupling_numden_generate(real *y_nd, real *nd, dom_struct *dom)
{
	int tk = blockIdx.x * blockDim.x + threadIdx.x;
	int ti = blockIdx.y * blockDim.y + threadIdx.y;
	for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
		if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
			int C  = ti + j      *dom->Gfy._s1b + tk*dom->Gfy._s2b;
			int C0 = ti + (j - 1)*dom->Gfy._s1b + tk*dom->Gfy._s2b;
			y_nd[C] = (nd[C] + nd[C0])/2.0;
		}
	}
}
*/

__global__ void kernel_fz_coupling_numden_generate(real *z_nd, real *nd, dom_struct *dom)
{
	int ti = blockIdx.x * blockDim.x + threadIdx.x;
	int tj = blockIdx.y * blockDim.y + threadIdx.y;
	if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
		for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
			int C  = ti + tj*dom->Gfz._s1b + k      *dom->Gfz._s2b;
			int C1 = ti + tj*dom->Gcc._s1b + k      *dom->Gcc._s2b;
			int C0 = ti + tj*dom->Gcc._s1b + (k - 1)*dom->Gcc._s2b;
			z_nd[C] = (nd[C1] + nd[C0])/2.0;
		}
	}
}

__global__ void kernel_march_concentration(real dt,
										   dom_struct *dom,
										   real *conc,
										   real *conc1,
										   real *ux,
										   real *uy,
										   real *uz,
										   real *mdot,
										   real *numd,
										   real ccdiff)
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
			
			real c_x1 = 0.5 * (1.0 - copysign(1.0, ux[fx1])) * conc[Cx1] + 0.5 * (1.0 + copysign(1.0, ux[fx1])) * conc[C];
			real c_y1 = 0.5 * (1.0 - copysign(1.0, uy[fy1])) * conc[Cy1] + 0.5 * (1.0 + copysign(1.0, uy[fy1])) * conc[C];
			real c_z1 = 0.5 * (1.0 - copysign(1.0, uz[fz1])) * conc[Cz1] + 0.5 * (1.0 + copysign(1.0, uz[fz1])) * conc[C];
			
			real c_x0 = 0.5 * (1.0 - copysign(1.0, ux[fx0])) * conc[C] + 0.5 * (1.0 + copysign(1.0, ux[fx0])) * conc[Cx0];
			real c_y0 = 0.5 * (1.0 - copysign(1.0, uy[fy0])) * conc[C] + 0.5 * (1.0 + copysign(1.0, uy[fy0])) * conc[Cy0];
			real c_z0 = 0.5 * (1.0 - copysign(1.0, uz[fz0])) * conc[C] + 0.5 * (1.0 + copysign(1.0, uz[fz0])) * conc[Cz0];
			
			real convec_x = dt / dom->dx * (ux[fx1] * c_x1 - ux[fx0] * c_x0);
			real convec_y = dt / dom->dy * (uy[fy1] * c_y1 - uy[fy0] * c_y0);
			real convec_z = dt / dom->dz * (uz[fz1] * c_z1 - uz[fz0] * c_z0);
			
			real diffu_x = dt * ccdiff * (conc[Cx0] - 2.0 * conc[C] + conc[Cx1]) / dom->dx / dom->dx;
			real diffu_y = dt * ccdiff * (conc[Cy0] - 2.0 * conc[C] + conc[Cy1]) / dom->dy / dom->dy;
			real diffu_z = dt * ccdiff * (conc[Cz0] - 2.0 * conc[C] + conc[Cz1]) / dom->dz / dom->dz;
			
			conc1[C] = conc[C] - convec_x - convec_y - convec_z + diffu_x + diffu_y + diffu_z - dt * numd[C] * mdot[C];
		}
	}
}

/*
__global__ void kernel_concentration_update_concen(dom_struct *dom, real *conc, real *conc1)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
	int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;
	
	if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
		for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
			int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
			conc[C] = conc1[C];
		}
	}
}
*/

__global__ void kernel_numden_inner_copy(dom_struct *dom, real *numd, real *numd_tmp)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int ti = i + DOM_BUF;
	int tj = j + DOM_BUF;
	
	if(ti < dom->Gcc._ie && tj < dom->Gcc._je) {
		for(int tk = dom->Gcc._ks; tk < dom->Gcc._ke; tk++) {
			int k = tk - DOM_BUF;
			int c = i + j*dom->Gcc._s1 + k*dom->Gcc._s2;
			int tc = ti + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
			numd_tmp[c] = numd[tc];
		}
	}
}

__global__ void kernel_compute_mdot(dom_struct *dom,
                                    real *numd,
                                    real *conc,
                                    real *dia,
                                    real *mdot,
                                    real scale,
                                    real ccdiss,
                                    real ccdiff,
                                    real nu)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	int tk = blockIdx.y * blockDim.y + threadIdx.y;
	
	if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
		for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
			
			int C   = i       + tj      *dom->Gcc._s1b + tk      *dom->Gcc._s2b;
			//int fx0 = i       + tj      *dom->Gfx._s1b + tk      *dom->Gfx._s2b;
			//int fx1 = (i + 1) + tj      *dom->Gfx._s1b + tk      *dom->Gfx._s2b;
			//int fy0 = i       + tj      *dom->Gfy._s1b + tk      *dom->Gfy._s2b;
			//int fy1 = i       + (tj + 1)*dom->Gfy._s1b + tk      *dom->Gfy._s2b;
			//int fz0 = i       + tj      *dom->Gfz._s1b + tk      *dom->Gfz._s2b;
			//int fz1 = i       + tj      *dom->Gfz._s1b + (tk + 1)*dom->Gfz._s2b;
			
			if(numd[C] > 0) {
				real Dia = dia[C];
				real tervel = scale * Dia * Dia;
				real Re = tervel * Dia / nu;
				real Sc = nu / ccdiff;
				real Sh = 2.0 + 0.6 * pow(Re, 0.5) * pow(Sc, 1.0/3.0);
				real h = Sh * ccdiff / Dia;
				real A = PI * Dia * Dia;
				mdot[C] = A * h * (conc[C] - ccdiss);
			} else {
				mdot[C] = 0;
			}
			//printf("i==%d, j==%d, k==%d, dia==%f, mdot==%f\n", i, tj, tk, dia[C], mdot[C]);
		}
	}
}

// pressure(cell-centered value); bottom; dirichlet
__global__ void BC_p_B_D(real *p, dom_struct *dom, real bc)
{
	int ti = blockDim.x*blockIdx.x + threadIdx.x;
	int tj = blockDim.y*blockIdx.y + threadIdx.y;
	
	int s1b = dom->Gcc._s1b;
	int s2b = dom->Gcc._s2b;
	
	if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
		p[ti + tj*s1b + dom->Gcc._ksb*s2b] = 2 * bc - p[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// pressure(cell-centered value); top; dirichlet
__global__ void BC_p_T_D(real *p, dom_struct *dom, real bc)
{
	int ti = blockDim.x*blockIdx.x + threadIdx.x;
	int tj = blockDim.y*blockIdx.y + threadIdx.y;
	
	int s1b = dom->Gcc._s1b;
	int s2b = dom->Gcc._s2b;
	
	if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
		p[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = 2 * bc - p[ti + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

__global__ void kernel_compute_bubble_diameter(dom_struct *dom,
                                               real *mas,
                                               real *numd,
                                               real *dia,
                                               real rho_fluid,
                                               real pre_a,
                                               real rho_a,
                                               real gravacc)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	int tk = blockIdx.y * blockDim.y + threadIdx.y;
	if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
		for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
			
			int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
			
			if(numd[C] > 0) {
				real h = ((real)(dom->Gcc._keb - tk) - 1.5) * dom->dz;
				real rho = rho_a * (1.0 + rho_fluid * gravacc * h / pre_a);
				real vol = mas[C] / numd[C] / rho;
				dia[C] = pow(6.0 * vol / PI, 1.0/3.0);
			} else {
				dia[C] = 0;
			}
			//if(i==32 && tj==32 && tk<10) printf("compute bub dia: k==%d, diafz==%f\n", tk, dia[C]);
		}
	}
}

__global__ void kernel_compute_bubble_diameterfz(dom_struct *dom,
                                                 real *mas,
                                                 real *numfz,
                                                 real *diafz,
                                                 real rho_fluid,
                                                 real pre_a,
                                                 real rho_a,
                                                 real gravacc)
{
	int ti = blockIdx.x * blockDim.x + threadIdx.x;
	int tj = blockIdx.y * blockDim.y + threadIdx.y;
	if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
		for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
			
			int C  = ti + tj*dom->Gfz._s1b + k      *dom->Gfz._s2b;
			int C0 = ti + tj*dom->Gcc._s1b + (k - 1)*dom->Gcc._s2b;
			int C1 = ti + tj*dom->Gcc._s1b + k      *dom->Gcc._s2b;
			
			real mass_fz = 0.5 * (mas[C0] + mas[C1]);
			
			if(mass_fz > 0 && numfz[C] > 0) {
				real h = ((real)(dom->Gcc._keb - k) - 2.0) * dom->dz;
				real rho = rho_a * (1.0 + rho_fluid * gravacc * h / pre_a);
				real vol = mass_fz / numfz[C] / rho;
				diafz[C] = pow(6.0 * vol / PI, 1.0/3.0);
			} else {
				diafz[C] = 0;
			}
			//if(ti==32 && tj==32 && k<10) printf("compute bub diafz: k==%d, diafz==%f\n", k, diafz[C]);
		}
	}
}

__global__ void kernel_forcing_add_z_field_bubble(real scale,
                                                  real *ndfz,
                                                  real *diafz,
                                                  real *fz,
                                                  dom_struct *dom)
{
	int ti = blockIdx.x * blockDim.x + threadIdx.x;
	int tj = blockIdx.y * blockDim.y + threadIdx.y;
	
	if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
		for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
			int C = ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b;
			fz[C] += scale * ndfz[C] * diafz[C] * diafz[C] * diafz[C];
			
			//if(ti==32 && tj==32 && k<10)
			//printf("forcing: k==%d, fz==%f, nd==%f, diafz==%f\n", k, fz[C], ndfz[C], diafz[C]);
		}
	}
}

__global__ void kernel_inner_scalarfield_update_x(dom_struct *dom,
                                                  real *left,
                                                  real *right)
{
	int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
	int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;
	
	if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
		for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
			int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
			left[C] = right[C];
		}
	}
}

__global__ void kernel_march_bubblemass(real dt,
										dom_struct *dom,
										real *bubm,
										real *bubm1,
										real *ux,
										real *uy,
										real *uz,
										real *numd,
										real *mdot)
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
			/*
			real ux_c = 0.5 * (ux[fx0] + ux[fx1]);
			real uy_c = 0.5 * (uy[fy0] + uy[fy1]);
			real uz_c = 0.5 * (uz[fz0] + uz[fz1]);
			
			real pv_x = 0.5 * (1.0 - copysign(1.0, ux_c)) * (bubm[Cx1] - bubm[C]) + 0.5 * (1.0 + copysign(1.0, ux_c)) * (bubm[C] - bubm[Cx0]);
			real pv_y = 0.5 * (1.0 - copysign(1.0, uy_c)) * (bubm[Cy1] - bubm[C]) + 0.5 * (1.0 + copysign(1.0, uy_c)) * (bubm[C] - bubm[Cy0]);
			real pv_z = 0.5 * (1.0 - copysign(1.0, uz_c)) * (bubm[Cz1] - bubm[C]) + 0.5 * (1.0 + copysign(1.0, uz_c)) * (bubm[C] - bubm[Cz0]);
			
			pv_x = pv_x / dom->dx * ux_c * dt;
			pv_y = pv_y / dom->dy * uy_c * dt;
			pv_z = pv_z / dom->dz * uz_c * dt;
			*/
			
			real m_x1 = 0.5 * (1.0 - copysign(1.0, ux[fx1])) * bubm[Cx1] + 0.5 * (1.0 + copysign(1.0, ux[fx1])) * bubm[C];
			real m_y1 = 0.5 * (1.0 - copysign(1.0, uy[fy1])) * bubm[Cy1] + 0.5 * (1.0 + copysign(1.0, uy[fy1])) * bubm[C];
			real m_z1 = 0.5 * (1.0 - copysign(1.0, uz[fz1])) * bubm[Cz1] + 0.5 * (1.0 + copysign(1.0, uz[fz1])) * bubm[C];
			
			real m_x0 = 0.5 * (1.0 - copysign(1.0, ux[fx0])) * bubm[C] + 0.5 * (1.0 + copysign(1.0, ux[fx0])) * bubm[Cx0];
			real m_y0 = 0.5 * (1.0 - copysign(1.0, uy[fy0])) * bubm[C] + 0.5 * (1.0 + copysign(1.0, uy[fy0])) * bubm[Cy0];
			real m_z0 = 0.5 * (1.0 - copysign(1.0, uz[fz0])) * bubm[C] + 0.5 * (1.0 + copysign(1.0, uz[fz0])) * bubm[Cz0];
			
			real mx = dt / dom->dx * (ux[fx1] * m_x1 - ux[fx0] * m_x0);
			real my = dt / dom->dy * (uy[fy1] * m_y1 - uy[fy0] * m_y0);
			real mz = dt / dom->dz * (uz[fz1] * m_z1 - uz[fz0] * m_z0);
			
			bubm1[C] = bubm[C] - mx - my - mz + numd[C] * mdot[C] * dt;
		}
	}
}
