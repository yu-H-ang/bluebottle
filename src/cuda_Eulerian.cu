#include "cuda_Eulerian.h"
#include "cuda_bluebottle.h"
#include "bluebottle.h"

#include <cuda.h>
#include <helper_cuda.h>

extern "C"
void cuda_numberdensity_malloc(void)
{
	// allocate device memory on host
	_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_nextnumden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_u_p = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_v_p = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_w_p = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	// allocate device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		checkCudaErrors(cudaMalloc((void**) &(_numden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _numden_next is to store _numden in next time step, only exist on deivce.
		checkCudaErrors(cudaMalloc((void**) &(_nextnumden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_u_p[dev]), sizeof(real) * dom[dev].Gfx.s3b));
		gpumem += dom[dev].Gfx.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_v_p[dev]), sizeof(real) * dom[dev].Gfy.s3b));
		gpumem += dom[dev].Gfy.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_w_p[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
	}
}

extern "C"
void cuda_numberdensity_push(void)
{
	int i, j, k;          // iterators
	int ii, jj, kk;       // helper iterators
	int C, CC;            // cell references
	
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		// set up host working arrays for subdomain copy from host to device
		real *nn = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *uu_p = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
		real *vv_p = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
		real *ww_p = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
		
		// number density
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					nn[CC] = numden[C];
				}
			}
		}
		
		// u_p
		for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
			for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
				for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
					ii = i - dom[dev].Gfx.isb;
					jj = j - dom[dev].Gfx.jsb;
					kk = k - dom[dev].Gfx.ksb;
					C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
					CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
					uu_p[CC] = u_p[C];
				}
			}
		}
		// v_p
		for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
			for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
				for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
					ii = i - dom[dev].Gfy.isb;
					jj = j - dom[dev].Gfy.jsb;
					kk = k - dom[dev].Gfy.ksb;
					C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
					CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
					vv_p[CC] = v_p[C];
				}
			}
		}
		// w_p
		for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
			for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
				for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
					ii = i - dom[dev].Gfz.isb;
					jj = j - dom[dev].Gfz.jsb;
					kk = k - dom[dev].Gfz.ksb;
					C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
					CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
					ww_p[CC] = w_p[C];
				}
			}
		}
		
		// copy from host to device
		checkCudaErrors(cudaMemcpy(_numden[dev], nn, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(_u_p[dev], uu_p, sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(_v_p[dev], vv_p, sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(_w_p[dev], ww_p, sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
		
		// free host subdomain working arrays
		free(nn);
		free(uu_p);
		free(vv_p);
		free(ww_p);
	}
}

void cuda_numberdensity_pull(void)
{
	// copy device data to host
	#pragma omp parallel num_threads(nsubdom)
	{
		int i, j, k;          // iterators
		int ii, jj, kk;       // helper iterators
		int C, CC;            // cell references
		
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));

		// host working arrays for subdomain copy from device to host
		real *nn = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *uu_p = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
		real *vv_p = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
		real *ww_p = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
		
		// copy from device to host
		checkCudaErrors(cudaMemcpy(nn, _numden[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(uu_p, _u_p[dev], sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(vv_p, _v_p[dev], sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(ww_p, _w_p[dev], sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));
		
		// fill in apropriate subdomain (copy back ghost cells)
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					numden[C] = nn[CC];
				}
			}
		}
		// u
		for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
			for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
				for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
					ii = i - dom[dev].Gfx.isb;
					jj = j - dom[dev].Gfx.jsb;
					kk = k - dom[dev].Gfx.ksb;
					C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
					CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
					u_p[C] = uu_p[CC];
				}
			}
		}
		// v
		for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
			for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
				for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
					ii = i - dom[dev].Gfy.isb;
					jj = j - dom[dev].Gfy.jsb;
					kk = k - dom[dev].Gfy.ksb;
					C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
					CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
					v_p[C] = vv_p[CC];
				}
			}
		}
		// w
		for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
			for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
				for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
					ii = i - dom[dev].Gfz.isb;
					jj = j - dom[dev].Gfz.jsb;
					kk = k - dom[dev].Gfz.ksb;
					C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
					CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
					w_p[C] = ww_p[CC];
				}
			}
		}
		
		// free host subdomain working arrays
		free(nn);
		free(uu_p);
		free(vv_p);
		free(ww_p);
	}
}

extern "C"
void cuda_numberdensity_free(void)
{
	// free device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		checkCudaErrors(cudaFree(_numden[dev]));
		checkCudaErrors(cudaFree(_nextnumden[dev]));
	}
	
	// free device memory on host
	free(_numden);
	free(_u_p);
	free(_v_p);
	free(_w_p);
}

extern "C"
void cuda_numberdensity_BC(void)
{
	// CPU threading for multi-GPU
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int threads_y = 0;
		int threads_z = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		int blocks_z = 0;
		
		// check whether each subdomain boundary (E, W, N, S, T, B) is
		// an external boundary
		// *************************************************************
		if(dom[dev].W == -1) {
			// set up kernel call
			if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
				threads_y = dom[dev].Gcc.jnb;
			else
				threads_y = MAX_THREADS_DIM;
				
			if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
				threads_z = dom[dev].Gcc.knb;
			else
				threads_z = MAX_THREADS_DIM;
			
			blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
			blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
			
			dim3 dimBlocks_n(threads_y, threads_z);
			dim3 numBlocks_n(blocks_y, blocks_z);
			
			// apply BC to numberdensity field for this face
			switch(numdenBC.nW) {
				case PERIODIC:
					BC_p_W_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_W_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
			}
		}
		// *************************************************************
		if(dom[dev].E == -1) {
			// set up kernel call
			if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
				threads_y = dom[dev].Gcc.jnb;
			else
				threads_y = MAX_THREADS_DIM;
				
			if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
				threads_z = dom[dev].Gcc.knb;
			else
				threads_z = MAX_THREADS_DIM;
			
			blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
			blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
			
			dim3 dimBlocks_n(threads_y, threads_z);
			dim3 numBlocks_n(blocks_y, blocks_z);
			// apply BC to numberdensity field for this face
			switch(numdenBC.nE) {
				case PERIODIC:
					BC_p_E_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_E_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
			}
		}
		// *************************************************************
		if(dom[dev].S == -1) {
			// set up kernel call
			if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
				threads_z = dom[dev].Gcc.knb;
			else
				threads_z = MAX_THREADS_DIM;

			if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
				threads_x = dom[dev].Gcc.inb;
			else
				threads_x = MAX_THREADS_DIM;

			blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
			blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

			dim3 dimBlocks_n(threads_z, threads_x);
			dim3 numBlocks_n(blocks_z, blocks_x);
			// apply BC to numberdensity field for this face
			switch(numdenBC.nS) {
				case PERIODIC:
					BC_p_S_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_S_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
			}
		}
		// *************************************************************
		if(dom[dev].N == -1) {
			// set up kernel call
			if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
				threads_z = dom[dev].Gcc.knb;
			else
				threads_z = MAX_THREADS_DIM;

			if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
				threads_x = dom[dev].Gcc.inb;
			else
				threads_x = MAX_THREADS_DIM;

			blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
			blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

			dim3 dimBlocks_n(threads_z, threads_x);
			dim3 numBlocks_n(blocks_z, blocks_x);

			// apply BC to numberdensity field for this face
			switch(numdenBC.nN) {
					case PERIODIC:
					BC_p_N_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
				break;
					case NEUMANN:
					BC_p_N_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
				break;
			}
		}
		// *************************************************************
		if(dom[dev].B == -1) {
			// set up kernel call
			if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
				threads_x = dom[dev].Gcc.inb;
			else
				threads_x = MAX_THREADS_DIM;

			if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
				threads_y = dom[dev].Gcc.jnb;
			else
				threads_y = MAX_THREADS_DIM;

			blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
			blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

			dim3 dimBlocks_n(threads_x, threads_y);
			dim3 numBlocks_n(blocks_x, blocks_y);
			// apply BC to numberdensity field for this face
			switch(numdenBC.nB) {
				case PERIODIC:
					BC_p_B_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_B_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
			}
		}
		// *************************************************************
		if(dom[dev].T == -1) {
			// set up kernel call
			if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
				threads_x = dom[dev].Gcc.inb;
			else
				threads_x = MAX_THREADS_DIM;

			if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
				threads_y = dom[dev].Gcc.jnb;
			else
				threads_y = MAX_THREADS_DIM;

			blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
			blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

			dim3 dimBlocks_n(threads_x, threads_y);
			dim3 numBlocks_n(blocks_x, blocks_y);
			// apply BC to numberdensity field for this face
			switch(numdenBC.nT) {
				case PERIODIC:
					BC_p_T_P<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_T_N<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
					break;
			}
		}
		// *************************************************************
		
		
	}
}

extern "C"
void cuda_numberdensity_march(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		int threads_y = 0;
		int threads_z = 0;
		int blocks_y = 0;
		int blocks_z = 0;
		
		if(dom[dev].Gcc._jn < MAX_THREADS_DIM)
			threads_y = dom[dev].Gcc._jn;
		else
			threads_y = MAX_THREADS_DIM;
		
		if(dom[dev].Gcc._kn < MAX_THREADS_DIM)
			threads_z = dom[dev].Gcc._kn;
		else
			threads_z = MAX_THREADS_DIM;
		
		blocks_y = (int)ceil((real) dom[dev].Gcc._jn / (real) threads_y);
		blocks_z = (int)ceil((real) dom[dev].Gcc._kn / (real) threads_z);
		
		dim3 dimBlocks_n(threads_y, threads_z);
		dim3 numBlocks_n(blocks_y, blocks_z);
		
		kernel_numberdensity_march<<<numBlocks_n, dimBlocks_n>>>(dt, _dom[dev], _numden[dev], _nextnumden[dev], _u[dev], _v[dev], _w_p[dev]);
		kernel_numberdensity_update_numden<<<numBlocks_n, dimBlocks_n>>>(_dom[dev], _numden[dev], _nextnumden[dev]);
	}
}


extern "C"
void cuda_numberdensity_particle_velz(void)
{
	// parallelize over CPU threads
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int threads_y = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		
		// z-component
		if(dom[dev].Gfz._inb < MAX_THREADS_DIM)
			threads_x = dom[dev].Gfz._inb;
		else
			threads_x = MAX_THREADS_DIM;
		if(dom[dev].Gfz._jnb < MAX_THREADS_DIM)
			threads_y = dom[dev].Gfz._jnb;
		else
			threads_y = MAX_THREADS_DIM;
		
		blocks_x = (int)ceil((real) dom[dev].Gfz._inb / (real) threads_x);
		blocks_y = (int)ceil((real) dom[dev].Gfz._jnb / (real) threads_y);
		dim3 dimBlocks_z(threads_x, threads_y);
		dim3 numBlocks_z(blocks_x, blocks_y);
		
		real drift_vel = -2.0/9.0*(rho_f - bubble_density)/mu*g.zm*bubble_radius*bubble_radius;
		//printf("terminal velocity = %f \n", drift_vel);
		
		kernel_numberdensity_particle_velz<<<numBlocks_z, dimBlocks_z>>>(drift_vel, _w_p[dev], _w[dev], _dom[dev]);
	}
}
