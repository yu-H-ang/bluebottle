#include "cuda_Eulerian.h"
#include "cuda_bluebottle.h"
#include "bluebottle.h"
#include "Eulerian.h"
#include "entrySearch.h"

#include <cuda.h>
#include <helper_cuda.h>

extern "C"
void cuda_Eulerian_push(void)
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
		//real *uu_p = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
		//real *vv_p = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
		//real *ww_p = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
		real *bubv = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *cc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		
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
		/*
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
		*/
		// bubble volume
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubv[CC] = bubvol[C];
				}
			}
		}
		// concentration
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					cc[CC] = concen[C];
				}
			}
		}
		
		// copy from host to device
		checkCudaErrors(cudaMemcpy(_numden[dev], nn, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		//checkCudaErrors(cudaMemcpy(_u_p[dev], uu_p, sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
		//checkCudaErrors(cudaMemcpy(_v_p[dev], vv_p, sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
		//checkCudaErrors(cudaMemcpy(_w_p[dev], ww_p, sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(_bubvol[dev], bubv, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(_concen[dev], cc, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		
		// free host subdomain working arrays
		free(nn);
		//free(uu_p);
		//free(vv_p);
		//free(ww_p);
		free(bubv);
		free(cc);
	}
}

void cuda_Eulerian_pull(void)
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
		//real *uu_p = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
		//real *vv_p = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
		//real *ww_p = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
		real *bubv = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *cc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		
		// copy from device to host
		checkCudaErrors(cudaMemcpy(nn, _numden[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		//checkCudaErrors(cudaMemcpy(uu_p, _u_p[dev], sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
		//checkCudaErrors(cudaMemcpy(vv_p, _v_p[dev], sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
		//checkCudaErrors(cudaMemcpy(ww_p, _w_p[dev], sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(bubv, _bubvol[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(cc, _concen[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		
		// fill in apropriate subdomain (copy back ghost cells)
		// numden
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
		/*
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
		*/
		// bubvol
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubvol[C] = bubv[CC];
				}
			}
		}
		// concen
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					concen[C] = cc[CC];
				}
			}
		}
		
		// free host subdomain working arrays
		free(nn);
		//free(uu_p);
		//free(vv_p);
		//free(ww_p);
		free(bubv);
		free(cc);
	}
}

extern "C"
void cuda_Eulerian_free(void)
{
	// free device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		checkCudaErrors(cudaFree(_numden[dev]));
		checkCudaErrors(cudaFree(_nextnumden[dev]));
		//checkCudaErrors(cudaFree(_u_p[dev]));
		//checkCudaErrors(cudaFree(_v_p[dev]));
		checkCudaErrors(cudaFree(_w_p[dev]));
		//checkCudaErrors(cudaFree(_f_x_coupling_numden[dev]));
		//checkCudaErrors(cudaFree(_f_y_coupling_numden[dev]));
		checkCudaErrors(cudaFree(_f_z_coupling_numden[dev]));
		
		checkCudaErrors(cudaFree(_bubvol[dev]));
		checkCudaErrors(cudaFree(_nextbubvol[dev]));
		checkCudaErrors(cudaFree(_bubdia[dev]));
		checkCudaErrors(cudaFree(_bubdiafz[dev]));
		
		checkCudaErrors(cudaFree(_concen[dev]));
		checkCudaErrors(cudaFree(_nextconcen[dev]));
		checkCudaErrors(cudaFree(_velmag[dev]));
		checkCudaErrors(cudaFree(_mdot[dev]));
	}
	
	// free device memory on host
	free(_numden);
	free(_nextnumden);
	//free(_u_p);
	//free(_v_p);
	free(_w_p);
	//free(_f_x_coupling_numden);
	//free(_f_y_coupling_numden);
	free(_f_z_coupling_numden);
	
	free(_bubvol);
	free(_nextbubvol);
	free(_bubdia);
	free(_bubdiafz);
	
	free(_concen);
	free(_nextconcen);
	free(_velmag);
	free(_mdot);
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
		// *********************************************************************
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
		// *********************************************************************
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
		// *********************************************************************
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
		// *********************************************************************
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
		// *********************************************************************
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
				case DIRICHLET:
					BC_p_B_D<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev], numdenBC.nBD);
					break;
			}
		}
		// *********************************************************************
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
				case DIRICHLET:
					BC_p_T_D<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev], numdenBC.nTD);
					break;
			}
		}
		// *********************************************************************
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
		
		kernel_march_numberdensity<<<numBlocks_n, dimBlocks_n>>>(dt,
		                                                         _dom[dev],
		                                                         _numden[dev],
		                                                         _nextnumden[dev],
		                                                         _u[dev],
		                                                         _v[dev],
		                                                         _w_p[dev]);
		
		kernel_inner_scalarfield_update_x<<<numBlocks_n, dimBlocks_n>>>(_dom[dev],
		                                                                _numden[dev],
		                                                                _nextnumden[dev]);
	}
}


extern "C"
void cuda_compute_particle_velz(void)
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
		
		// since bubble diameter is a field, here combine all the values in terminal velocity except bubble diameter
		real cons = -1.0/18.0*(rho_f - bubble_density)/mu*g.zm;
		
		kernel_numberdensity_particle_velz<<<numBlocks_z, dimBlocks_z>>>(cons, _w_p[dev], _w[dev], _bubdia[dev], _dom[dev]);
	}
}

extern "C"
void cuda_Eulerian_malloc(void)
{
	// allocate device memory on host
	_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_nextnumden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	//_u_p = (real**) malloc(nsubdom * sizeof(real*));
	//cpumem += nsubdom * sizeof(real*);
	//_v_p = (real**) malloc(nsubdom * sizeof(real*));
	//cpumem += nsubdom * sizeof(real*);
	_w_p = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	//_f_x_coupling_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	//_f_y_coupling_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_f_z_coupling_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	_bubvol = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_nextbubvol = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_bubdia = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_bubdiafz = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	_concen = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_nextconcen = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_velmag = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_mdot = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	// allocate device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		// number density equation
		checkCudaErrors(cudaMalloc((void**) &(_numden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _numden_next is to store _numden in next time step, only exist on deivce.
		checkCudaErrors(cudaMalloc((void**) &(_nextnumden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		//checkCudaErrors(cudaMalloc((void**) &(_u_p[dev]), sizeof(real) * dom[dev].Gfx.s3b));
		//gpumem += dom[dev].Gfx.s3b * sizeof(real);
		//checkCudaErrors(cudaMalloc((void**) &(_v_p[dev]), sizeof(real) * dom[dev].Gfy.s3b));
		//gpumem += dom[dev].Gfy.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_w_p[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		//checkCudaErrors(cudaMalloc((void**) &(_f_x_coupling_numden[dev]), sizeof(real) * dom[dev].Gfx.s3b));
		//gpumem += dom[dev].Gfx.s3b * sizeof(real);
		//checkCudaErrors(cudaMalloc((void**) &(_f_y_coupling_numden[dev]), sizeof(real) * dom[dev].Gfy.s3b));
		//gpumem += dom[dev].Gfy.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_f_z_coupling_numden[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		
		// bubble volume
		checkCudaErrors(cudaMalloc((void**) &(_bubvol[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_nextbubvol[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_bubdia[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_bubdiafz[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		
		// concentration
		checkCudaErrors(cudaMalloc((void**) &(_concen[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _nextconcen is to store _concen in next time step, only exist on deivce.
		checkCudaErrors(cudaMalloc((void**) &(_nextconcen[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _velmag is to store the velocity magnitude of the fluid flow,
		// which is used in the source term of mass transfer, only exist on device.
		checkCudaErrors(cudaMalloc((void**) &(_velmag[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		checkCudaErrors(cudaMalloc((void**) &(_mdot[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
	}
}

extern "C"
void cuda_concentration_BC(void)
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
		// *********************************************************************
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
			
			// apply BC to concentration field for this face
			switch(concenBC.nW) {
				case PERIODIC:
					BC_p_W_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_W_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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
			// apply BC to concentration field for this face
			switch(concenBC.nE) {
				case PERIODIC:
					BC_p_E_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_E_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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
			// apply BC to concentration field for this face
			switch(concenBC.nS) {
				case PERIODIC:
					BC_p_S_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_S_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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

			// apply BC to concentration field for this face
			switch(concenBC.nN) {
					case PERIODIC:
					BC_p_N_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
				break;
					case NEUMANN:
					BC_p_N_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
				break;
			}
		}
		// *********************************************************************
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
			// apply BC to concentration field for this face
			switch(concenBC.nB) {
				case PERIODIC:
					BC_p_B_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_B_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_B_D<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev], concenBC.nBD);
					break;
			}
		}
		// *********************************************************************
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
			// apply BC to concentration field for this face
			switch(concenBC.nT) {
				case PERIODIC:
					BC_p_T_P<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_T_N<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_T_D<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev], concenBC.nTD);
					break;
			}
		}
		// *********************************************************************
	}
}

extern "C"
void cuda_compute_coupling_forcing(void)
{
	// parallelize over CPU threads
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

		// x-component
		if(dom[dev].Gfx._jnb < MAX_THREADS_DIM)
			threads_y = dom[dev].Gfx._jnb;
		else
			threads_y = MAX_THREADS_DIM;

		if(dom[dev].Gfx._knb < MAX_THREADS_DIM)
			threads_z = dom[dev].Gfx._knb;
		else
			threads_z = MAX_THREADS_DIM;

		blocks_y = (int)ceil((real) dom[dev].Gfx._jnb / (real) threads_y);
		blocks_z = (int)ceil((real) dom[dev].Gfx._knb / (real) threads_z);

		dim3 dimBlocks_x(threads_y, threads_z);
		dim3 numBlocks_x(blocks_y, blocks_z);

		// y-component
		if(dom[dev].Gfy._knb < MAX_THREADS_DIM)
			threads_z = dom[dev].Gfy._knb;
		else
			threads_z = MAX_THREADS_DIM;

		if(dom[dev].Gfy._inb < MAX_THREADS_DIM)
			threads_x = dom[dev].Gfy._inb;
		else
			threads_x = MAX_THREADS_DIM;

		blocks_z = (int)ceil((real) dom[dev].Gfy._knb / (real) threads_z);
		blocks_x = (int)ceil((real) dom[dev].Gfy._inb / (real) threads_x);

		dim3 dimBlocks_y(threads_z, threads_x);
		dim3 numBlocks_y(blocks_z, blocks_x);
    
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

		//need to do a linear interpolation to calculate the number density field on cell faces
		kernel_fz_coupling_numden_generate<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _numden[dev], _dom[dev]);

		//calculate numden[] on boundaries, this part will fail when multi-device domain decomposition is used.
		if(dom[dev].W == -1) {
			switch(numdenBC.nW) {
				case PERIODIC:
				BC_w_W_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].E == -1) {
			switch(numdenBC.nE) {
				case PERIODIC:
				BC_w_E_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].S == -1) {
			switch(numdenBC.nS) {
				case PERIODIC:
				BC_w_S_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].N == -1) {
			switch(numdenBC.nN) {
				case PERIODIC:
				BC_w_N_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].B == -1) {
			switch(numdenBC.nB) {
				case PERIODIC:
				BC_w_B_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
				case DIRICHLET:
				BC_w_B_D<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev], numdenBC.nBD);
				break;
				case NEUMANN:
				BC_w_B_N<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].T == -1) {
			switch(numdenBC.nT) {
				case PERIODIC:
				BC_w_T_P<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
				case DIRICHLET:
				BC_w_T_D<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev], numdenBC.nTD);
				break;
				case NEUMANN:
				BC_w_T_N<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev], _dom[dev]);
				break;
			}
		}
		
		// reset forcing arrays
		forcing_reset_x<<<numBlocks_x, dimBlocks_x>>>(_f_x[dev], _dom[dev]);
		forcing_reset_y<<<numBlocks_y, dimBlocks_y>>>(_f_y[dev], _dom[dev]);
		forcing_reset_z<<<numBlocks_z, dimBlocks_z>>>(_f_z[dev], _dom[dev]);
		
		// now add in the forcing
		real forcing_scale = -1.0/6.0*PI*(rho_f - bubble_density)*g.zm/rho_f;
		
		kernel_forcing_add_z_field_bubble<<<numBlocks_z, dimBlocks_z>>>(forcing_scale,
		                                                                _f_z_coupling_numden[dev],
		                                                                _bubdia[dev],
		                                                                _f_z[dev],
		                                                                _dom[dev]);
	}
}

extern "C"
void cuda_concentration_march(void)
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
		
		dim3 dimBlocks_c(threads_y, threads_z);
		dim3 numBlocks_c(blocks_y, blocks_z);
		
		// march the concentration equation
		kernel_march_concentration<<<numBlocks_c, dimBlocks_c>>>(dt,
                                                                 _dom[dev],
                                                                 _concen[dev],
                                                                 _nextconcen[dev],
                                                                 _u[dev],
                                                                 _v[dev],
                                                                 _w[dev],
                                                                 _mdot[dev],
                                                                 _numden[dev],
                                                                 concen_diff);
		
		// update concentration field
		kernel_inner_scalarfield_update_x<<<numBlocks_c, dimBlocks_c>>>(_dom[dev],
		                                                                 _concen[dev],
		                                                                 _nextconcen[dev]);
	}
}

extern "C"
void cuda_numberdensity_compute_totalnumden(void)
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
		
		// N is the number of inner cells
		int N = dom[dev].Gcc.s3;
		
		// create temporary storage for reduction algorithms
		real *_numden_tmp;
		checkCudaErrors(cudaMalloc((void**) &_numden_tmp, sizeof(real) * N));
		gpumem += sizeof(real) * N;
		
		// set up kernel call
		if(dom[dev].Gcc.in < MAX_THREADS_DIM)
			threads_x = dom[dev].Gcc.in;
		else
			threads_x = MAX_THREADS_DIM;
			
		if(dom[dev].Gcc.jn < MAX_THREADS_DIM)
			threads_y = dom[dev].Gcc.jn;
		else
			threads_y = MAX_THREADS_DIM;
		
		blocks_x = (int)ceil((real) dom[dev].Gcc.in / (real) threads_x);
		blocks_y = (int)ceil((real) dom[dev].Gcc.jn / (real) threads_y);
		
		dim3 dimBlocks_z(threads_x, threads_y);
		dim3 numBlocks_z(blocks_x, blocks_y);
		
		// calculate summation
		kernel_numden_inner_copy<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
		                                                       _numden[dev],
		                                                       _numden_tmp);
		totalnumden = sum_entries(N, _numden_tmp);
		
		// clean up
		checkCudaErrors(cudaFree(_numden_tmp));
		}
}

extern "C"
void cuda_bubblevolume_BC(void)
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
		// *********************************************************************
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
			
			dim3 dimBlocks(threads_y, threads_z);
			dim3 numBlocks(blocks_y, blocks_z);
			
			// apply BC to bubble volume field for this face
			switch(bubvolBC.nW) {
				case PERIODIC:
					BC_p_W_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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
			
			dim3 dimBlocks(threads_y, threads_z);
			dim3 numBlocks(blocks_y, blocks_z);
			// apply BC to bubble volume for this face
			switch(bubvolBC.nE) {
				case PERIODIC:
					BC_p_E_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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

			dim3 dimBlocks(threads_z, threads_x);
			dim3 numBlocks(blocks_z, blocks_x);
			// apply BC to bubble volume for this face
			switch(bubvolBC.nS) {
				case PERIODIC:
					BC_p_S_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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

			dim3 dimBlocks(threads_z, threads_x);
			dim3 numBlocks(blocks_z, blocks_x);

			// apply BC to bubble volume for this face
			switch(bubvolBC.nN) {
					case PERIODIC:
					BC_p_N_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
				break;
			}
		}
		// *********************************************************************
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

			dim3 dimBlocks(threads_x, threads_y);
			dim3 numBlocks(blocks_x, blocks_y);
			// apply BC to bubble volume for this face
			switch(bubvolBC.nB) {
				case PERIODIC:
					BC_p_B_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_B_D<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev], bubvolBC.nBD);
					break;
				case NEUMANN:
					BC_p_B_N<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
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

			dim3 dimBlocks(threads_x, threads_y);
			dim3 numBlocks(blocks_x, blocks_y);
			// apply BC to bubble volume for this face
			switch(bubvolBC.nT) {
				case PERIODIC:
					BC_p_T_P<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_T_D<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev], bubvolBC.nTD);
					break;
				case NEUMANN:
					BC_p_T_N<<<numBlocks, dimBlocks>>>(_bubvol[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
		
		
	}
}

extern "C"
void cuda_compute_bubble_diameter(void)
{
	// parallelize over CPU threads
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
		
		// generate bubdia using bubvol(cell-centered field)
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
		
		dim3 dimBlocks(threads_y, threads_z);
		dim3 numBlocks(blocks_y, blocks_z);
		
		kernel_compute_bubble_diameter<<<numBlocks, dimBlocks>>>(_dom[dev],
		                                                         _bubvol[dev],
		                                                         _bubdia[dev]);
		
		// interpolate bubdia to generate a z-face-centered field, which is needed in terminal velocity
		if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
			threads_x = dom[dev].Gfz.inb;
		else
			threads_x = MAX_THREADS_DIM;
		
		if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
			threads_y = dom[dev].Gfz.jnb;
		else
			threads_y = MAX_THREADS_DIM;
		
		blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
		blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
		
		dim3 dimBlocks_z(threads_x, threads_y);
		dim3 numBlocks_z(blocks_x, blocks_y);
		
		kernel_compute_bubble_diameterfz<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
		                                                               _bubdia[dev],
		                                                               _bubdiafz[dev]);
		
		if(dom[dev].W == -1) {
			switch(bubvolBC.nW) {
				case PERIODIC:
				BC_w_W_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].E == -1) {
			switch(bubvolBC.nE) {
				case PERIODIC:
				BC_w_E_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].S == -1) {
			switch(bubvolBC.nS) {
				case PERIODIC:
				BC_w_S_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].N == -1) {
			switch(bubvolBC.nN) {
				case PERIODIC:
				BC_w_N_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].B == -1) {
			switch(bubvolBC.nB) {
				case PERIODIC:
				BC_w_B_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
				case DIRICHLET:
				BC_w_B_D<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev], pow(6.0*bubvolBC.nBD/PI, 1.0/3.0));
				break;
				case NEUMANN:
				BC_w_B_N<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
		if(dom[dev].T == -1) {
			switch(bubvolBC.nT) {
				case PERIODIC:
				BC_w_T_P<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
				case DIRICHLET:
				BC_w_T_D<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev], pow(6.0*bubvolBC.nTD/PI, 1.0/3.0));
				break;
				case NEUMANN:
				BC_w_T_N<<<numBlocks_z, dimBlocks_z>>>(_bubdiafz[dev], _dom[dev]);
				break;
			}
		}
	}
}

extern "C"
void cuda_compute_mdot(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		checkCudaErrors(cudaSetDevice(dev + dev_start));
		
		int threads_y = 0;
		int threads_z = 0;
		int blocks_y = 0;
		int blocks_z = 0;
		
		if(dom[dev].Gcc._jnb < MAX_THREADS_DIM)
			threads_y = dom[dev].Gcc._jnb;
		else
			threads_y = MAX_THREADS_DIM;
		
		if(dom[dev].Gcc._knb < MAX_THREADS_DIM)
			threads_z = dom[dev].Gcc._knb;
		else
			threads_z = MAX_THREADS_DIM;
		
		blocks_y = (int)ceil((real) dom[dev].Gcc._jnb / (real) threads_y);
		blocks_z = (int)ceil((real) dom[dev].Gcc._knb / (real) threads_z);
		
		dim3 dimBlocks(threads_y, threads_z);
		dim3 numBlocks(blocks_y, blocks_z);
		
		kernel_compute_mdot<<<numBlocks, dimBlocks>>>(_dom[dev],
													  _concen[dev],
													  _bubdia[dev],
													  _mdot[dev],
													  _velmag[dev],
													  _u[dev],
													  _v[dev],
													  _w[dev],
													  concen_diss,
													  concen_diff,
													  nu);
	}
}

extern "C"
void cuda_bubblevolume_march(void)
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
		
		dim3 dimBlocks(threads_y, threads_z);
		dim3 numBlocks(blocks_y, blocks_z);
		
		// march the bubble volume equation
		kernel_march_bubblevolume<<<numBlocks, dimBlocks>>>(dt,
                                                            _dom[dev],
                                                            _bubvol[dev],
                                                            _nextbubvol[dev],
                                                            _u[dev],
                                                            _v[dev],
                                                            _w_p[dev],
                                                            _mdot[dev],
                                                            bubble_density);
		
		// update concentration field
		kernel_inner_scalarfield_update_x<<<numBlocks, dimBlocks>>>(_dom[dev],
		                                                            _bubvol[dev],
		                                                            _nextbubvol[dev]);
	}
}
