#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/sort.h>

#include "Eulerian_cuda.h"
#include "cuda_bluebottle.h"
#include "bluebottle.h"
#include "Eulerian.h"
#include "entrySearch.h"

extern "C"
void cuda_Eulerian_push(void)
{
	int i, j, k;          // iterators
	int ii, jj, kk;       // helper iterators
	int C, CC;            // cell references
	
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		// set up host working arrays for subdomain copy from host to device
		real *nn = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *bubm = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *cc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		//real *uu_p = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
		//real *vv_p = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
		//real *ww_b = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
		real *bden = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *ccsat = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *bubgen_num = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *bubgen_mas = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		
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
		
		// bubble generator
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubgen_num[CC] = BGndot[C];
				}
			}
		}
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubgen_mas[CC] = BGmdot[C];
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
		// bubble mass
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubm[CC] = bubmas[C];
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
		// bubble density
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bden[CC] = bubden[C];
				}
			}
		}
		// saturation concentration at different depth
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					ccsat[CC] = concen_sat[C];
				}
			}
		}
		
		// face-centered bubble density
		if(gravity_direction == GRAVITY_X) {
			
		}
		else if(gravity_direction == GRAVITY_Y) {
			
		}
		else if(gravity_direction == GRAVITY_Z) {
			
			real *bdenf = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
			for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
				for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
					for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
						ii = i - dom[dev].Gfz.isb;
						jj = j - dom[dev].Gfz.jsb;
						kk = k - dom[dev].Gfz.ksb;
						C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
						CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
						bdenf[CC] = bubden_face[C];
					}
				}
			}
			(cudaMemcpy(_bubden_face[dev], bdenf, sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
			free(bdenf);
			
			/*
			real *nGau = (real*) malloc(dom[dev].Gfz.s2b * sizeof(real));
			for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
				for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
					ii = i - dom[dev].Gfz.isb;
					jj = j - dom[dev].Gfz.jsb;
					C = i + j * Dom.Gfz.s1b;
					CC = ii + jj * dom[dev].Gfz.s1b;
					nGau[CC] = numGau[C];
				}
			}
			(cudaMemcpy(_numGau[dev], nGau, sizeof(real) * dom[dev].Gfz.s2b, cudaMemcpyHostToDevice));
			free(nGau);
			
			real *mGau = (real*) malloc(dom[dev].Gfz.s2b * sizeof(real));
			for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
				for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
					ii = i - dom[dev].Gfz.isb;
					jj = j - dom[dev].Gfz.jsb;
					C = i + j * Dom.Gfz.s1b;
					CC = ii + jj * dom[dev].Gfz.s1b;
					mGau[CC] = masGau[C];
				}
			}
			(cudaMemcpy(_masGau[dev], mGau, sizeof(real) * dom[dev].Gfz.s2b, cudaMemcpyHostToDevice));
			free(mGau);
			*/
		}
		
		
		
		// copy from host to device
		(cudaMemcpy(_numden[dev], nn, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_bubmas[dev], bubm, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_concen[dev], cc, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_bubden[dev], bden, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_concen_sat[dev], ccsat, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_BGndot[dev], bubgen_num, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_BGmdot[dev], bubgen_mas, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		
		// free host subdomain working arrays
		free(nn);
		free(bubm);
		free(cc);
		free(bden);
		free(ccsat);
		free(bubgen_num);
		free(bubgen_mas);
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
		(cudaSetDevice(dev + dev_start));

		// host working arrays for subdomain copy from device to host
		real *nn = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *cc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *ter = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *bubm = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		real *bdia = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
		
		// copy from device to host
		(cudaMemcpy(nn, _numden[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		(cudaMemcpy(cc, _concen[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		(cudaMemcpy(bdia, _bubdia[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		(cudaMemcpy(bubm, _bubmas[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		(cudaMemcpy(ter, _ter_cell[dev], sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost));
		
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
		// w_b
		for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
			for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
				for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
					ii = i - dom[dev].Gfz.isb;
					jj = j - dom[dev].Gfz.jsb;
					kk = k - dom[dev].Gfz.ksb;
					C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
					CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
					ter_cell[C] = ter[CC];
				}
			}
		}
		// bubmas
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubmas[C] = bubm[CC];
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
		for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
			for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
				for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
					ii = i - dom[dev].Gcc.isb;
					jj = j - dom[dev].Gcc.jsb;
					kk = k - dom[dev].Gcc.ksb;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
					bubdia[C] = bdia[CC];
				}
			}
		}
		
		// free host subdomain working arrays
		free(nn);
		free(cc);
		free(ter);
		free(bubm);
		free(bdia);
	}
}

extern "C"
void cuda_Eulerian_free(void)
{
	// free device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		(cudaFree(_numden[dev]));
		(cudaFree(_nextnumden[dev]));
		(cudaFree(_w_b[dev]));
		(cudaFree(_ter_cell[dev]));
		(cudaFree(_f_z_coupling_numden[dev]));
		(cudaFree(_BGndot[dev]));
		(cudaFree(_BGmdot[dev]));
		
		(cudaFree(_bubmas[dev]));
		(cudaFree(_nextbubmas[dev]));
		(cudaFree(_bubdia[dev]));
		(cudaFree(_bubdiafz[dev]));
		
		(cudaFree(_concen[dev]));
		(cudaFree(_nextconcen[dev]));
		(cudaFree(_velmag[dev]));
		(cudaFree(_mdot[dev]));
		(cudaFree(_bubden[dev]));
		(cudaFree(_bubden_face[dev]));
		(cudaFree(_concen_sat[dev]));
	}
	
	// free device memory on host
	free(_numden);
	free(_nextnumden);
	free(_w_b);
	free(_ter_cell);
	free(_f_z_coupling_numden);
	free(_BGndot);
	free(_BGmdot);
	
	free(_bubmas);
	free(_nextbubmas);
	free(_bubdia);
	free(_bubdiafz);
	
	free(_concen);
	free(_nextconcen);
	free(_velmag);
	free(_mdot);
	free(_bubden);
	free(_bubden_face);
	free(_concen_sat);
}

extern "C"
void cuda_numberdensity_BC(void)
{
	// CPU threading for multi-GPU
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
					BC_p_B_ghost_D<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
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
					BC_p_T_ghost_D<<<numBlocks_n, dimBlocks_n>>>(_numden[dev], _dom[dev]);
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
		(cudaSetDevice(dev + dev_start));
		
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
		                                                         _w_b[dev],
		                                                         _BGndot[dev]);
		
		kernel_inner_scalarfield_update_x<<<numBlocks_n, dimBlocks_n>>>(_dom[dev],
		                                                                _numden[dev],
		                                                                _nextnumden[dev]);
	}
}


extern "C"
void cuda_add_terminal_velz(void)
{
	// parallelize over CPU threads
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int blocks_x = 0;
		int threads_y = 0;
		int blocks_y = 0;
		int threads_z = 0;
		int blocks_z = 0;
		
		//======================================================================
		// generate cell-centered bubdia using bubmas
		
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
		
		dim3 dimBlocks_Gcc_x(threads_y, threads_z);
		dim3 numBlocks_Gcc_x(blocks_y, blocks_z);
		
		kernel_compute_bubble_diameter<<<numBlocks_Gcc_x, dimBlocks_Gcc_x>>>(_dom[dev],
		                                                                     _bubmas[dev],
		                                                                     _numden[dev],
		                                                                     _bubdia[dev],
		                                                                     _bubden[dev],
		                                                                     rho_f,
		                                                                     pressure_atm,
		                                                                     rho_atm,
		                                                                     grav_acc);
		
		//======================================================================
		// use cell-centered bubble diameter to calculate the terminal velocity
		// at cell centers
		
		// since bubble diameter is a field, here combine all the values
		// in terminal velocity except bubble diameter
		real cons = grav_acc / 18.0 / mu;
		
		kernel_compute_terminal_velocity<<<numBlocks_Gcc_x, dimBlocks_Gcc_x>>>(_dom[dev],
		                                                                       _ter_cell[dev],
		                                                                       _bubdia[dev],
		                                                                       _bubden[dev],
		                                                                       rho_f,
		                                                                       cons,
		                                                                       TerminalVelocityLimit);
		
		//======================================================================
		// Add the terminal velocity to fluid velocity field to get bubble
		// velocity
		
		if(dom[dev].Gfz._in < MAX_THREADS_DIM)
			threads_x = dom[dev].Gfz._in;
		else
			threads_x = MAX_THREADS_DIM;
		if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
			threads_y = dom[dev].Gfz._jn;
		else
			threads_y = MAX_THREADS_DIM;
		
		blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
		blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);
		
		dim3 dimBlocks_Gfz_z(threads_x, threads_y);
		dim3 numBlocks_Gfz_z(blocks_x, blocks_y);
		/*
		real rho_bot = rho_atm * (1.0 + rho_f * grav_acc * Dom.zl / pressure_atm);
		real mass_init = bubmasBC.nBD / numdenBC.nBD;
		real v_init = mass_init / rho_bot;
		real d_init = pow(6.0 * v_init / PI, 1./3.);
		real u_ter_init = d_init * d_init * (rho_f - rho_bot) * grav_acc / 18. / mu;
		*/
		kernel_numberdensity_particle_velz<<<numBlocks_Gfz_z, dimBlocks_Gfz_z>>>(_dom[dev],
		                                                                         _w_b[dev],
		                                                                         _w[dev],
		                                                                         _ter_cell[dev]);
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
	_w_b = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_ter_cell = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_f_z_coupling_numden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_BGndot = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_BGmdot = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	_bubmas = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_nextbubmas = (real**) malloc(nsubdom * sizeof(real*));
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
	_bubden = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_bubden_face = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	_concen_sat = (real**) malloc(nsubdom * sizeof(real*));
	cpumem += nsubdom * sizeof(real*);
	
	// allocate device memory on device
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		// number density equation
		(cudaMalloc((void**) &(_numden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _numden_next is to store _numden in next time step, only exist on deivce.
		(cudaMalloc((void**) &(_nextnumden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		//(cudaMalloc((void**) &(_u_p[dev]), sizeof(real) * dom[dev].Gfx.s3b));
		//gpumem += dom[dev].Gfx.s3b * sizeof(real);
		//(cudaMalloc((void**) &(_v_p[dev]), sizeof(real) * dom[dev].Gfy.s3b));
		//gpumem += dom[dev].Gfy.s3b * sizeof(real);
		(cudaMalloc((void**) &(_w_b[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		(cudaMalloc((void**) &(_ter_cell[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		//(cudaMalloc((void**) &(_f_x_coupling_numden[dev]), sizeof(real) * dom[dev].Gfx.s3b));
		//gpumem += dom[dev].Gfx.s3b * sizeof(real);
		//(cudaMalloc((void**) &(_f_y_coupling_numden[dev]), sizeof(real) * dom[dev].Gfy.s3b));
		//gpumem += dom[dev].Gfy.s3b * sizeof(real);
		(cudaMalloc((void**) &(_f_z_coupling_numden[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		/*
		(cudaMalloc((void**) &(_numGau[dev]), sizeof(real) * dom[dev].Gfz.s2b));
		gpumem += dom[dev].Gfz.s2b * sizeof(real);
		(cudaMalloc((void**) &(_masGau[dev]), sizeof(real) * dom[dev].Gfz.s2b));
		gpumem += dom[dev].Gfz.s2b * sizeof(real);
		(cudaMalloc((void**) &(_numGau_temp[dev]), sizeof(real) * dom[dev].Gfz.s2b));
		gpumem += dom[dev].Gfz.s2b * sizeof(real);
		(cudaMalloc((void**) &(_masGau_temp[dev]), sizeof(real) * dom[dev].Gfz.s2b));
		gpumem += dom[dev].Gfz.s2b * sizeof(real);
		*/
		(cudaMalloc((void**) &(_BGndot[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_BGmdot[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		
		// bubble mass
		(cudaMalloc((void**) &(_bubmas[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_nextbubmas[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_bubdia[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_bubdiafz[dev]), sizeof(real) * dom[dev].Gfz.s3b));
		gpumem += dom[dev].Gfz.s3b * sizeof(real);
		
		// concentration
		(cudaMalloc((void**) &(_concen[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _nextconcen is to store _concen in next time step, only exist on deivce.
		(cudaMalloc((void**) &(_nextconcen[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		// _velmag is to store the velocity magnitude of the fluid flow,
		// which is used in the source term of mass transfer, only exist on device.
		(cudaMalloc((void**) &(_velmag[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_mdot[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_bubden[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		(cudaMalloc((void**) &(_concen_sat[dev]), sizeof(real) * dom[dev].Gcc.s3b));
		gpumem += dom[dev].Gcc.s3b * sizeof(real);
		
		if(gravity_direction == GRAVITY_X) {
			(cudaMalloc((void**) &(_bubden_face[dev]), sizeof(real) * dom[dev].Gfx.s3b));
			gpumem += dom[dev].Gfx.s3b * sizeof(real);
		} else if(gravity_direction == GRAVITY_Y) {
			(cudaMalloc((void**) &(_bubden_face[dev]), sizeof(real) * dom[dev].Gfy.s3b));
			gpumem += dom[dev].Gfy.s3b * sizeof(real);
		} else if(gravity_direction == GRAVITY_Z) {
			(cudaMalloc((void**) &(_bubden_face[dev]), sizeof(real) * dom[dev].Gfz.s3b));
			gpumem += dom[dev].Gfz.s3b * sizeof(real);
		}
	}
}

extern "C"
void cuda_concentration_BC(void)
{
	// CPU threading for multi-GPU
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
					//BC_p_B_D<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev], concenBC.nBD);
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
					//BC_p_T_D<<<numBlocks_n, dimBlocks_n>>>(_concen[dev], _dom[dev], concenBC.nTD);
					break;
			}
		}
		// *********************************************************************
	}
}

extern "C"
void cuda_compute_coupling_force(void)
{
	// parallelize over CPU threads
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));

		int threads_x = 0;
		int threads_y = 0;
		int threads_z = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		int blocks_z = 0;

		// x-component
		if(dom[dev].Gfz._jnb < MAX_THREADS_DIM)
			threads_y = dom[dev].Gfz._jnb;
		else
			threads_y = MAX_THREADS_DIM;

		if(dom[dev].Gfz._knb < MAX_THREADS_DIM)
			threads_z = dom[dev].Gfz._knb;
		else
			threads_z = MAX_THREADS_DIM;

		blocks_y = (int)ceil((real) dom[dev].Gfz._jnb / (real) threads_y);
		blocks_z = (int)ceil((real) dom[dev].Gfz._knb / (real) threads_z);

		dim3 dimBlocks_x(threads_y, threads_z);
		dim3 numBlocks_x(blocks_y, blocks_z);

		// y-component
		if(dom[dev].Gfz._knb < MAX_THREADS_DIM)
			threads_z = dom[dev].Gfz._knb;
		else
			threads_z = MAX_THREADS_DIM;

		if(dom[dev].Gfz._inb < MAX_THREADS_DIM)
			threads_x = dom[dev].Gfz._inb;
		else
			threads_x = MAX_THREADS_DIM;

		blocks_z = (int)ceil((real) dom[dev].Gfz._knb / (real) threads_z);
		blocks_x = (int)ceil((real) dom[dev].Gfz._inb / (real) threads_x);

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
		
		// reset forcing arrays
		forcing_reset_x<<<numBlocks_x, dimBlocks_x>>>(_f_x[dev], _dom[dev]);
		forcing_reset_y<<<numBlocks_y, dimBlocks_y>>>(_f_y[dev], _dom[dev]);
		forcing_reset_z<<<numBlocks_z, dimBlocks_z>>>(_f_z[dev], _dom[dev]);
		
		// now add in the forcing
		/*
		if(quiescent_fluid == OFF) {
			real forcing_scale = 1.0 / 6.0 * PI * grav_acc / rho_f;
			kernel_forcing_add_z_field_bubble<<<numBlocks_z, dimBlocks_z>>>(forcing_scale,
		                                                                    _f_z_coupling_numden[dev],
		                                                                    _bubdiafz[dev],
		                                                                    _f_z[dev],
		                                                                    _dom[dev],
		                                                                    _bubden_face[dev],
		                                                                    rho_f);
		}
		*/
		// change the numerical implementation of forceing
		if(quiescent_fluid == OFF) {
			real forcing_scale = 1.0 / 6.0 * PI * grav_acc / rho_f;
			kernel_forcing_add_z_field_bubble<<<numBlocks_z, dimBlocks_z>>>(forcing_scale,
		                                                                    _numden[dev],
		                                                                    _bubdia[dev],
		                                                                    _f_z[dev],
		                                                                    _dom[dev],
		                                                                    _bubden[dev],
		                                                                    rho_f);
		}
	}
}

extern "C"
void cuda_concentration_march(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
		(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int threads_y = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		
		// N is the number of inner cells
		int N = dom[dev].Gcc.s3;
		
		// create temporary storage for reduction algorithms
		real *_numden_tmp;
		(cudaMalloc((void**) &_numden_tmp, sizeof(real) * N));
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
		totalnumden = totalnumden * Dom.dx * Dom.dy * Dom.dz;
		printf("Total bubble number(on device) = %f\n\n", totalnumden);
		
		// clean up
		(cudaFree(_numden_tmp));
	}
}

extern "C"
void cuda_bubblemass_BC(void)
{
	// CPU threading for multi-GPU
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
			switch(bubmasBC.nW) {
				case PERIODIC:
					BC_p_W_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
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
			switch(bubmasBC.nE) {
				case PERIODIC:
					BC_p_E_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
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
			switch(bubmasBC.nS) {
				case PERIODIC:
					BC_p_S_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
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
			switch(bubmasBC.nN) {
					case PERIODIC:
					BC_p_N_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
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
			switch(bubmasBC.nB) {
				case PERIODIC:
					BC_p_B_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_B_ghost_D<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_B_N<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
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
			switch(bubmasBC.nT) {
				case PERIODIC:
					BC_p_T_P<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
					break;
				case DIRICHLET:
					BC_p_T_ghost_D<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
					break;
				case NEUMANN:
					BC_p_T_N<<<numBlocks, dimBlocks>>>(_bubmas[dev], _dom[dev]);
					break;
			}
		}
		// *********************************************************************
		
		
	}
}

extern "C"
void cuda_coupling_force_preparation(void)
{
	// parallelize over CPU threads
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int threads_y = 0;
		//int threads_z = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		//int blocks_z = 0;
		
		/*//======================================================================
		// Generate bubdia using bubmas(cell-centered field)
		
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
		
		dim3 dimBlocks_Gcc_x(threads_y, threads_z);
		dim3 numBlocks_Gcc_x(blocks_y, blocks_z);
		
		kernel_compute_bubble_diameter<<<numBlocks_Gcc_x, dimBlocks_Gcc_x>>>(_dom[dev],
		                                                                     _bubmas[dev],
		                                                                     _numden[dev],
		                                                                     _bubdia[dev],
		                                                                     _bubden[dev],
		                                                                     rho_f,
		                                                                     pressure_atm,
		                                                                     rho_atm,
		                                                                     grav_acc);
		
		*///======================================================================
		// Calculate the number density field on cell faces.
		
		// set up kernel call: Gfz-z-ghost
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
		
		kernel_fz_coupling_numden_generate<<<numBlocks_z, dimBlocks_z>>>(_f_z_coupling_numden[dev],
		                                                                 _numden[dev],
		                                                                 _w_b[dev],
		                                                                 _dom[dev]);
		
		//======================================================================
		// Calculate z-face-centered bubble diameter field, which is
		// needed in coupling force of momentum equation.
		
		// set up kernel call: z
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
		
		dim3 dimBlocks_Gfz_z(threads_x, threads_y);
		dim3 numBlocks_Gfz_z(blocks_x, blocks_y);
		
		//real bottom_vol = bubmasBC.nBD / numdenBC.nBD / rho_atm / (1.0 + rho_f * grav_acc * dom[dev].zl / pressure_atm);
		//real diafz_init = pow(6.0 * bottom_vol / PI, 1.0/3.0);
		
		kernel_compute_bubble_diameterfz<<<numBlocks_Gfz_z, dimBlocks_Gfz_z>>>(_dom[dev],
		                                                                       _bubdiafz[dev],
		                                                                       _bubdia[dev],
		                                                                       _w_b[dev]);
	}
}

extern "C"
void cuda_compute_mdot(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
                                                      _numden[dev],
                                                      _bubmas[dev],
                                                      _concen[dev],
                                                      _bubdia[dev],
                                                      _ter_cell[dev],
                                                      _mdot[dev],
                                                      _concen_sat[dev],
                                                      concen_diff,
                                                      nu);
	}
}

extern "C"
void cuda_bubblemass_march(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
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
		
		// march the bubble mass equation
		kernel_march_bubblemass<<<numBlocks, dimBlocks>>>(dt,
                                                          _dom[dev],
                                                          _bubmas[dev],
                                                          _nextbubmas[dev],
                                                          _u[dev],
                                                          _v[dev],
                                                          _w_b[dev],
                                                          _numden[dev],
                                                          _mdot[dev],
                                                          _BGmdot[dev]);
		
		// update bubmas field
		kernel_inner_scalarfield_update_x<<<numBlocks, dimBlocks>>>(_dom[dev],
		                                                            _bubmas[dev],
		                                                            _nextbubmas[dev]);
	}
}

extern "C"
real cuda_Eulerian_find_dt(void)
{
  // results from all devices
  real *dts0 = (real*) malloc(nsubdom * sizeof(real));
    // cpumem += nsubdom * sizeof(real);
  real *dts1 = (real*) malloc(nsubdom * sizeof(real));

  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    // search
    real u_max = find_max_mag(dom[dev].Gfx.s3b, _u[dev]);
    real v_max = find_max_mag(dom[dev].Gfy.s3b, _v[dev]);
    real w_max = find_max_mag(dom[dev].Gfz.s3b, _w[dev]);
    real wp_max = find_max_mag(dom[dev].Gfz.s3b, _w_b[dev]);
    
    dts0[dev] = (u_max + 2. * nu / dom[dev].dx) / dom[dev].dx;
    dts0[dev] += (v_max + 2. * nu / dom[dev].dy) / dom[dev].dy;
    dts0[dev] += (w_max + 2. * nu / dom[dev].dz) / dom[dev].dz;
    //dts0[dev] += u_max * u_max / (2.0 * nu);
    //dts0[dev] += v_max * v_max / (2.0 * nu);
    //dts0[dev] += w_max * w_max / (2.0 * nu);
    dts0[dev] = CFL / dts0[dev];
    
    dts1[dev] = (u_max + 2. * concen_diff / dom[dev].dx) / dom[dev].dx;
    dts1[dev] += (v_max + 2. * concen_diff / dom[dev].dy) / dom[dev].dy;
    dts1[dev] += (wp_max + 2. * concen_diff / dom[dev].dz) / dom[dev].dz;
    
    //dts1[dev] += u_max * u_max / (2.0 * concen_diff);
    //dts1[dev] += v_max * v_max / (2.0 * concen_diff);
    //dts1[dev] += wp_max * wp_max / (2.0 * concen_diff);
    dts1[dev] = CFL / dts1[dev];
    
    if(dts1[dev] < dts0[dev]) dts0[dev] = dts1[dev];
  }

  // find min of all devices
  real min = FLT_MAX;
  for(int i = 0; i < nsubdom; i++)
    if(dts0[i] < min) min = dts0[i];

  // clean up
  free(dts0);
  free(dts1);

  return min;
}

extern "C"
void cuda_Eulerian_compute_forcing(real *pid_int, real *pid_back, real Kp, real Ki,
  real Kd)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

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

    //##########################################################################
    // reset forcing arrays
    //forcing_reset_x<<<numBlocks_x, dimBlocks_x>>>(_f_x[dev], _dom[dev]);
    //forcing_reset_y<<<numBlocks_y, dimBlocks_y>>>(_f_y[dev], _dom[dev]);
    //forcing_reset_z<<<numBlocks_z, dimBlocks_z>>>(_f_z[dev], _dom[dev]);
    //##########################################################################

    // linearly accelerate pressure gradient from zero
    real delta = ttime - p_tDelay;
    if (delta >= 0 ) {
      if(gradP.xa == 0) gradP.x = gradP.xm;
      else if(fabs(delta*gradP.xa) > fabs(gradP.xm)) gradP.x = gradP.xm;
      else gradP.x = delta*gradP.xa;

      if(gradP.ya == 0) gradP.y = gradP.ym;
      else if(fabs(delta*gradP.ya) > fabs(gradP.ym)) gradP.y = gradP.ym;
      else gradP.y = delta*gradP.ya;

      // turn off if PID controller is being used
      if(!(Kp > 0 || Ki > 0 || Kd > 0)) {
        if(gradP.za == 0) gradP.z = gradP.zm;
        else if(fabs(delta*gradP.za) > fabs(gradP.zm)) gradP.z = gradP.zm;
        else gradP.z = delta*gradP.za;
      }
    }

    // linearly accelerate gravitational acceleration from zero
    delta = ttime - g_tDelay;
    if (delta >= 0) {
      if(g.xa == 0) g.x = g.xm;
      else if(fabs(delta*g.xa) > fabs(g.xm)) g.x = g.xm;
      else g.x = delta*g.xa;

      if(g.ya == 0) g.y = g.ym;
      else if(fabs(delta*g.ya) > fabs(g.ym)) g.y = g.ym;
      else g.y = delta*g.ya;

      if(g.za == 0) g.z = g.zm;
      else if(fabs(delta*g.za) > fabs(g.zm)) g.z = g.zm;
      else g.z = delta*g.za;
    }

    delta = ttime - p_tDelay;
    // PID controller  TODO: make this into a kernel function
    if (delta >= 0) {
      if(Kp > 0 || Ki > 0 || Kd > 0) {
        cuda_part_pull();
        real acc_z = 0;
        real volp = 0;
        real massp = 0;
        for(int i = 0; i < nparts; i++) {
          volp += parts[i].r*parts[i].r*parts[i].r;
          massp += parts[i].rho*parts[i].r*parts[i].r*parts[i].r;
          acc_z += parts[i].wdot;
        }
        volp *= 4./3.*PI;
        massp *= 4./3.*PI;
        real volfrac = volp / (Dom.xl * Dom.yl * Dom.zl);
        real rho_avg = massp/volp*volfrac + rho_f*(1.-volfrac);
        acc_z /= nparts;

        *pid_int = *pid_int + acc_z*dt;
        gradP.z = gradP.z
          + (Kp*acc_z + Ki*(*pid_int)/ttime + (Kd)*(acc_z-*pid_back))*rho_avg;
        *pid_back = acc_z;
      }
    }
    forcing_add_x_const<<<numBlocks_x, dimBlocks_x>>>(-gradP.x / rho_f,
      _f_x[dev], _dom[dev]);
    forcing_add_y_const<<<numBlocks_y, dimBlocks_y>>>(-gradP.y / rho_f,
      _f_y[dev], _dom[dev]);
    forcing_add_z_const<<<numBlocks_z, dimBlocks_z>>>(-gradP.z / rho_f,
      _f_z[dev], _dom[dev]);

  }
}

extern "C"
void cuda_botbubbvel_BC(void)
{
	// CPU threading for multi-GPU
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		int threads_x = 0;
		int threads_y = 0;
		int blocks_x = 0;
		int blocks_y = 0;
		
		// check whether bottom boundary is an external boundary
		if(dom[dev].B == -1) {
			// set up kernel call
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
			
			dim3 dimBlocks_w(threads_x, threads_y);
			dim3 numBlocks_w(blocks_x, blocks_y);
			
			// apply zero-Dirichlet BC to bubble velocity for this face
			BC_w_B_D<<<numBlocks_w, dimBlocks_w>>>(_w_b[dev], _dom[dev], 0.);
		}
	}
}

extern "C"
void cuda_compute_HYPERBOLICTANRANDOM(void)
{
	#pragma omp parallel num_threads(nsubdom)
	{
		int dev = omp_get_thread_num();
		(cudaSetDevice(dev + dev_start));
		
		const int bubble_number_this_step = lrint(dt * BubbleGenerator.BubGen_bubblenumber);
		int C, i, j, k;
		int *keys;
		int *keys_out;
		int *values;
		int *values_out;
		
		// allocate memory on host or device
		bub_gen_pos = (bubble_struct*) malloc(bubble_number_this_step * sizeof(bubble_struct));
		keys     = (int*) malloc(bubble_number_this_step * sizeof(int));
		keys_out = (int*) malloc(bubble_number_this_step * sizeof(int));
		values     = (int*) malloc(bubble_number_this_step * sizeof(int));
		values_out = (int*) malloc(bubble_number_this_step * sizeof(int));
		
		//_bub_gen_pos = (bubble_struct**) malloc(nsubdom * sizeof(bubble_struct*));
		for(i = 0; i < bubble_number_this_step; i++) {
			keys[i] = 0;
			keys_out[i] = 0;
			values[i] = 1;
			values_out[i] = 1;
		}
		
		for(C = 0; C < bubble_number_this_step; C++) {
			
			real temp_x = (real)rand()/(real)(RAND_MAX);
			real temp_y = (real)rand()/(real)(RAND_MAX);
			real temp_z = (real)rand()/(real)(RAND_MAX);
			
			for(i = 0; i < Dom.xn && temp_x > Ix[i]; i++);
			if(i > 0) {
				float x1 = Dom.xs + Dom.dx * (i - 0.5);
				float x2 = x1 + Dom.dx;
				bub_gen_pos[C].x = (temp_x - Ix[i-1]) / (Ix[i] - Ix[i-1]) * x2 + (Ix[i] - temp_x) / (Ix[i] - Ix[i-1]) * x1;
			} else {
				float x1 = Dom.xs;
				float x2 = Dom.xs + Dom.dx * 0.5;
				bub_gen_pos[C].x = (temp_x - 0.) / Ix[0] * x2 + (Ix[0] - temp_x) / Ix[0] * x1;
			}
			
			for(j = 0; j < Dom.yn && temp_y > Iy[j]; j++);
			if(j > 0) {
				float y1 = Dom.ys + Dom.dy * (j - 0.5);
				float y2 = y1 + Dom.dy;
				bub_gen_pos[C].y = (temp_y - Iy[j-1]) / (Iy[j] - Iy[j-1]) * y2 + (Iy[j] - temp_y) / (Iy[j] - Iy[j-1]) * y1;
			} else {
				float y1 = Dom.ys;
				float y2 = Dom.ys + Dom.dy * 0.5;
				bub_gen_pos[C].y = (temp_y - 0.) / Iy[0] * y2 + (Iy[0] - temp_y) / Iy[0] * y1;
			}
			
			for(k = 0; k < Dom.zn && temp_z > Iz[k]; k++);
			if(k > 0) {
				float z1 = Dom.zs + Dom.dz * (k - 0.5);
				float z2 = z1 + Dom.dz;
				bub_gen_pos[C].z = (temp_z - Iz[k-1]) / (Iz[k] - Iz[k-1]) * z2 + (Iz[k] - temp_z) / (Iz[k] - Iz[k-1]) * z1;
			} else {
				float z1 = Dom.zs;
				float z2 = Dom.zs + Dom.dz * 0.5;
				bub_gen_pos[C].z = (temp_z - 0.) / Iz[0] * z2 + (Iz[0] - temp_z) / Iz[0] * z1;
			}
		}
		
		// calculate keys
		for(C = 0; C < bubble_number_this_step; C++) {
			i = (int)ceil(((bub_gen_pos[C].x - Dom.xs) / Dom.dx));
			j = (int)ceil(((bub_gen_pos[C].y - Dom.ys) / Dom.dy));
			k = (int)ceil(((bub_gen_pos[C].z - Dom.zs) / Dom.dz));
			keys[C] = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
			//printf("number = %d, i = %d, j = %d, k = %d, key = %d, value = %d\n", bubble_number_this_step, i, j, k, keys[C], values[C]);
		}
		
		// sort
		thrust::sort(keys, keys + bubble_number_this_step);
		
		// reduce by key
		thrust::reduce_by_key(keys, keys + bubble_number_this_step, values, keys_out, values_out);
		
		// reset bubble injection rate
		for(C = 0; C < Dom.Gcc.s3b; C++) {
			BGndot[C] = 0.0;
			BGmdot[C] = 0.0;
		}
		
		// calculate
		for(C = 0; C < bubble_number_this_step; C++) {
			if(keys_out[C] > 0) {
				BGndot[keys_out[C]] = values_out[C] / (Dom.dx * Dom.dy * Dom.dz) / dt;
				BGmdot[keys_out[C]] = PI / 6. * pow(BubbleGenerator.BubGen_dia, 3.0) * bubden[keys_out[C]] * BGndot[keys_out[C]];
				//printf("position: %d, ndot: %f, mdot: %f\n", keys_out[C], BGndot[keys_out[C]], BGmdot[keys_out[C]]);
			}
		}
		
		// push memory
		(cudaMemcpy(_BGndot[dev], BGndot, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		(cudaMemcpy(_BGmdot[dev], BGmdot, sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyHostToDevice));
		
		// free memory
		free(keys);
		free(values);
		free(keys_out);
		free(values_out);
	}
}
