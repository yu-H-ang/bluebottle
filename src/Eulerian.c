#include "Eulerian.h"

void numberdensity_read_input(void)
{
	int fret = 0;
	fret = fret; // prevent compiler warning
	
	// open configuration file for reading
	char fname[FILE_NAME_SIZE];
	sprintf(fname, "%s/input/numden.config", ROOT_DIR);
	FILE *infile = fopen(fname, "r");
	if(infile == NULL) {
		fprintf(stderr, "Could not open file %s\n", fname);
		exit(EXIT_FAILURE);
	}
	
	char buf[CHAR_BUF_SIZE];  // character read buffer
	
	fret = fscanf(infile, "BUBBLE PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "bubble_radius %lf\n", &bubble_radius);
		fret = fscanf(infile, "bubble_density %lf\n", &bubble_density);
	#else // single
		fret = fscanf(infile, "bubble_radius %f\n", &bubble_radius);
		fret = fscanf(infile, "bubble_density %f\n", &bubble_density);
	#endif
	
	fret = fscanf(infile, "\n");//======================================
	
	fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
	fret = fscanf(infile, "NUMBER DENSITY\n");

		fret = fscanf(infile, "numdenBC.nW %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nW = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "numdenBC.nE %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nE = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "numdenBC.nS %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nS = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "numdenBC.nN %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nN = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "numdenBC.nB %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nB = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "numdenBC.nT %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			numdenBC.nT = PERIODIC;
		} else {
			fprintf(stderr, "numden.config read error.\n");
			exit(EXIT_FAILURE);
		}
	
	fret = fscanf(infile, "\n");//======================================
	
	fret = fscanf(infile, "BUBBLE INITIAL CONDITION\n");
	fret = fscanf(infile, "bubble_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		bubble_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &bubble_init_cond_uniform_m);
	} else {
		fprintf(stderr, "numden.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		bubble_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &bubble_init_cond_uniform_m);
	} else {
		fprintf(stderr, "numden.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#endif
	
	fclose(infile);
	
	/*
	printf("bubble_radius = %f\n" ,bubble_radius);
	printf("bubble_density = %f\n" ,bubble_density);
	printf("bubble_init_cond = %d\n", UNIFORM);
	printf("bubble_init_cond_uniform_m = %f\n", bubble_init_cond_uniform_m);
	printf("numdenBC.nW = %d\n", numdenBC.nW);
	printf("numdenBC.nE = %d\n", numdenBC.nE);
	printf("numdenBC.nS = %d\n", numdenBC.nS);
	printf("numdenBC.nN = %d\n", numdenBC.nN);
	printf("numdenBC.nB = %d\n", numdenBC.nB);
	printf("numdenBC.nT = %d\n", numdenBC.nT);
	*/
}

int numberdensity_init(void)
{
	int i, j, k, C;  // iterators

	// make sure there are enough GPU devices in the given range
	if(nsubdom > dev_end - dev_start + 1) {
		return EXIT_FAILURE;
	}
	
	// allocate and initialize number density vectors
	numden = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	u_p = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
	cpumem += Dom.Gfx.s3b * sizeof(real);
	v_p = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
	cpumem += Dom.Gfy.s3b * sizeof(real);
	w_p = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
	cpumem += Dom.Gfz.s3b * sizeof(real);
  
	// initialize UNIFORM number density field
	if(bubble_init_cond == UNIFORM) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					numden[C] = 10 * (k + 1);
				}
			}
		}
	}
	
	// initialize particle velocity field to be quiescent;
	for(i = 0; i < Dom.Gfx.s3b; i++) {
		u_p[i] = 0.;
	}
	for(j = 0; j < Dom.Gfy.s3b; j++) {
		v_p[i] = 0.;
	}
	for(k = 0; k < Dom.Gfz.s3b; k++) {
		w_p[i] = 0.;
	}
	
	// initialize RANDOM number density field //TODO
	if(bubble_density == 1) {
		
	}
	
	return EXIT_SUCCESS;
}

void numberdensity_clean(void)
{
	
	free(numden);
	free(u_p);
	free(v_p);
	free(w_p);
 }

void concentration_read_input(void)
{
	int fret = 0;
	fret = fret; // prevent compiler warning
	
	// open configuration file for reading
	char fname[FILE_NAME_SIZE];
	sprintf(fname, "%s/input/concen.config", ROOT_DIR);
	FILE *infile = fopen(fname, "r");
	if(infile == NULL) {
		fprintf(stderr, "Could not open file %s\n", fname);
		exit(EXIT_FAILURE);
	}
	
	char buf[CHAR_BUF_SIZE];  // character read buffer
	
	fret = fscanf(infile, "DISSOLUTION PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "concen_diff %lf\n", &concen_diff);
		fret = fscanf(infile, "concen_diss %lf\n", &concen_diss);
	#else // single
		fret = fscanf(infile, "concen_diff %f\n", &concen_diff);
		fret = fscanf(infile, "concen_diss %f\n", &concen_diss);
	#endif
	
	fret = fscanf(infile, "\n");//======================================
	
	fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
	fret = fscanf(infile, "CONCENTRATION\n");

		fret = fscanf(infile, "concenBC.nW %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nW = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "concenBC.nE %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nE = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "concenBC.nS %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nS = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "concenBC.nN %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nN = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "concenBC.nB %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nB = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
		
		fret = fscanf(infile, "concenBC.nT %s\n", buf);
		if(strcmp(buf, "PERIODIC") == 0) {
			concenBC.nT = PERIODIC;
		} else {
			fprintf(stderr, "concen.config read error.\n");
			exit(EXIT_FAILURE);
		}
	
	fret = fscanf(infile, "\n");//======================================
	
	fret = fscanf(infile, "CONCENTRATION INITIAL CONDITION\n");
	fret = fscanf(infile, "concen_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		concen_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &concen_init_cond_uniform_m);
	} else {
		fprintf(stderr, "concen.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		concen_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &concen_init_cond_uniform_m);
	} else {
		fprintf(stderr, "concen.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#endif
	
	fclose(infile);
	
	/*
	printf("concen_diff = %f\n" ,concen_diff);
	printf("concen_diss = %f\n" ,concen_diss);
	printf("concen_init_cond = %d\n", concen_init_cond);
	printf("concen_init_cond_uniform_m = %f\n", concen_init_cond_uniform_m);
	printf("concenBC.nW = %d\n", concenBC.nW);
	printf("concenBC.nE = %d\n", concenBC.nE);
	printf("concenBC.nS = %d\n", concenBC.nS);
	printf("concenBC.nN = %d\n", concenBC.nN);
	printf("concenBC.nB = %d\n", concenBC.nB);
	printf("concenBC.nT = %d\n", concenBC.nT);
	*/
}
