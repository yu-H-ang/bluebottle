#include "Eulerian.h"

void Eulerian_read_input(void)
{
	int fret = 0;
	fret = fret; // prevent compiler warning
	
	// open configuration file for reading
	char fname[FILE_NAME_SIZE];
	sprintf(fname, "%s/input/Eulerian.config", ROOT_DIR);
	FILE *infile = fopen(fname, "r");
	if(infile == NULL) {
		fprintf(stderr, "Could not open file %s\n", fname);
		exit(EXIT_FAILURE);
	}
	
	char buf[CHAR_BUF_SIZE];  // character read buffer
	
	//==========================================================================
	
	fret = fscanf(infile, "BUBBLE PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "bubble_radius %lf\n", &bubble_radius);
		fret = fscanf(infile, "bubble_density %lf\n", &bubble_density);
	#else // single
		fret = fscanf(infile, "bubble_radius %f\n", &bubble_radius);
		fret = fscanf(infile, "bubble_density %f\n", &bubble_density);
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "DISSOLUTION PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "concen_diff %lf\n", &concen_diff);
		fret = fscanf(infile, "concen_diss %lf\n", &concen_diss);
	#else // single
		fret = fscanf(infile, "concen_diff %f\n", &concen_diff);
		fret = fscanf(infile, "concen_diss %f\n", &concen_diss);
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "NUMBER DENSITY BOUNDARY CONDITIONS\n");
	
	fret = fscanf(infile, "numdenBC.nW %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nW = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "numdenBC.nE %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nE = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "numdenBC.nS %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nS = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "numdenBC.nN %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nN = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	
	fret = fscanf(infile, "numdenBC.nB %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nB = PERIODIC;
	} else if(strcmp(buf, "NEUMANN") == 0) {
		numdenBC.nB = NEUMANN;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		numdenBC.nB = DIRICHLET;
		fret = fscanf(infile, " %lf", &numdenBC.nBD);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	
	fret = fscanf(infile, "numdenBC.nT %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		numdenBC.nT = PERIODIC;
	} else if(strcmp(buf, "NEUMANN") == 0) {
		numdenBC.nT = NEUMANN;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		numdenBC.nT = DIRICHLET;
		fret = fscanf(infile, " %lf", &numdenBC.nTD);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "BUBBLE INITIAL CONDITION\n");
	fret = fscanf(infile, "bubble_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		bubble_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &bubble_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
		bubble_init_cond = RANDOM;
		fret = fscanf(infile, " %d", &bubble_init_cond_random_min);
		fret = fscanf(infile, " %d", &bubble_init_cond_random_max);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		bubble_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &bubble_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
		fret = fscanf(infile, " %d", &bubble_init_cond_random_min);
		fret = fscanf(infile, " %d", &bubble_init_cond_random_max);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "BUBBLE VOLUME BOUNDARY CONDITIONS\n");
	
	fret = fscanf(infile, "bubvolBC.nW %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nW = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubvolBC.nE %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nE = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubvolBC.nS %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nS = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubvolBC.nN %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nN = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubvolBC.nB %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nB = PERIODIC;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		bubvolBC.nB = DIRICHLET;
		fret = fscanf(infile, " %lf", &bubvolBC.nBD);
	} else if(strcmp(buf, "NEUMANN") == 0) {
		bubvolBC.nB = NEUMANN;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	fret = fscanf(infile, "bubvolBC.nT %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubvolBC.nT = PERIODIC;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		bubvolBC.nT = DIRICHLET;
		fret = fscanf(infile, " %lf", &bubvolBC.nTD);
	} else if(strcmp(buf, "NEUMANN") == 0) {
		bubvolBC.nT = NEUMANN;
	} else {
	fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "BUBBLE VOLUME INITIAL CONDITION\n");
	fret = fscanf(infile, "bubvol_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		bubvol_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &bubvol_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		bubvol_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &bubvol_init_cond_uniform_m;);
	} else if(strcmp(buf, "RANDOM") == 0) {
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "CONCENTRATION BOUNDARY CONDITIONS\n");
	fret = fscanf(infile, "concenBC.nW %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nW = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "concenBC.nE %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nE = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "concenBC.nS %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nS = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "concenBC.nN %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nN = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "concenBC.nB %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nB = PERIODIC;
	} else if(strcmp(buf, "NEUMANN") == 0) {
		concenBC.nB = NEUMANN;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		concenBC.nB = DIRICHLET;
		fret = fscanf(infile, " %lf", &concenBC.nBD);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "concenBC.nT %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		concenBC.nT = PERIODIC;
	} else if(strcmp(buf, "NEUMANN") == 0) {
		concenBC.nT = NEUMANN;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		concenBC.nT = DIRICHLET;
		fret = fscanf(infile, " %lf", &concenBC.nTD);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "CONCENTRATION INITIAL CONDITION\n");
	fret = fscanf(infile, "concen_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		concen_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &concen_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
		concen_init_cond = RANDOM;
		fret = fscanf(infile, " %d", &concen_init_cond_random_min);
		fret = fscanf(infile, " %d", &concen_init_cond_random_max);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		concen_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &concen_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
		concen_init_cond = RANDOM;
		fret = fscanf(infile, " %d", &concen_init_cond_random_min);
		fret = fscanf(infile, " %d", &concen_init_cond_random_max);
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#endif
	
	
	fclose(infile);
}

void Eulerian_clean(void)
{
	free(numden);
	free(bubvol);
	free(concen);
	//free(u_p);
	//free(v_p);
	free(w_p);
 }

int Eulerian_init(void)
{
	int i, j, k, C;  // iterators

	// make sure there are enough GPU devices in the given range
	if(nsubdom > dev_end - dev_start + 1) {
		return EXIT_FAILURE;
	}
	
	// allocate memory
	numden = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	//u_p = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
	//cpumem += Dom.Gfx.s3b * sizeof(real);
	//v_p = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
	//cpumem += Dom.Gfy.s3b * sizeof(real);
	w_p = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
	cpumem += Dom.Gfz.s3b * sizeof(real);
	concen = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	bubvol = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	
	// initialize UNIFORM number density field
	if(bubble_init_cond == UNIFORM) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					numden[C] = bubble_init_cond_uniform_m * k;
				}
			}
		}
	}
	// initialize RANDOM number density field
	if(bubble_init_cond == RANDOM) {
		srand(time(NULL));
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					numden[C] = (rand()%(bubble_init_cond_random_max - bubble_init_cond_random_min)) + bubble_init_cond_random_min;
				}
			}
		}
	}
	// initialize UNIFORM bubble volume field
	if(bubvol_init_cond == UNIFORM) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					bubvol[C] = bubvol_init_cond_uniform_m;
				}
			}
		}
	}
	// initialize UNIFORM concentration field
	if(concen_init_cond == UNIFORM) {
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					concen[C] = concen_init_cond_uniform_m;
				}
			}
		}
	}
	// initialize RANDOM concentration field
	if(concen_init_cond == RANDOM) {
		srand(time(NULL));
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					concen[C] = (rand()%(concen_init_cond_random_max - concen_init_cond_random_min)) + concen_init_cond_random_min;
				}
			}
		}
	}
	
	// initialize particle velocity field to be quiescent;
	/*
	for(i = 0; i < Dom.Gfx.s3b; i++) {
		u_p[i] = 0.;
	}
	for(j = 0; j < Dom.Gfy.s3b; j++) {
		v_p[i] = 0.;
	}
	*/
	for(k = 0; k < Dom.Gfz.s3b; k++) {
		w_p[k] = 0.;
	}
	
	// initialize some values
	// totalnumden is the sum of number density inside the domain(without ghost cell)
	totalnumden = 0;
	for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
		for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
				C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
				totalnumden += numden[C];
			}
		}
	}
	
	return EXIT_SUCCESS;
}

void Eulerian_show_config()
{
	printf("########################################################################\n");
	printf("BUBBLE PARAMETERS\n");
	printf("bubble_radius ");
	printf("%f\n", bubble_radius);
	printf("bubble_density ");
	printf("%f\n", bubble_density);
	printf("\n");
	
	printf("DISSOLUTION PARAMETERS\n");
	printf("concen_diff ");
	printf("%f\n", concen_diff);
	printf("concen_diss ");
	printf("%f\n", concen_diss);
	printf("\n");
	
	printf("NUMBER DENSITY BOUNDARY CONDITIONS\n");
	printf("numdenBC.nW ");
	if(numdenBC.nW == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("numdenBC.nE ");
	if(numdenBC.nE == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("numdenBC.nS ");
	if(numdenBC.nS == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("numdenBC.nN ");
	if(numdenBC.nN == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("numdenBC.nB ");
	if(numdenBC.nB == PERIODIC) {
		printf("PERIODIC\n");
	} else if(numdenBC.nB == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", numdenBC.nBD);
	} else if(numdenBC.nB == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("numdenBC.nT ");
	if(numdenBC.nT == PERIODIC) {
		printf("PERIODIC\n");
	} else if(numdenBC.nT == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", numdenBC.nTD);
	} else if(numdenBC.nT == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("\n");
	
	printf("BUBBLE INITIAL CONDITION\n");
	printf("bubble_init_cond ");
	if(bubble_init_cond == RANDOM) {
		printf("RANDOM ");
		printf("%d ", bubble_init_cond_random_min);
		printf("%d\n", bubble_init_cond_random_max);
	} else if(bubble_init_cond == UNIFORM) {
		printf("UNIFORM ");
		printf("%f\n", bubble_init_cond_uniform_m);
	}
	printf("\n");
	
	printf("BUBBLE VOLUME BOUNDARY CONDITIONS\n");
	printf("bubvolBC.nW ");
	if(bubvolBC.nW == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubvolBC.nE ");
	if(bubvolBC.nE == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubvolBC.nS ");
	if(bubvolBC.nS == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubvolBC.nN ");
	if(bubvolBC.nN == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubvolBC.nB ");
	if(bubvolBC.nB == PERIODIC) {
		printf("PERIODIC\n");
	} else if(bubvolBC.nB == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", bubvolBC.nBD);
	} else if(bubvolBC.nB == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("bubvolBC.nT ");
	if(bubvolBC.nT == PERIODIC) {
		printf("PERIODIC\n");
	} else if(bubvolBC.nT == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", bubvolBC.nTD);
	} else if(bubvolBC.nT == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("\n");
	
	printf("BUBBLE VOLUME INITIAL CONDITION\n");
	printf("bubvol_init_cond ");
	if(bubvol_init_cond == UNIFORM) {
		printf("UNIFORM ");
		printf("%f\n", bubvol_init_cond_uniform_m);
	}
	printf("\n");
	
	printf("CONCENTRATION BOUNDARY CONDITIONS\n");
	printf("concenBC.nW ");
	if(concenBC.nW == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("concenBC.nE ");
	if(concenBC.nE == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("concenBC.nS ");
	if(concenBC.nS == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("concenBC.nN ");
	if(concenBC.nN == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("concenBC.nB ");
	if(concenBC.nB == PERIODIC) {
		printf("PERIODIC\n");
	} else if(concenBC.nB == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", concenBC.nBD);
	} else if(concenBC.nB == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("concenBC.nT ");
	if(concenBC.nT == PERIODIC) {
		printf("PERIODIC\n");
	} else if(concenBC.nT == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", concenBC.nTD);
	} else if(concenBC.nT == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("\n");
	
	printf("CONCENTRATION INITIAL CONDITION\n");
	printf("concen_init_cond ");
	if(concen_init_cond == UNIFORM) {
		printf("UNIFORM ");
		printf("%f\n", concen_init_cond_uniform_m);
	} else if(concen_init_cond == RANDOM) {
		printf("RANDOM ");
		printf("%d ", concen_init_cond_random_min);
		printf("%d\n", concen_init_cond_random_max);
	}
	printf("\n");
	
	printf("########################################################################\n");
}
