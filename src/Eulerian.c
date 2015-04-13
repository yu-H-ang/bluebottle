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
	
	fret = fscanf(infile, "ATMOSHPERE PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "rho_atm %lf\n", &rho_atm);
		fret = fscanf(infile, "pressure_atm %lf\n", &pressure_atm);
		fret = fscanf(infile, "grav_acc %lf\n", &grav_acc);
	#else // single
		fret = fscanf(infile, "rho_atm %f\n", &rho_atm);
		fret = fscanf(infile, "pressure_atm %f\n", &pressure_atm);
		fret = fscanf(infile, "grav_acc %f\n", &grav_acc);
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
	
	fret = fscanf(infile, "BUBBLE MASS BOUNDARY CONDITIONS\n");
	
	fret = fscanf(infile, "bubmasBC.nW %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nW = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubmasBC.nE %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nE = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubmasBC.nS %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nS = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubmasBC.nN %s\n", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nN = PERIODIC;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "bubmasBC.nB %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nB = PERIODIC;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		bubmasBC.nB = DIRICHLET;
		fret = fscanf(infile, " %lf", &bubmasBC.nBD);
	} else if(strcmp(buf, "NEUMANN") == 0) {
		bubmasBC.nB = NEUMANN;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	fret = fscanf(infile, "bubmasBC.nT %s", buf);
	if(strcmp(buf, "PERIODIC") == 0) {
		bubmasBC.nT = PERIODIC;
	} else if(strcmp(buf, "DIRICHLET") == 0) {
		bubmasBC.nT = DIRICHLET;
		fret = fscanf(infile, " %lf", &bubmasBC.nTD);
	} else if(strcmp(buf, "NEUMANN") == 0) {
		bubmasBC.nT = NEUMANN;
	} else {
	fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	fret = fscanf(infile, "\n");
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "BUBBLE MASS INITIAL CONDITION\n");
	fret = fscanf(infile, "bubmas_init_cond %s", buf);
	#ifdef DOUBLE
	if(strcmp(buf, "UNIFORM") == 0) {
		bubmas_init_cond = UNIFORM;
		fret = fscanf(infile, " %lf\n", &bubmas_init_cond_uniform_m);
	} else if(strcmp(buf, "RANDOM") == 0) {
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#else
	if(strcmp(buf, "UNIFORM") == 0) {
		bubmas_init_cond = UNIFORM;
		fret = fscanf(infile, " %f\n", &bubmas_init_cond_uniform_m;);
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
	free(bubmas);
	free(concen);
	free(w_b);
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
	w_b = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
	cpumem += Dom.Gfz.s3b * sizeof(real);
	concen = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	bubmas = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	
	// initialize UNIFORM number density field
	if(bubble_init_cond == UNIFORM) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					numden[C] = bubble_init_cond_uniform_m;
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
	// initialize UNIFORM bubble mass field
	if(bubmas_init_cond == UNIFORM) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					bubmas[C] = bubmas_init_cond_uniform_m;
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
	for(k = 0; k < Dom.Gfz.s3b; k++) {
		w_b[k] = 0.;
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
	printf("ATMOSHPERE PARAMETERS\n");
	printf("rho_atm ");
	printf("%f\n", rho_atm);
	printf("pressure_atm ");
	printf("%f\n", pressure_atm);
	printf("grav_acc ");
	printf("%f\n", grav_acc);
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
	
	printf("BUBBLE MASS BOUNDARY CONDITIONS\n");
	printf("bubmasBC.nW ");
	if(bubmasBC.nW == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubmasBC.nE ");
	if(bubmasBC.nE == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubmasBC.nS ");
	if(bubmasBC.nS == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubmasBC.nN ");
	if(bubmasBC.nN == PERIODIC) {
		printf("PERIODIC\n");
	}
	printf("bubmasBC.nB ");
	if(bubmasBC.nB == PERIODIC) {
		printf("PERIODIC\n");
	} else if(bubmasBC.nB == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", bubmasBC.nBD);
	} else if(bubmasBC.nB == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("bubmasBC.nT ");
	if(bubmasBC.nT == PERIODIC) {
		printf("PERIODIC\n");
	} else if(bubmasBC.nT == DIRICHLET) {
		printf("DIRICHLET ");
		printf("%f\n", bubmasBC.nTD);
	} else if(bubmasBC.nT == NEUMANN) {
		printf("NEUMANN\n");
	}
	printf("\n");
	
	printf("BUBBLE MASS INITIAL CONDITION\n");
	printf("bubmas_init_cond ");
	if(bubmas_init_cond == UNIFORM) {
		printf("UNIFORM ");
		printf("%f\n", bubmas_init_cond_uniform_m);
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
	printf("Some important information\n");
	//printf("%lf\n",Dom.zl);
	real rho_bot = rho_atm * (1.0 + rho_f * grav_acc * Dom.zl / pressure_atm);
	real mass_init = bubmasBC.nBD / numdenBC.nBD;
	real v_init = mass_init / rho_bot;
	real d_init = pow(6.0 * v_init / PI, 1./3.);
	real u_ter_init = d_init * d_init * (rho_f - rho_bot) * grav_acc / 18 / mu;
	printf("bubble density at bottom: %f\n", rho_bot);
	printf("terminal velocity: %f\n", u_ter_init);
	printf("########################################################################\n");
}
