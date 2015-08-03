#include "Eulerian.h"
#include "recorder.h"

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
	fret = fscanf(infile, "gravity_direction %s\n", buf);
	if(strcmp(buf, "GRAVITY_X") == 0) {
		gravity_direction = GRAVITY_X;
	} else if(strcmp(buf, "GRAVITY_Y") == 0) {
		gravity_direction = GRAVITY_Y;
	} else if(strcmp(buf, "GRAVITY_Z") == 0) {
		gravity_direction = GRAVITY_Z;
	}
	#ifdef DOUBLE
		fret = fscanf(infile, "TerminalVelocityLimit %lf\n", &TerminalVelocityLimit);
	#else // single
		fret = fscanf(infile, "TerminalVelocityLimit %f\n", &TerminalVelocityLimit);
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "DISSOLUTION PARAMETERS\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "concen_diff %lf\n", &concen_diff);
		fret = fscanf(infile, "concen_atm %lf\n", &concen_atm);
	#else // single
		fret = fscanf(infile, "concen_diff %f\n", &concen_diff);
		fret = fscanf(infile, "concen_atm %f\n", &concen_atm);
	#endif
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "BUBBLE GENERATOR\n");
	#ifdef DOUBLE
		fret = fscanf(infile, "(BGis, BGjs, BGks)");
		fret = fscanf(infile, " %d", &BubbleGenerator.BGis);
		fret = fscanf(infile, " %d", &BubbleGenerator.BGjs);
		fret = fscanf(infile, " %d\n", &BubbleGenerator.BGks);
		fret = fscanf(infile, "(BGie, BGje, BGke)");
		fret = fscanf(infile, " %d", &BubbleGenerator.BGie);
		fret = fscanf(infile, " %d", &BubbleGenerator.BGje);
		fret = fscanf(infile, " %d\n", &BubbleGenerator.BGke);
		fret = fscanf(infile, "BubGen_type");
		fret = fscanf(infile, " %s", buf);
		if(strcmp(buf, "GAUSSIAN") == 0) {
			BubbleGenerator.BubGen_type = GAUSSIAN;
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_amplitude);
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_sigmaX);
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_sigmaY);
			fret = fscanf(infile, " %lf\n", &BubbleGenerator.BubGen_sigmaZ);
			fret = fscanf(infile, "BubGen_dia %lf\n", &BubbleGenerator.BubGen_dia);
		} else if(strcmp(buf, "HYPERBOLICTAN") == 0) {
			BubbleGenerator.BubGen_type = HYPERBOLICTAN;
			fret = fscanf(infile, " %lf\n", &BubbleGenerator.BubGen_bubblenumber);
			fret = fscanf(infile, "(Lx1, Lx2, epsx)");
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_Lx1);
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_Lx2);
			fret = fscanf(infile, " %lf\n", &BubbleGenerator.BubGen_epsx);
			fret = fscanf(infile, "(Ly1, Ly2, epsy)");
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_Ly1);
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_Ly2);
			fret = fscanf(infile, " %lf\n", &BubbleGenerator.BubGen_epsy);
			fret = fscanf(infile, "(z0, sigmaZ)");
			fret = fscanf(infile, " %lf", &BubbleGenerator.BubGen_z0);
			fret = fscanf(infile, " %lf\n", &BubbleGenerator.BubGen_sigmaZ);
			fret = fscanf(infile, "BubGen_dia %lf\n", &BubbleGenerator.BubGen_dia);
		} else {
			fprintf(stderr, "Eulerian.config read error.\n");
			exit(EXIT_FAILURE);
		}
	#else // single precision
		fret = fscanf(infile, "(BGis, BGjs, BGks)");
		fret = fscanf(infile, " %d", &BubbleGenerator.BGis);
		fret = fscanf(infile, " %d", &BubbleGenerator.BGjs);
		fret = fscanf(infile, " %d\n", &BubbleGenerator.BGks);
		fret = fscanf(infile, "(BGie, BGje, BGke)");
		fret = fscanf(infile, " %d", &BubbleGenerator.BGie);
		fret = fscanf(infile, " %d", &BubbleGenerator.BGje);
		fret = fscanf(infile, " %d\n", &BubbleGenerator.BGke);
		fret = fscanf(infile, "BubGen_type");
		fret = fscanf(infile, " %s", buf);
		if(strcmp(buf, "GAUSSIAN") == 0) {
			BubbleGenerator.BubGen_type = GAUSSIAN;
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_amplitude);
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_sigmaX);
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_sigmaY);
			fret = fscanf(infile, " %f\n", &BubbleGenerator.BubGen_sigmaZ);
			fret = fscanf(infile, "BubGen_dia %f\n", &BubbleGenerator.BubGen_dia);
		} else if(strcmp(buf, "HYPERBOLICTAN") == 0) {
			BubbleGenerator.BubGen_type = HYPERBOLICTAN;
			fret = fscanf(infile, " %f\n", &BubbleGenerator.BubGen_bubblenumber);
			fret = fscanf(infile, "(Lx1, Lx2, epsx)");
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_Lx1);
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_Lx2);
			fret = fscanf(infile, " %f\n", &BubbleGenerator.BubGen_epsx);
			fret = fscanf(infile, "(Ly1, Ly2, epsy)");
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_Ly1);
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_Ly2);
			fret = fscanf(infile, " %f\n", &BubbleGenerator.BubGen_epsy);
			fret = fscanf(infile, "(z0, sigmaZ)");
			fret = fscanf(infile, " %f", &BubbleGenerator.BubGen_z0);
			fret = fscanf(infile, " %f\n", &BubbleGenerator.BubGen_sigmaZ);
			fret = fscanf(infile, "BubGen_dia %f\n", &BubbleGenerator.BubGen_dia);
		} else {
			fprintf(stderr, "Eulerian.config read error.\n");
			exit(EXIT_FAILURE);
		}
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
		fret = fscanf(infile, " %lf %lf", &numdenBC.nBD, &numdenBC.nBDt);
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
		fret = fscanf(infile, " %lf %lf", &bubmasBC.nBD, &bubmasBC.nBDt);
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
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "TURBULENCE DRIVING CONDITIONS\n");
	fret = fscanf(infile, "turb_switch %s\n", buf);
	if(strcmp(buf, "ON") == 0) {
		turb_switch = ON;
	} else if(strcmp(buf, "OFF") == 0) {
		turb_switch = OFF;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	#ifdef DOUBLE
	fret = fscanf(infile, "turbA %lf\n", &turbA);
	#else
	fret = fscanf(infile, "turbA %f\n", &turbA);
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
	bubdia = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
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
	
	// initialize bubble diameter field
	for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
		for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
			for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
				C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
				bubdia[C] = 0;
			}
		}
	}
	
	/*
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
	*/
	
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
	printf("gravity_direction ");
	printf("%d\n", gravity_direction);
	printf("TerminalVelocityLimit ");
	printf("%f\n", TerminalVelocityLimit);
	printf("\n");
	
	printf("DISSOLUTION PARAMETERS\n");
	printf("concen_diff ");
	printf("%f\n", concen_diff);
	printf("concen_atm ");
	printf("%f\n", concen_atm);
	printf("\n");
	
	printf("BUBBLE GENERATOR\n");
	printf("(BGis, BGjs, BGks)");
	printf(" %d", BubbleGenerator.BGis);
	printf(" %d", BubbleGenerator.BGjs);
	printf(" %d\n", BubbleGenerator.BGks);
	printf("(BGie, BGje, BGke)");
	printf(" %d", BubbleGenerator.BGie);
	printf(" %d", BubbleGenerator.BGje);
	printf(" %d\n", BubbleGenerator.BGke);
	printf("BubGen_type ");
	if(BubbleGenerator.BubGen_type == GAUSSIAN) {
		printf("GAUSSIAN");
		printf(" %f", BubbleGenerator.BubGen_amplitude);
		printf(" %f", BubbleGenerator.BubGen_sigmaX);
		printf(" %f", BubbleGenerator.BubGen_sigmaX);
		printf(" %f\n", BubbleGenerator.BubGen_sigmaX);
	} else if(BubbleGenerator.BubGen_type == HYPERBOLICTAN) {
		printf("HYPERBOLICTAN");
		printf(" %f\n", BubbleGenerator.BubGen_bubblenumber);
		printf("(Lx1, Lx2, epsx)");
		printf(" %f", BubbleGenerator.BubGen_Lx1);
		printf(" %f", BubbleGenerator.BubGen_Lx2);
		printf(" %f\n", BubbleGenerator.BubGen_epsx);
		printf("(Ly1, Ly2, epsy)");
		printf(" %f", BubbleGenerator.BubGen_Ly1);
		printf(" %f", BubbleGenerator.BubGen_Ly2);
		printf(" %f\n", BubbleGenerator.BubGen_epsy);
		printf("(z0, sigmaZ)");
		printf(" %f", BubbleGenerator.BubGen_z0);
		printf(" %f\n",  BubbleGenerator.BubGen_sigmaZ);
	}
	printf("BubGen_dia");
	printf(" %f\n", BubbleGenerator.BubGen_dia);
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
	printf("TURBULENCE DRIVING CONDITIONS\n");
	printf("turb_switch %d\n", turb_switch);
	printf("turbA %f\n", turbA);
	
	printf("\n");
	/*
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
	*/
}

void cgns_grid_Eulerian(void)
{
	// create the file
	char fname[FILE_NAME_SIZE];
	sprintf(fname, "%s/output/%s", ROOT_DIR, "grid.cgns");
	
	int fn;
	int bn;
	int zn;
	int gn;
	int cn;
	
	// check that grid.cgns exists, for restart process
	cg_open(fname, CG_MODE_WRITE, &fn);
	
	cg_base_write(fn, "Base", 3, 3, &bn);
	cgsize_t size[9];
	size[0] = Dom.xn+1; // cells -> vertices
	size[1] = Dom.yn+1;
	size[2] = Dom.zn+1;
	size[3] = Dom.xn;
	size[4] = Dom.yn;
	size[5] = Dom.zn;
	size[6] = 0;
	size[7] = 0;
	size[8] = 0;
	//############################################################################
	cg_biter_write(fn, bn, "BaseIterativeData", 0);
	cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
	cg_ziter_write(fn, bn, zn, "ZoneIterativeData");
	//############################################################################
	cg_grid_write(fn, bn, zn, "GridCoordinates", &gn);
	real *x = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
	// cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
	real *y = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
	// cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
	real *z = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
	// cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke+1; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je+1; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie+1; i++) {
				int C = (i-1) + (j-1)*(Dom.xn+1) + (k-1)*(Dom.xn+1)*(Dom.yn+1);
				x[C] = Dom.xs + (i-1)*Dom.dx;
				y[C] = Dom.ys + (j-1)*Dom.dy;
				z[C] = Dom.zs + (k-1)*Dom.dz;
			}
		}
	}
	
	cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &cn);
	cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &cn);
	cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &cn);
	
	free(x);
	free(y);
	free(z);
	
	cg_close(fn);
}

void cgns_flow_field_Eulerian(void)
{
	// create the solution file
	//char fname[FILE_NAME_SIZE];
	//char fname2[FILE_NAME_SIZE];
	//char fnameall[FILE_NAME_SIZE];
	//char fnameall2[FILE_NAME_SIZE];
	char gname[FILE_NAME_SIZE];
	char gnameall[FILE_NAME_SIZE];
	//real tout = ttime; //  = rec_flow_field_stepnum_out * dtout;
	//char format[CHAR_BUF_SIZE];
	char snodename[CHAR_BUF_SIZE];
	//char snodenameall[CHAR_BUF_SIZE];
	//int sigfigs = ceil(log10(1. / dtout));
	//if(sigfigs < 1) sigfigs = 1;
	//sprintf(format, "%%.%df", sigfigs);
	
	//sprintf(fname2, "flow-%s.cgns", format);
	//sprintf(fnameall2, "%s/output/flow-%s.cgns", ROOT_DIR, format);
	sprintf(snodename, "Solution-");
	//sprintf(snodenameall, "/Base/Zone0/Solution-");
	sprintf(snodename, "%s%f", snodename, ttime);
	//sprintf(snodenameall, "%s%s", snodenameall, format);
	//sprintf(fname, fname2, tout);
	//sprintf(fnameall, fnameall2, tout);
	//sprintf(snodename, snodename, ttime);
	//sprintf(snodenameall, snodenameall, tout);

	// grid file name
	sprintf(gname, "grid.cgns");
	// grid file name and path
	sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
	
	int fn;
	// in our case we have only one base and one zone, need to check in future
	int bn = 1;
	int zn = 1;
	int sn;
	int fnpress;
	int fnu;
	int fnv;
	int fnw;
	//int fnwb;
	int fnnumden;
	//int fnbubmas;
	int fnconcen;
	int fnbubdia;

	// check that grid.cgns exists
	if(cg_open(gnameall, CG_MODE_MODIFY, &fn) != 0) {
		fprintf(stderr, "CGNS flow field write failure: no grid.cgns\n");
		exit(EXIT_FAILURE);
	}
	
	// write the solution node for this time step
	cg_sol_write(fn, bn, zn, snodename, CellCenter, &sn);
	
	// write the fields
	//==========================================================================
	real *pout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				pout[C] = p[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Pressure", pout, &fnpress);
	
	real *uout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
		for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
			for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = (i-1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC1 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC2 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC3 = (i+2) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				uout[C] = -0.0625*u[CC0] + 0.5625*u[CC1] + 0.5625*u[CC2] - 0.0625*u[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", uout, &fnu);
	
	real *vout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
		for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
			for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC1 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC2 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC3 = i + (j+2)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				vout[C] = -0.0625*v[CC0] + 0.5625*v[CC1] + 0.5625*v[CC2] - 0.0625*v[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", vout, &fnv);
	
	real *wout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
		for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
			for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
				int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
				int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
				int CC3 = i + j*Dom.Gfz.s1b + (k+2)*Dom.Gfz.s2b;
				wout[C] = -0.0625*w[CC0] + 0.5625*w[CC1] + 0.5625*w[CC2] - 0.0625*w[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", wout, &fnw);
	/*
	real *wbout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
		for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
			for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
				int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
				wbout[C] = 0.5 * (w_b[CC1] + w_b[CC2]);
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "BubbleVelocityZ", wbout, &fnwb);
	*/
	real *numdenout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				numdenout[C] = numden[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "NumberDensity", numdenout, &fnnumden);
	/*
	real *bubmasout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				bubmasout[C] = bubmas[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "BubbleMass", bubmasout, &fnbubmas);
	*/
	real *concenout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				concenout[C] = concen[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Concentration", concenout, &fnconcen);
	
	real *bubdiaout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				bubdiaout[C] = bubdia[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Bubble Diameter", bubdiaout, &fnbubdia);
	//==========================================================================
	
	free(pout);
	free(uout);
	free(vout);
	free(wout);
	//free(wbout);
	free(numdenout);
	//free(bubmasout);
	free(concenout);
	free(bubdiaout);
	
	// now add this timestep into time series
	int N;
	char buf[32];
	cgsize_t tsize[1];
	cgsize_t solsize[2];
	real *cgns_tseries;
	
	// get the number of solutions, remember, the new solution is already written
	cg_nsols(fn, bn, zn, &N);
	
	cgns_tseries = malloc(N * sizeof(real));
	cg_goto(fn, bn, "BaseIterativeData_t", 1, "end");
	cg_array_read(1, cgns_tseries);
	cg_biter_write(fn, bn, "BaseIterativeData", N);
	cg_goto(fn, bn, "BaseIterativeData_t", 1, "end");
	cgns_tseries[N - 1] = ttime;
	tsize[0] = N;
	cg_array_write("TimeValues", RealDouble, 1, tsize, cgns_tseries);
	
	const char *solname[N];
	GridLocation_t location;
	for(int i = 0; i < N - 1; i++) {
		cg_sol_info(fn, bn, zn, i + 1, buf, &location);
		//printf("===%s===\n", buf);
		solname[i] = buf;
	}
	solname[N - 1] = snodename;
	cg_goto(fn, bn, "Zone_t", zn, "ZoneIterativeData_t", 1, "end");
	solsize[0] = 32;
	solsize[1] = N;
	cg_array_write("FlowSolutionPointers", Character, 2, solsize, solname);
	
	cg_simulation_type_write(fn, bn, TimeAccurate);
	
	free(cgns_tseries);
	cg_close(fn);
}
/*
void cgns_finish_Eulerian(void)
{
	int i;
	int N;
	char buf[CHAR_BUF_SIZE];
	char gnameall[FILE_NAME_SIZE];
	int fn;
	// in our case we have only one base and one zone, need to check in future
	int bn = 1;
	int zn = 1;
	cgsize_t tsize[1];
	cgsize_t solsize[2];
	
	// grid file name and path
	sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
	cg_open(gnameall, CG_MODE_MODIFY, &fn);
	
	// get the number of solutions
	cg_nsols(fn, bn, zn, &N);
	
	cgns_tseries = malloc(N * sizeof(real));
	
	const char *cgns_solname[N];
	
	int fret = 0;
	fret = fret; // prevent compiler warning
	
	// open soluton name record file for reading
	char fname[FILE_NAME_SIZE];
	sprintf(fname, "%s/record/solname.rec", ROOT_DIR);
	FILE *infile = fopen(fname, "r");
	if(infile == NULL) {
		fprintf(stderr, "Could not open file %s\n", fname);
		exit(EXIT_FAILURE);
	}
	
	fret = fscanf(infile, "Solution names:\n");
	for(i = 0; i < N; i++) {
		fret = fscanf(infile, "%s\n", buf);
		cgns_solname[i] = buf;
		fret = fscanf(infile, "%lf\n", &cgns_tseries[i]);
		printf("+++%s, %lf\n", cgns_solname[i], cgns_tseries[i]);
	}
	
	fclose(infile);
	
	tsize[0] = N;
	cg_biter_write(fn, bn, "BaseIterativeData", N);
	cg_goto(fn, bn, "BaseIterativeData_t", 1, "end");
	cg_array_write("TimeValues", RealDouble, 1, tsize, cgns_tseries);
	cg_goto(fn, bn, "Zone_t", zn, "ZoneIterativeData_t", 1, "end");
	solsize[0] = 32;
	solsize[1] = N;
	cgns_solname[0] = "Solution-0.0000";
	cgns_solname[1] = "Solution-0.0022";
	cgns_solname[2] = "Solution-0.0044";
	cg_array_write("FlowSolutionPointers", Character, 2, solsize, cgns_solname);
	
	cg_simulation_type_write(fn, bn, TimeAccurate);
	
	free(cgns_tseries);
	cg_close(fn);
}
*/


int Eulerian_init_parameters(void)
{
	int i, j, k, C, CC;  // iterators
	real h; // depth
	
	// make sure there are enough GPU devices in the given range
	if(nsubdom > dev_end - dev_start + 1) {
		return EXIT_FAILURE;
	}
	
	// allocate memory
	// bubble density at different depth
	bubden = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	// numGau is used for a Gaussian number density BC
	numGau = (real*) malloc(Dom.Gfz.s2b * sizeof(real));
	cpumem += Dom.Gfz.s2b * sizeof(real);
	masGau = (real*) malloc(Dom.Gfz.s2b * sizeof(real));
	cpumem += Dom.Gfz.s2b * sizeof(real);
	// saturation concentration at different depth
	concen_sat = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	// bubble generator(source term in number density equation)
	BGndot = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	BGmdot = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	
	/*
	// Gaussian BC
	if(bubble_init_cond == UNIFORM) {
		// Gasussian number density
		for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
			for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
				C = i + j * Dom.Gcc.s1b;
				real xx = Dom.xs+Dom.dx*((real)(i-Dom.Gfz.is)+0.5);
				real yy = Dom.ys+Dom.dy*((real)(j-Dom.Gfz.js)+0.5);
				numGau[C] = numdenBC.nBD*exp(-0.5*(pow(xx,2.0)+pow(yy,2.0)));
			}
		}
		// Periodic BC in X and Y
		for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
			C  = i + Dom.Gfz.jsb    * Dom.Gfz.s1b;
			CC = i + (Dom.Gfz.je-1) * Dom.Gfz.s1b;
			numGau[C] = numGau[CC];
			C  = i + (Dom.Gfz.jeb-1) * Dom.Gfz.s1b;
			CC = i + Dom.Gfz.js      * Dom.Gfz.s1b;
			numGau[C] = numGau[CC];
		}
		for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
			C  = Dom.Gfz.isb    + j * Dom.Gfz.s1b;
			CC = (Dom.Gfz.ie-1) + j * Dom.Gfz.s1b;
			numGau[C] = numGau[CC];
			C  = (Dom.Gfz.ieb-1) + j * Dom.Gfz.s1b;
			CC = Dom.Gfz.is      + j * Dom.Gfz.s1b;
			numGau[C] = numGau[CC];
		}
		// Gaussian bubble mass
		for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
			for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
				C = i + j * Dom.Gcc.s1b;
				masGau[C] = bubmasBC.nBD / numdenBC.nBD * numGau[C];
			}
		}
	}
	*/
	
	// initialize cell-centered and face-centered bubble density field, 
	// need to define the gravity direction, note that POSITIVE direction have 
	// to be the water surface and NEGATIVE direction have to be the bottom.
	if(gravity_direction == GRAVITY_X) {
		/*
		bubden_face = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
		cpumem += Dom.Gfx.s3b * sizeof(real);
		
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb - 1; i++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					h = ((real)(Dom.Gcc.ieb - i) - 1.5) * Dom.dx;
					bubden[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
					i = Dom.Gcc.ieb - 1;
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					bubden[C] = rho_atm;
			}
		}
		for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
			for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
				for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb - 1; i++) {
					C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
					h = ((real)(Dom.Gfx.ieb - i) - 2.0) * Dom.dx;
					bubden_face[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
			for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
					i = Dom.Gfx.ieb - 1;
					C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
					bubden_face[C] = rho_atm;
			}
		}
		*/
	}
	else if(gravity_direction == GRAVITY_Y) {
		/*
		bubden_face = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
		cpumem += Dom.Gfy.s3b * sizeof(real);
		
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb - 1; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					h = ((real)(Dom.Gcc.jeb - j) - 1.5) * Dom.dy;
					bubden[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
				j = Dom.Gcc.jeb - 1;
				C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
				bubden[C] = rho_atm;
			}
		}
		for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
			for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb - 1; j++) {
				for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
					C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
					h = ((real)(Dom.Gfy.ieb - i) - 2.0) * Dom.dy;
					bubden_face[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
			for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
				j = Dom.Gfy.jeb - 1;
				C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
				bubden_face[C] = rho_atm;
			}
		}
		*/
	}
	else if(gravity_direction == GRAVITY_Z) {
		
		// Cell-centered bubble density, the values in ghost cells are equal to
		// those in neighbor cells, this is good for N-BC
		// And saturation concentration at different depth
		for(k = Dom.Gcc.ksb + 1; k < Dom.Gcc.keb - 1; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					h = ((real)(Dom.Gcc.keb - k) - 1.5) * Dom.dz;
					bubden[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
					concen_sat[C] = concen_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
			for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
				k = Dom.Gcc.keb - 1;
				C  = i + j * Dom.Gcc.s1b + k       * Dom.Gcc.s2b;
				CC = i + j * Dom.Gcc.s1b + (k - 1) * Dom.Gcc.s2b;
				bubden[C] = bubden[CC];
				concen_sat[C] = concen_sat[CC];
			}
		}
		for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
			for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
				k = Dom.Gcc.ksb;
				C  = i + j * Dom.Gcc.s1b + k       * Dom.Gcc.s2b;
				CC = i + j * Dom.Gcc.s1b + (k + 1) * Dom.Gcc.s2b;
				bubden[C] = bubden[CC];
				concen_sat[C] = concen_sat[CC];
			}
		}
		
		// Face-centered bubble density, the values in ghost cells are equal to
		// those in neighbor cells
		bubden_face = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
		cpumem += Dom.Gfz.s3b * sizeof(real);
		
		for(k = Dom.Gfz.ksb + 1; k < Dom.Gfz.keb - 1; k++) {
			for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
				for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
					C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
					h = ((real)(Dom.Gfz.keb - k) - 2.0) * Dom.dz;
					bubden_face[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
			for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
				k = Dom.Gfz.keb - 1;
				C  = i + j * Dom.Gfz.s1b + k       * Dom.Gfz.s2b;
				CC = i + j * Dom.Gfz.s1b + (k - 1) * Dom.Gfz.s2b;
				bubden_face[C] = bubden_face[CC];
			}
		}
		for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
			for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
				k = Dom.Gfz.ksb;
				C  = i + j * Dom.Gfz.s1b + k       * Dom.Gfz.s2b;
				CC = i + j * Dom.Gfz.s1b + (k + 1) * Dom.Gfz.s2b;
				bubden_face[C] = bubden_face[CC];
			}
		}
	}
	
	// initialize the bubble generator after density field being initialized
	if(BubbleGenerator.BubGen_type == GAUSSIAN) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					BGndot[C] = 0.;
				}
			}
		}
		for(i = BubbleGenerator.BGis; i < BubbleGenerator.BGie; i++) {
			for(j = BubbleGenerator.BGjs; j < BubbleGenerator.BGje; j++) {
				for(k = BubbleGenerator.BGks; k < BubbleGenerator.BGke; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					real xx = (i - 0.5) * Dom.dx + Dom.xs;
					real yy = (j - 0.5) * Dom.dy + Dom.ys;
					real zz = (k - 0.5) * Dom.dz + Dom.zs;
					real XX = - pow((xx - BubbleGenerator.BubGen_x0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaX, 2.0);
					real YY = - pow((yy - BubbleGenerator.BubGen_y0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaY, 2.0);
					real ZZ = - pow((zz - BubbleGenerator.BubGen_z0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaZ, 2.0);
					BGndot[C] = BubbleGenerator.BubGen_amplitude * exp(XX + YY + ZZ);
				}
			}
		}
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					BGmdot[C] = PI / 6 * pow(BubbleGenerator.BubGen_dia, 3.0) * bubden[C] * BGndot[C];
				}
			}
		}
	} else if(BubbleGenerator.BubGen_type == HYPERBOLICTAN) {
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					BGndot[C] = 0.;
				}
			}
		}
		if(BubbleGenerator.BubGen_Lx1 >= BubbleGenerator.BubGen_Lx2 || BubbleGenerator.BubGen_Ly1 >= BubbleGenerator.BubGen_Ly2) {
			fprintf(stderr, "Bubble generator L1 must be smaller than L2.\n");
			exit(EXIT_FAILURE);
		}
		// calculate the summation of number density, so we can know the amplitude
		real summation = 0.0;
		for(i = BubbleGenerator.BGis; i < BubbleGenerator.BGie; i++) {
			for(j = BubbleGenerator.BGjs; j < BubbleGenerator.BGje; j++) {
				for(k = BubbleGenerator.BGks; k < BubbleGenerator.BGke; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					real xx = (i - 0.5) * Dom.dx + Dom.xs;
					real yy = (j - 0.5) * Dom.dy + Dom.ys;
					real zz = (k - 0.5) * Dom.dz + Dom.zs;
					real XX = 0.5*(1+tanh((xx-BubbleGenerator.BubGen_Lx1)/BubbleGenerator.BubGen_epsx));
					XX = XX * 0.5*(1+tanh((BubbleGenerator.BubGen_Lx2-xx)/BubbleGenerator.BubGen_epsx));
					real YY = 0.5*(1+tanh((yy-BubbleGenerator.BubGen_Ly1)/BubbleGenerator.BubGen_epsy));
					YY = YY * 0.5*(1+tanh((BubbleGenerator.BubGen_Ly2-yy)/BubbleGenerator.BubGen_epsy));
					real ZZ = - pow((zz - BubbleGenerator.BubGen_z0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaZ, 2.0);
					BGndot[C] = XX * YY * exp(ZZ);
					summation = summation + BGndot[C];
				}
			}
		}
		BubbleGenerator.BubGen_amplitude = BubbleGenerator.BubGen_bubblenumber / (Dom.dx * Dom.dy * Dom.dz * summation);
		for(i = BubbleGenerator.BGis; i < BubbleGenerator.BGie; i++) {
			for(j = BubbleGenerator.BGjs; j < BubbleGenerator.BGje; j++) {
				for(k = BubbleGenerator.BGks; k < BubbleGenerator.BGke; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					BGndot[C] = BubbleGenerator.BubGen_amplitude * BGndot[C];
				}
			}
		}
		for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					BGmdot[C] = PI / 6 * pow(BubbleGenerator.BubGen_dia, 3.0) * bubden[C] * BGndot[C];
				}
			}
		}
	} else {
		fprintf(stderr, "Bubble generator type incorrect.\n");
		exit(EXIT_FAILURE);
	}
	
	return EXIT_SUCCESS;
}

void in_restart_Eulerian(void)
{
  int i, j, k;  // iterators
  int fret = 0;
  fret = fret; // prevent compiler warning
  // open configuration file for reading
  char fname[FILE_NAME_SIZE];
  sprintf(fname, "%s/input/restart.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // read domain
  // flow solver data
#ifdef DOUBLE
  fret = fscanf(infile, "%le %le %le %d %d\n", &ttime, &dt0, &dt, &stepnum,
    &rec_paraview_stepnum_out);
#else
  fret = fscanf(infile, "%e %e %e %d %d\n", &ttime, &dt0, &dt, &stepnum,
    &rec_paraview_stepnum_out);
#endif

  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
#ifdef DOUBLE
        fret = fscanf(infile, "%le %le %le %le %le %le %le ",
          &u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &u0[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &diff0_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &conv0_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &diff_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &conv_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &u_star[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
#else
        fret = fscanf(infile, "%e %e %e %e %e %e %e ",
          &u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &u0[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &diff0_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &conv0_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &diff_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &conv_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          &u_star[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
#endif
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
#ifdef DOUBLE
        fret = fscanf(infile, "%le %le %le %le %le %le %le ",
          &v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &v0[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &diff0_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &conv0_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &diff_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &conv_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &v_star[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b]);
#else
        fret = fscanf(infile, "%e %e %e %e %e %e %e ",
          &v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &v0[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &diff0_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &conv0_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &diff_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &conv_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          &v_star[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b]);
#endif
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
#ifdef DOUBLE
        fret = fscanf(infile, "%le %le %le %le %le %le %le ",
          &w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &w0[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &diff0_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &conv0_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &diff_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &conv_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &w_star[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b]);
#else
        fret = fscanf(infile, "%e %e %e %e %e %e %e ",
          &w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &w0[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &diff0_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &conv0_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &diff_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &conv_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          &w_star[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b]);
#endif
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        #ifdef DOUBLE
        fret = fscanf(infile, "%le ", &p[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%le ", &p0[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        #else
        fret = fscanf(infile, "%e ", &p[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%e ", &p0[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        #endif
        fret = fscanf(infile, "%d ", &phase[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%d ", &phase_shell[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        //######################################################################
        #ifdef DOUBLE
        fret = fscanf(infile, "%le ", &numden[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%le ", &bubmas[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%le ", &concen[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        #else
        fret = fscanf(infile, "%e ", &numden[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%e ", &bubmas[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fret = fscanf(infile, "%e ", &concen[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        #endif
        //######################################################################
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        fret = fscanf(infile, "%d ", &flag_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        fret = fscanf(infile, "%d ", &flag_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b]);
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        fret = fscanf(infile, "%d ", &flag_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b]);
      }
    }
  }
  fret = fscanf(infile, "\n");

  for(i = 0; i < 9; i++) {
#ifdef DOUBLE
    fret = fscanf(infile, "%le ", &bc_plane_pos[i]);
#else
    fret = fscanf(infile, "%e ", &bc_plane_pos[i]);
#endif
  }
  fprintf(infile, "\n");

  for(i = 0; i < nparts; i++) {
#ifdef DOUBLE
    fret = fscanf(infile, "%le ", &parts[i].r);
    fret = fscanf(infile, "%le ", &parts[i].x);
    fret = fscanf(infile, "%le ", &parts[i].y);
    fret = fscanf(infile, "%le ", &parts[i].z);
    fret = fscanf(infile, "%le ", &parts[i].u);
    fret = fscanf(infile, "%le ", &parts[i].v);
    fret = fscanf(infile, "%le ", &parts[i].w);
    fret = fscanf(infile, "%le ", &parts[i].u0);
    fret = fscanf(infile, "%le ", &parts[i].v0);
    fret = fscanf(infile, "%le ", &parts[i].w0);
    fret = fscanf(infile, "%le ", &parts[i].udot);
    fret = fscanf(infile, "%le ", &parts[i].vdot);
    fret = fscanf(infile, "%le ", &parts[i].wdot);
    fret = fscanf(infile, "%le ", &parts[i].axx);
    fret = fscanf(infile, "%le ", &parts[i].axy);
    fret = fscanf(infile, "%le ", &parts[i].axz);
    fret = fscanf(infile, "%le ", &parts[i].ayx);
    fret = fscanf(infile, "%le ", &parts[i].ayy);
    fret = fscanf(infile, "%le ", &parts[i].ayz);
    fret = fscanf(infile, "%le ", &parts[i].azx);
    fret = fscanf(infile, "%le ", &parts[i].azy);
    fret = fscanf(infile, "%le ", &parts[i].azz);
    fret = fscanf(infile, "%le ", &parts[i].ox);
    fret = fscanf(infile, "%le ", &parts[i].oy);
    fret = fscanf(infile, "%le ", &parts[i].oz);
    fret = fscanf(infile, "%le ", &parts[i].ox0);
    fret = fscanf(infile, "%le ", &parts[i].oy0);
    fret = fscanf(infile, "%le ", &parts[i].oz0);
    fret = fscanf(infile, "%le ", &parts[i].oxdot);
    fret = fscanf(infile, "%le ", &parts[i].oydot);
    fret = fscanf(infile, "%le ", &parts[i].ozdot);
    fret = fscanf(infile, "%le ", &parts[i].Fx);
    fret = fscanf(infile, "%le ", &parts[i].Fy);
    fret = fscanf(infile, "%le ", &parts[i].Fz);
    fret = fscanf(infile, "%le ", &parts[i].Lx);
    fret = fscanf(infile, "%le ", &parts[i].Ly);
    fret = fscanf(infile, "%le ", &parts[i].Lz);
    fret = fscanf(infile, "%le ", &parts[i].aFx);
    fret = fscanf(infile, "%le ", &parts[i].aFy);
    fret = fscanf(infile, "%le ", &parts[i].aFz);
    fret = fscanf(infile, "%le ", &parts[i].aLx);
    fret = fscanf(infile, "%le ", &parts[i].aLy);
    fret = fscanf(infile, "%le ", &parts[i].aLz);
    fret = fscanf(infile, "%le ", &parts[i].rho);
    for(j = 0; j<NNODES; j++){    
    fret = fscanf(infile, "%d ", &parts[i].nodes[j]);
    }
    fret = fscanf(infile, "%d ", &parts[i].order);
    fret = fscanf(infile, "%le ", &parts[i].rs);
    fret = fscanf(infile, "%d ", &parts[i].ncoeff);
    fret = fscanf(infile, "%le ", &parts[i].spring_k);
    fret = fscanf(infile, "%le ", &parts[i].spring_x);
    fret = fscanf(infile, "%le ", &parts[i].spring_y);
    fret = fscanf(infile, "%le ", &parts[i].spring_z);
    fret = fscanf(infile, "%le ", &parts[i].spring_l);
    fret = fscanf(infile, "%d ", &parts[i].translating);
    fret = fscanf(infile, "%d ", &parts[i].rotating);
    fret = fscanf(infile, "%d ", &parts[i].cage.cx);
    fret = fscanf(infile, "%d ", &parts[i].cage.cy);
    fret = fscanf(infile, "%d ", &parts[i].cage.cz);
    fret = fscanf(infile, "%d ", &parts[i].cage.is);
    fret = fscanf(infile, "%d ", &parts[i].cage.ibs);
    fret = fscanf(infile, "%d ", &parts[i].cage.ie);
    fret = fscanf(infile, "%d ", &parts[i].cage.ibe);
    fret = fscanf(infile, "%d ", &parts[i].cage.in);
    fret = fscanf(infile, "%d ", &parts[i].cage.js);
    fret = fscanf(infile, "%d ", &parts[i].cage.jbs);
    fret = fscanf(infile, "%d ", &parts[i].cage.je);
    fret = fscanf(infile, "%d ", &parts[i].cage.jbe);
    fret = fscanf(infile, "%d ", &parts[i].cage.jn);
    fret = fscanf(infile, "%d ", &parts[i].cage.ks);
    fret = fscanf(infile, "%d ", &parts[i].cage.kbs);
    fret = fscanf(infile, "%d ", &parts[i].cage.ke);
    fret = fscanf(infile, "%d ", &parts[i].cage.kbe);
    fret = fscanf(infile, "%d ", &parts[i].cage.kn);
#else
    fret = fscanf(infile, "%e ", &parts[i].r);
    fret = fscanf(infile, "%e ", &parts[i].x);
    fret = fscanf(infile, "%e ", &parts[i].y);
    fret = fscanf(infile, "%e ", &parts[i].z);
    fret = fscanf(infile, "%e ", &parts[i].u);
    fret = fscanf(infile, "%e ", &parts[i].v);
    fret = fscanf(infile, "%e ", &parts[i].w);
    fret = fscanf(infile, "%e ", &parts[i].u0);
    fret = fscanf(infile, "%e ", &parts[i].v0);
    fret = fscanf(infile, "%e ", &parts[i].w0);
    fret = fscanf(infile, "%e ", &parts[i].udot);
    fret = fscanf(infile, "%e ", &parts[i].vdot);
    fret = fscanf(infile, "%e ", &parts[i].wdot);
    fret = fscanf(infile, "%e ", &parts[i].axx);
    fret = fscanf(infile, "%e ", &parts[i].axy);
    fret = fscanf(infile, "%e ", &parts[i].axz);
    fret = fscanf(infile, "%e ", &parts[i].ayx);
    fret = fscanf(infile, "%e ", &parts[i].ayy);
    fret = fscanf(infile, "%e ", &parts[i].ayz);
    fret = fscanf(infile, "%e ", &parts[i].azx);
    fret = fscanf(infile, "%e ", &parts[i].azy);
    fret = fscanf(infile, "%e ", &parts[i].azz);
    fret = fscanf(infile, "%e ", &parts[i].ox);
    fret = fscanf(infile, "%e ", &parts[i].oy);
    fret = fscanf(infile, "%e ", &parts[i].oz);
    fret = fscanf(infile, "%e ", &parts[i].ox0);
    fret = fscanf(infile, "%e ", &parts[i].oy0);
    fret = fscanf(infile, "%e ", &parts[i].oz0);
    fret = fscanf(infile, "%e ", &parts[i].oxdot);
    fret = fscanf(infile, "%e ", &parts[i].oydot);
    fret = fscanf(infile, "%e ", &parts[i].ozdot);
    fret = fscanf(infile, "%e ", &parts[i].Fx);
    fret = fscanf(infile, "%e ", &parts[i].Fy);
    fret = fscanf(infile, "%e ", &parts[i].Fz);
    fret = fscanf(infile, "%e ", &parts[i].Lx);
    fret = fscanf(infile, "%e ", &parts[i].Ly);
    fret = fscanf(infile, "%e ", &parts[i].Lz);
    fret = fscanf(infile, "%e ", &parts[i].aFx);
    fret = fscanf(infile, "%e ", &parts[i].aFy);
    fret = fscanf(infile, "%e ", &parts[i].aFz);
    fret = fscanf(infile, "%e ", &parts[i].aLx);
    fret = fscanf(infile, "%e ", &parts[i].aLy);
    fret = fscanf(infile, "%e ", &parts[i].aLz);
    fret = fscanf(infile, "%e ", &parts[i].rho);
    for(j = 0; j<NNODES; j++){    
    fret = fscanf(infile, "%d ", &parts[i].nodes[j]);
    }
    fret = fscanf(infile, "%d ", &parts[i].order);
    fret = fscanf(infile, "%e ", &parts[i].rs);
    fret = fscanf(infile, "%d ", &parts[i].ncoeff);
    fret = fscanf(infile, "%e ", &parts[i].spring_k);
    fret = fscanf(infile, "%e ", &parts[i].spring_x);
    fret = fscanf(infile, "%e ", &parts[i].spring_y);
    fret = fscanf(infile, "%e ", &parts[i].spring_z);
    fret = fscanf(infile, "%e ", &parts[i].spring_l);
    fret = fscanf(infile, "%d ", &parts[i].translating);
    fret = fscanf(infile, "%d ", &parts[i].rotating);
    fret = fscanf(infile, "%d ", &parts[i].cage.cz);
    fret = fscanf(infile, "%d ", &parts[i].cage.is);
    fret = fscanf(infile, "%d ", &parts[i].cage.ibs);
    fret = fscanf(infile, "%d ", &parts[i].cage.ie);
    fret = fscanf(infile, "%d ", &parts[i].cage.ibe);
    fret = fscanf(infile, "%d ", &parts[i].cage.in);
    fret = fscanf(infile, "%d ", &parts[i].cage.js);
    fret = fscanf(infile, "%d ", &parts[i].cage.jbs);
    fret = fscanf(infile, "%d ", &parts[i].cage.je);
    fret = fscanf(infile, "%d ", &parts[i].cage.jbe);
    fret = fscanf(infile, "%d ", &parts[i].cage.jn);
    fret = fscanf(infile, "%d ", &parts[i].cage.ks);
    fret = fscanf(infile, "%d ", &parts[i].cage.kbs);
    fret = fscanf(infile, "%d ", &parts[i].cage.ke);
    fret = fscanf(infile, "%d ", &parts[i].cage.kbe);
    fret = fscanf(infile, "%d ", &parts[i].cage.kn);
#endif
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "%d\n", &coeff_stride);

  for(i = 0; i < nparts; i++) {
    for(j = 0; j < coeff_stride; j++) {
#ifdef DOUBLE
      fret = fscanf(infile, "%le ", &pnm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &pnm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &pnm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &pnm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &pnm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &pnm_im00[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &phinm_im00[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%le ", &chinm_im00[j + i*coeff_stride]);

#else
      fret = fscanf(infile, "%e ", &pnm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &pnm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &pnm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &pnm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &pnm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &pnm_im00[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &phinm_im00[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_re[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_im[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_re0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_im0[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_re00[j + i*coeff_stride]);
      fret = fscanf(infile, "%e ", &chinm_im00[j + i*coeff_stride]);
#endif
    }
  }

#ifdef DOUBLE
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "%le %le %le %le %le\n", &rec_flow_field_ttime_out,
    &rec_paraview_ttime_out, &rec_particle_ttime_out, &rec_restart_ttime_out,
    &rec_precursor_ttime_out);
#else
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "%e %e %e %e %e\n", &rec_flow_field_ttime_out,
    &rec_paraview_ttime_out, &rec_particle_ttime_out, &rec_restart_ttime_out,
    &rec_precursor_ttime_out);
#endif

#ifdef DOUBLE
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "%le %le %le\n", &pid_int, &pid_back, &gradP.z);
#else
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "%e %e %e\n", &pid_int, &pid_back, &gradP.z);
#endif

  // close file
  fclose(infile);
}

void out_restart_Eulerian(void)
{
  int i, j, k;  // iterators

  // create the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/input/restart.config", ROOT_DIR);
  FILE *rest = fopen(path, "w");
  if(rest == NULL) {
    fprintf(stderr, "Could not open file restart.input.\n");
    exit(EXIT_FAILURE);
  }

  // write current timestep information (uvw0 is previous timestep info)
  // flow solver data
  fprintf(rest, "%e %e %e %d %d\n", ttime, dt0, dt, stepnum, rec_paraview_stepnum_out);

  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        fprintf(rest, "%e %e %e %e %e %e %e ",
          u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          u0[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          diff0_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          conv0_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          diff_u[i+j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          conv_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b],
          u_star[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        fprintf(rest, "%e %e %e %e %e %e %e ",
          v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          v0[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          diff0_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          conv0_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          diff_v[i+j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          conv_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b],
          v_star[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        fprintf(rest, "%e %e %e %e %e %e %e ",
          w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          w0[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          diff0_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b], 
          conv0_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          diff_w[i+j*Dom.Gfz.s1b + k*Dom.Gfz.s2b], 
          conv_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b],
          w_star[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        fprintf(rest, "%e ", p[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fprintf(rest, "%e ", p0[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fprintf(rest, "%d ", phase[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fprintf(rest, "%d ", phase_shell[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        //######################################################################
        fprintf(rest, "%e ", numden[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fprintf(rest, "%e ", bubmas[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        fprintf(rest, "%e ", concen[i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b]);
        //######################################################################
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        fprintf(rest, "%d ", flag_u[i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        fprintf(rest, "%d ", flag_v[i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        fprintf(rest, "%d ", flag_w[i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b]);
      }
    }
  }
  fprintf(rest, "\n");

  for(i = 0; i < 9; i++) {
    fprintf(rest, "%e ", bc_plane_pos[i]);
  }
  fprintf(rest, "\n");

  for(i = 0; i < nparts; i++) {
    fprintf(rest, "%e ", parts[i].r);
    fprintf(rest, "%e ", parts[i].x);
    fprintf(rest, "%e ", parts[i].y);
    fprintf(rest, "%e ", parts[i].z);
    fprintf(rest, "%e ", parts[i].u);
    fprintf(rest, "%e ", parts[i].v);
    fprintf(rest, "%e ", parts[i].w);
    fprintf(rest, "%e ", parts[i].u0);
    fprintf(rest, "%e ", parts[i].v0);
    fprintf(rest, "%e ", parts[i].w0);
    fprintf(rest, "%e ", parts[i].udot);
    fprintf(rest, "%e ", parts[i].vdot);
    fprintf(rest, "%e ", parts[i].wdot);
    fprintf(rest, "%e ", parts[i].axx);
    fprintf(rest, "%e ", parts[i].axy);
    fprintf(rest, "%e ", parts[i].axz);
    fprintf(rest, "%e ", parts[i].ayx);
    fprintf(rest, "%e ", parts[i].ayy);
    fprintf(rest, "%e ", parts[i].ayz);
    fprintf(rest, "%e ", parts[i].azx);
    fprintf(rest, "%e ", parts[i].azy);
    fprintf(rest, "%e ", parts[i].azz);
    fprintf(rest, "%e ", parts[i].ox);
    fprintf(rest, "%e ", parts[i].oy);
    fprintf(rest, "%e ", parts[i].oz);
    fprintf(rest, "%e ", parts[i].ox0);
    fprintf(rest, "%e ", parts[i].oy0);
    fprintf(rest, "%e ", parts[i].oz0);
    fprintf(rest, "%e ", parts[i].oxdot);
    fprintf(rest, "%e ", parts[i].oydot);
    fprintf(rest, "%e ", parts[i].ozdot);
    fprintf(rest, "%e ", parts[i].Fx);
    fprintf(rest, "%e ", parts[i].Fy);
    fprintf(rest, "%e ", parts[i].Fz);
    fprintf(rest, "%e ", parts[i].Lx);
    fprintf(rest, "%e ", parts[i].Ly);
    fprintf(rest, "%e ", parts[i].Lz);
    fprintf(rest, "%e ", parts[i].aFx);
    fprintf(rest, "%e ", parts[i].aFy);
    fprintf(rest, "%e ", parts[i].aFz);
    fprintf(rest, "%e ", parts[i].aLx);
    fprintf(rest, "%e ", parts[i].aLy);
    fprintf(rest, "%e ", parts[i].aLz);
    fprintf(rest, "%e ", parts[i].rho);
    for(j = 0; j<NNODES; j++){    
      fprintf(rest, "%d ", parts[i].nodes[j]);
    }
    fprintf(rest, "%d ", parts[i].order);
    fprintf(rest, "%e ", parts[i].rs);
    fprintf(rest, "%d ", parts[i].ncoeff);
    fprintf(rest, "%e ", parts[i].spring_k);
    fprintf(rest, "%e ", parts[i].spring_x);
    fprintf(rest, "%e ", parts[i].spring_y);
    fprintf(rest, "%e ", parts[i].spring_z);
    fprintf(rest, "%e ", parts[i].spring_l);
    fprintf(rest, "%d ", parts[i].translating);
    fprintf(rest, "%d ", parts[i].rotating);
    fprintf(rest, "%d ", parts[i].cage.cx);
    fprintf(rest, "%d ", parts[i].cage.cy);
    fprintf(rest, "%d ", parts[i].cage.cz);
    fprintf(rest, "%d ", parts[i].cage.is);
    fprintf(rest, "%d ", parts[i].cage.ibs);
    fprintf(rest, "%d ", parts[i].cage.ie);
    fprintf(rest, "%d ", parts[i].cage.ibe);
    fprintf(rest, "%d ", parts[i].cage.in);
    fprintf(rest, "%d ", parts[i].cage.js);
    fprintf(rest, "%d ", parts[i].cage.jbs);
    fprintf(rest, "%d ", parts[i].cage.je);
    fprintf(rest, "%d ", parts[i].cage.jbe);
    fprintf(rest, "%d ", parts[i].cage.jn);
    fprintf(rest, "%d ", parts[i].cage.ks);
    fprintf(rest, "%d ", parts[i].cage.kbs);
    fprintf(rest, "%d ", parts[i].cage.ke);
    fprintf(rest, "%d ", parts[i].cage.kbe);
    fprintf(rest, "%d ", parts[i].cage.kn);
  }
  fprintf(rest, "\n");

  fprintf(rest, "%d\n", coeff_stride);

  for(i = 0; i < nparts; i++) {
    for(j = 0; j < coeff_stride; j++) {
      fprintf(rest, "%e ", pnm_re[j + i*coeff_stride]);
      fprintf(rest, "%e ", pnm_im[j + i*coeff_stride]);
      fprintf(rest, "%e ", pnm_re0[j + i*coeff_stride]);
      fprintf(rest, "%e ", pnm_im0[j + i*coeff_stride]);
      fprintf(rest, "%e ", pnm_re00[j + i*coeff_stride]);
      fprintf(rest, "%e ", pnm_im00[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_re[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_im[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_re0[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_im0[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_re00[j + i*coeff_stride]);
      fprintf(rest, "%e ", phinm_im00[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_re[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_im[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_re0[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_im0[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_re00[j + i*coeff_stride]);
      fprintf(rest, "%e ", chinm_im00[j + i*coeff_stride]);
    }
  }

  fprintf(rest, "\n");
  fprintf(rest, "%e %e %e %e %e\n", rec_flow_field_ttime_out,
    rec_paraview_ttime_out, rec_particle_ttime_out, rec_restart_ttime_out,
    rec_precursor_ttime_out);
  fprintf(rest, "\n");
  fprintf(rest, "%e %e %e\n", pid_int, pid_back, gradP.z);

  // close the file
  fclose(rest);
}


int domain_init_Eulerian(void)
{
  int i, j, k;    // iterator
  int C, W, E, S, N, B, T;
  real tmp;

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }

  // calculate domain sizes
  Dom.xl = Dom.xe - Dom.xs;
  Dom.yl = Dom.ye - Dom.ys;
  Dom.zl = Dom.ze - Dom.zs;

  // calculate cell sizes
  Dom.dx = Dom.xl / Dom.xn;
  Dom.dy = Dom.yl / Dom.yn;
  Dom.dz = Dom.zl / Dom.zn;

  // set up grids
  // Gcc
  Dom.Gcc.is = DOM_BUF;
  Dom.Gcc.isb = Dom.Gcc.is - DOM_BUF;
  Dom.Gcc.in = Dom.xn;
  Dom.Gcc.inb = Dom.Gcc.in + 2 * DOM_BUF;
  Dom.Gcc.ie = Dom.Gcc.is + Dom.Gcc.in;
  Dom.Gcc.ieb = Dom.Gcc.ie + DOM_BUF;

  Dom.Gcc.js = DOM_BUF;
  Dom.Gcc.jsb = Dom.Gcc.js - DOM_BUF;
  Dom.Gcc.jn = Dom.yn;
  Dom.Gcc.jnb = Dom.Gcc.jn + 2 * DOM_BUF;
  Dom.Gcc.je = Dom.Gcc.js + Dom.Gcc.jn;
  Dom.Gcc.jeb = Dom.Gcc.je + DOM_BUF;

  Dom.Gcc.ks = DOM_BUF;
  Dom.Gcc.ksb = Dom.Gcc.ks - DOM_BUF;
  Dom.Gcc.kn = Dom.zn;
  Dom.Gcc.knb = Dom.Gcc.kn + 2 * DOM_BUF;
  Dom.Gcc.ke = DOM_BUF + Dom.Gcc.kn;
  Dom.Gcc.keb = Dom.Gcc.ke + DOM_BUF;

  Dom.Gcc.s1 = Dom.Gcc.in;
  Dom.Gcc.s2 = Dom.Gcc.s1 * Dom.Gcc.jn;
  Dom.Gcc.s3 = Dom.Gcc.s2 * Dom.Gcc.kn;
  Dom.Gcc.s1b = Dom.Gcc.inb;
  Dom.Gcc.s2b = Dom.Gcc.s1b * Dom.Gcc.jnb;
  Dom.Gcc.s3b = Dom.Gcc.s2b * Dom.Gcc.knb;

  // Gfx
  Dom.Gfx.is = DOM_BUF;
  Dom.Gfx.isb = Dom.Gfx.is - DOM_BUF;
  Dom.Gfx.in = Dom.xn + 1;
  Dom.Gfx.inb = Dom.Gfx.in + 2 * DOM_BUF;
  Dom.Gfx.ie = Dom.Gfx.is + Dom.Gfx.in;
  Dom.Gfx.ieb = Dom.Gfx.ie + DOM_BUF;

  Dom.Gfx.js = DOM_BUF;
  Dom.Gfx.jsb = Dom.Gfx.js - DOM_BUF;
  Dom.Gfx.jn = Dom.yn;
  Dom.Gfx.jnb = Dom.Gfx.jn + 2 * DOM_BUF;
  Dom.Gfx.je = Dom.Gfx.js + Dom.Gfx.jn;
  Dom.Gfx.jeb = Dom.Gfx.je + DOM_BUF;

  Dom.Gfx.ks = DOM_BUF;
  Dom.Gfx.ksb = Dom.Gfx.ks - DOM_BUF;
  Dom.Gfx.kn = Dom.zn;
  Dom.Gfx.knb = Dom.Gfx.kn + 2 * DOM_BUF;
  Dom.Gfx.ke = Dom.Gfx.ks + Dom.Gfx.kn;
  Dom.Gfx.keb = Dom.Gfx.ke + DOM_BUF;

  Dom.Gfx.s1 = Dom.Gfx.in;
  Dom.Gfx.s2 = Dom.Gfx.s1 * Dom.Gfx.jn;
  Dom.Gfx.s3 = Dom.Gfx.s2 * Dom.Gfx.kn;
  Dom.Gfx.s1b = Dom.Gfx.inb;
  Dom.Gfx.s2b = Dom.Gfx.s1b * Dom.Gfx.jnb;
  Dom.Gfx.s3b = Dom.Gfx.s2b * Dom.Gfx.knb;

  // Gfy
  Dom.Gfy.is = DOM_BUF;
  Dom.Gfy.isb = Dom.Gfy.is - DOM_BUF;
  Dom.Gfy.in = Dom.xn;
  Dom.Gfy.inb = Dom.Gfy.in + 2 * DOM_BUF;
  Dom.Gfy.ie = Dom.Gfy.is + Dom.Gfy.in;
  Dom.Gfy.ieb = Dom.Gfy.ie + DOM_BUF;

  Dom.Gfy.js = DOM_BUF;
  Dom.Gfy.jsb = Dom.Gfy.js - DOM_BUF;
  Dom.Gfy.jn = Dom.yn + 1;
  Dom.Gfy.jnb = Dom.Gfy.jn + 2 * DOM_BUF;
  Dom.Gfy.je = Dom.Gfy.js + Dom.Gfy.jn;
  Dom.Gfy.jeb = Dom.Gfy.je + DOM_BUF;

  Dom.Gfy.ks = DOM_BUF;
  Dom.Gfy.ksb = Dom.Gfy.ks - DOM_BUF;
  Dom.Gfy.kn = Dom.zn;
  Dom.Gfy.knb = Dom.Gfy.kn + 2 * DOM_BUF;
  Dom.Gfy.ke = Dom.Gfy.ks + Dom.Gfy.kn;
  Dom.Gfy.keb = Dom.Gfy.ke + DOM_BUF;

  Dom.Gfy.s1 = Dom.Gfy.in;
  Dom.Gfy.s2 = Dom.Gfy.s1 * Dom.Gfy.jn;
  Dom.Gfy.s3 = Dom.Gfy.s2 * Dom.Gfy.kn;
  Dom.Gfy.s1b = Dom.Gfy.inb;
  Dom.Gfy.s2b = Dom.Gfy.s1b * Dom.Gfy.jnb;
  Dom.Gfy.s3b = Dom.Gfy.s2b * Dom.Gfy.knb;

  // Gfz
  Dom.Gfz.is = DOM_BUF;
  Dom.Gfz.isb = Dom.Gfz.is - DOM_BUF;
  Dom.Gfz.in = Dom.xn;
  Dom.Gfz.inb = Dom.Gfz.in + 2 * DOM_BUF;
  Dom.Gfz.ie = Dom.Gfz.is + Dom.Gfz.in;
  Dom.Gfz.ieb = Dom.Gfz.ie + DOM_BUF;

  Dom.Gfz.js = DOM_BUF;
  Dom.Gfz.jsb = Dom.Gfz.js - DOM_BUF;
  Dom.Gfz.jn = Dom.yn;
  Dom.Gfz.jnb = Dom.Gfz.jn + 2 * DOM_BUF;
  Dom.Gfz.je = Dom.Gfz.js + Dom.Gfz.jn;
  Dom.Gfz.jeb = Dom.Gfz.je + DOM_BUF;

  Dom.Gfz.ks = DOM_BUF;
  Dom.Gfz.ksb = Dom.Gfz.ks - DOM_BUF;
  Dom.Gfz.kn = Dom.zn + 1;
  Dom.Gfz.knb = Dom.Gfz.kn + 2 * DOM_BUF;
  Dom.Gfz.ke = Dom.Gfz.ks + Dom.Gfz.kn;
  Dom.Gfz.keb = Dom.Gfz.ke + DOM_BUF;

  Dom.Gfz.s1 = Dom.Gfz.in;
  Dom.Gfz.s2 = Dom.Gfz.s1 * Dom.Gfz.jn;
  Dom.Gfz.s3 = Dom.Gfz.s2 * Dom.Gfz.kn;
  Dom.Gfz.s1b = Dom.Gfz.inb;
  Dom.Gfz.s2b = Dom.Gfz.s1b * Dom.Gfz.jnb;
  Dom.Gfz.s3b = Dom.Gfz.s2b * Dom.Gfz.knb;

  // initialize subdomains
  for(i = 0; i < nsubdom; i++) {
    dom[i].xl = dom[i].xe - dom[i].xs;
    dom[i].yl = dom[i].ye - dom[i].ys;
    dom[i].zl = dom[i].ze - dom[i].zs;
    dom[i].dx = dom[i].xl / dom[i].xn;
    dom[i].dy = dom[i].yl / dom[i].yn;
    dom[i].dz = dom[i].zl / dom[i].zn;

    // TODO: this algorithm will fail if subdomains are not numbered in
    // increasing order.  This will need to be fixed before going to a machine
    // with more than two GPU's.
    // Gcc
    if(dom[i].W > -1)
      dom[i].Gcc.is = dom[dom[i].W].Gcc.ie;
    else
      dom[i].Gcc.is = DOM_BUF;
    dom[i].Gcc.isb = dom[i].Gcc.is - DOM_BUF;
    dom[i].Gcc.in = dom[i].xn;
    dom[i].Gcc.inb = dom[i].Gcc.in + 2 * DOM_BUF;
    dom[i].Gcc.ie = dom[i].Gcc.is + dom[i].Gcc.in;
    dom[i].Gcc.ieb = dom[i].Gcc.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gcc.js = dom[dom[i].S].Gcc.je;
    else
      dom[i].Gcc.js = DOM_BUF;
    dom[i].Gcc.jsb = dom[i].Gcc.js - DOM_BUF;
    dom[i].Gcc.jn = dom[i].yn;
    dom[i].Gcc.jnb = dom[i].Gcc.jn + 2 * DOM_BUF;
    dom[i].Gcc.je = dom[i].Gcc.js + dom[i].Gcc.jn;
    dom[i].Gcc.jeb = dom[i].Gcc.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gcc.ks = dom[dom[i].B].Gcc.ke;
    else
      dom[i].Gcc.ks = DOM_BUF;
    dom[i].Gcc.ksb = dom[i].Gcc.ks - DOM_BUF;
    dom[i].Gcc.kn = dom[i].zn;
    dom[i].Gcc.knb = dom[i].Gcc.kn + 2 * DOM_BUF;
    dom[i].Gcc.ke = DOM_BUF + dom[i].Gcc.kn;
    dom[i].Gcc.keb = dom[i].Gcc.ke + DOM_BUF;

    dom[i].Gcc.s1 = dom[i].Gcc.in;
    dom[i].Gcc.s2 = dom[i].Gcc.s1 * dom[i].Gcc.jn;
    dom[i].Gcc.s3 = dom[i].Gcc.s2 * dom[i].Gcc.kn;
    dom[i].Gcc.s1b = dom[i].Gcc.inb;
    dom[i].Gcc.s2b = dom[i].Gcc.s1b * dom[i].Gcc.jnb;
    dom[i].Gcc.s3b = dom[i].Gcc.s2b * dom[i].Gcc.knb;

    dom[i].Gcc._is = DOM_BUF;
    dom[i].Gcc._isb = dom[i].Gcc._is - DOM_BUF;
    dom[i].Gcc._in = dom[i].xn;
    dom[i].Gcc._inb = dom[i].Gcc._in + 2 * DOM_BUF;
    dom[i].Gcc._ie = dom[i].Gcc._is + dom[i].Gcc._in;
    dom[i].Gcc._ieb = dom[i].Gcc._ie + DOM_BUF;

    dom[i].Gcc._js = DOM_BUF;
    dom[i].Gcc._jsb = dom[i].Gcc._js - DOM_BUF;
    dom[i].Gcc._jn = dom[i].yn;
    dom[i].Gcc._jnb = dom[i].Gcc._jn + 2 * DOM_BUF;
    dom[i].Gcc._je = dom[i].Gcc._js + dom[i].Gcc._jn;
    dom[i].Gcc._jeb = dom[i].Gcc._je + DOM_BUF;

    dom[i].Gcc._ks = DOM_BUF;
    dom[i].Gcc._ksb = dom[i].Gcc._ks - DOM_BUF;
    dom[i].Gcc._kn = dom[i].zn;
    dom[i].Gcc._knb = dom[i].Gcc._kn + 2 * DOM_BUF;
    dom[i].Gcc._ke = dom[i].Gcc._ks + dom[i].Gcc._kn;
    dom[i].Gcc._keb = dom[i].Gcc._ke + DOM_BUF;

    dom[i].Gcc._s1 = dom[i].Gcc._in;
    dom[i].Gcc._s2 = dom[i].Gcc._s1 * dom[i].Gcc._jn;
    dom[i].Gcc._s3 = dom[i].Gcc._s2 * dom[i].Gcc._kn;
    dom[i].Gcc._s1b = dom[i].Gcc._inb;
    dom[i].Gcc._s2b = dom[i].Gcc._s1b * dom[i].Gcc._jnb;
    dom[i].Gcc._s3b = dom[i].Gcc._s2b * dom[i].Gcc._knb;

    // Gfx
    if(dom[i].W > -1)
      dom[i].Gfx.is = dom[dom[i].W].Gfx.ie - 1;
    else
      dom[i].Gfx.is = DOM_BUF;
    dom[i].Gfx.isb = dom[i].Gfx.is - DOM_BUF;
    dom[i].Gfx.in = dom[i].xn + 1;
    dom[i].Gfx.inb = dom[i].Gfx.in + 2 * DOM_BUF;
    dom[i].Gfx.ie = dom[i].Gfx.is + dom[i].Gfx.in;
    dom[i].Gfx.ieb = dom[i].Gfx.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfx.js = dom[dom[i].S].Gfx.je - 1;
    else
      dom[i].Gfx.js = DOM_BUF;
    dom[i].Gfx.jsb = dom[i].Gfx.js - DOM_BUF;
    dom[i].Gfx.jn = dom[i].yn;
    dom[i].Gfx.jnb = dom[i].Gfx.jn + 2 * DOM_BUF;
    dom[i].Gfx.je = dom[i].Gfx.js + dom[i].Gfx.jn;
    dom[i].Gfx.jeb = dom[i].Gfx.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfx.ks = dom[dom[i].B].Gfx.ke - 1;
    else
      dom[i].Gfx.ks = DOM_BUF;
    dom[i].Gfx.ksb = dom[i].Gfx.ks - DOM_BUF;
    dom[i].Gfx.kn = dom[i].zn;
    dom[i].Gfx.knb = dom[i].Gfx.kn + 2 * DOM_BUF;
    dom[i].Gfx.ke = dom[i].Gfx.ks + dom[i].Gfx.kn;
    dom[i].Gfx.keb = dom[i].Gfx.ke + DOM_BUF;

    dom[i].Gfx.s1 = dom[i].Gfx.in;
    dom[i].Gfx.s2 = dom[i].Gfx.s1 * dom[i].Gfx.jn;
    dom[i].Gfx.s3 = dom[i].Gfx.s2 * dom[i].Gfx.kn;
    dom[i].Gfx.s1b = dom[i].Gfx.inb;
    dom[i].Gfx.s2b = dom[i].Gfx.s1b * dom[i].Gfx.jnb;
    dom[i].Gfx.s3b = dom[i].Gfx.s2b * dom[i].Gfx.knb;

    dom[i].Gfx._is = DOM_BUF;
    dom[i].Gfx._isb = dom[i].Gfx._is - DOM_BUF;
    dom[i].Gfx._in = dom[i].xn + 1;
    dom[i].Gfx._inb = dom[i].Gfx._in + 2 * DOM_BUF;
    dom[i].Gfx._ie = dom[i].Gfx._is + dom[i].Gfx._in;
    dom[i].Gfx._ieb = dom[i].Gfx._ie + DOM_BUF;

    dom[i].Gfx._js = DOM_BUF;
    dom[i].Gfx._jsb = dom[i].Gfx._js - DOM_BUF;
    dom[i].Gfx._jn = dom[i].yn;
    dom[i].Gfx._jnb = dom[i].Gfx._jn + 2 * DOM_BUF;
    dom[i].Gfx._je = dom[i].Gfx._js + dom[i].Gfx._jn;
    dom[i].Gfx._jeb = dom[i].Gfx._je + DOM_BUF;

    dom[i].Gfx._ks = DOM_BUF;
    dom[i].Gfx._ksb = dom[i].Gfx._ks - DOM_BUF;
    dom[i].Gfx._kn = dom[i].zn;
    dom[i].Gfx._knb = dom[i].Gfx._kn + 2 * DOM_BUF;
    dom[i].Gfx._ke = dom[i].Gfx._ks + dom[i].Gfx._kn;
    dom[i].Gfx._keb = dom[i].Gfx._ke + DOM_BUF;

    dom[i].Gfx._s1 = dom[i].Gfx._in;
    dom[i].Gfx._s2 = dom[i].Gfx._s1 * dom[i].Gfx._jn;
    dom[i].Gfx._s3 = dom[i].Gfx._s2 * dom[i].Gfx._kn;
    dom[i].Gfx._s1b = dom[i].Gfx._inb;
    dom[i].Gfx._s2b = dom[i].Gfx._s1b * dom[i].Gfx._jnb;
    dom[i].Gfx._s3b = dom[i].Gfx._s2b * dom[i].Gfx._knb;

    // Gfy
    if(dom[i].W > -1)
      dom[i].Gfy.is = dom[dom[i].W].Gfy.ie;
    else
      dom[i].Gfy.is = DOM_BUF;
    dom[i].Gfy.isb = dom[i].Gfy.is - DOM_BUF;
    dom[i].Gfy.in = dom[i].xn;
    dom[i].Gfy.inb = dom[i].Gfy.in + 2 * DOM_BUF;
    dom[i].Gfy.ie = dom[i].Gfy.is + dom[i].Gfy.in;
    dom[i].Gfy.ieb = dom[i].Gfy.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfy.js = dom[dom[i].S].Gfy.je;
    else
      dom[i].Gfy.js = DOM_BUF;
    dom[i].Gfy.jsb = dom[i].Gfy.js - DOM_BUF;
    dom[i].Gfy.jn = dom[i].yn + 1;
    dom[i].Gfy.jnb = dom[i].Gfy.jn + 2 * DOM_BUF;
    dom[i].Gfy.je = dom[i].Gfy.js + dom[i].Gfy.jn;
    dom[i].Gfy.jeb = dom[i].Gfy.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfy.ks = dom[dom[i].B].Gfy.ke;
    else
      dom[i].Gfy.ks = DOM_BUF;
    dom[i].Gfy.ksb = dom[i].Gfy.ks - DOM_BUF;
    dom[i].Gfy.kn = dom[i].zn;
    dom[i].Gfy.knb = dom[i].Gfy.kn + 2 * DOM_BUF;
    dom[i].Gfy.ke = dom[i].Gfy.ks + dom[i].Gfy.kn;
    dom[i].Gfy.keb = dom[i].Gfy.ke + DOM_BUF;

    dom[i].Gfy.s1 = dom[i].Gfy.in;
    dom[i].Gfy.s2 = dom[i].Gfy.s1 * dom[i].Gfy.jn;
    dom[i].Gfy.s3 = dom[i].Gfy.s2 * dom[i].Gfy.kn;
    dom[i].Gfy.s1b = dom[i].Gfy.inb;
    dom[i].Gfy.s2b = dom[i].Gfy.s1b * dom[i].Gfy.jnb;
    dom[i].Gfy.s3b = dom[i].Gfy.s2b * dom[i].Gfy.knb;

    dom[i].Gfy._is = DOM_BUF;
    dom[i].Gfy._isb = dom[i].Gfy._is - DOM_BUF;
    dom[i].Gfy._in = dom[i].xn;
    dom[i].Gfy._inb = dom[i].Gfy._in + 2 * DOM_BUF;
    dom[i].Gfy._ie = dom[i].Gfy._is + dom[i].Gfy._in;
    dom[i].Gfy._ieb = dom[i].Gfy._ie + DOM_BUF;

    dom[i].Gfy._js = DOM_BUF;
    dom[i].Gfy._jsb = dom[i].Gfy._js - DOM_BUF;
    dom[i].Gfy._jn = dom[i].yn + 1;
    dom[i].Gfy._jnb = dom[i].Gfy._jn + 2 * DOM_BUF;
    dom[i].Gfy._je = dom[i].Gfy._js + dom[i].Gfy._jn;
    dom[i].Gfy._jeb = dom[i].Gfy._je + DOM_BUF;

    dom[i].Gfy._ks = DOM_BUF;
    dom[i].Gfy._ksb = dom[i].Gfy._ks - DOM_BUF;
    dom[i].Gfy._kn = dom[i].zn;
    dom[i].Gfy._knb = dom[i].Gfy._kn + 2 * DOM_BUF;
    dom[i].Gfy._ke = dom[i].Gfy._ks + dom[i].Gfy._kn;
    dom[i].Gfy._keb = dom[i].Gfy._ke + DOM_BUF;

    dom[i].Gfy._s1 = dom[i].Gfy._in;
    dom[i].Gfy._s2 = dom[i].Gfy._s1 * dom[i].Gfy._jn;
    dom[i].Gfy._s3 = dom[i].Gfy._s2 * dom[i].Gfy._kn;
    dom[i].Gfy._s1b = dom[i].Gfy._inb;
    dom[i].Gfy._s2b = dom[i].Gfy._s1b * dom[i].Gfy._jnb;
    dom[i].Gfy._s3b = dom[i].Gfy._s2b * dom[i].Gfy._knb;

    // Gfz
    if(dom[i].W > -1)
      dom[i].Gfz.is = dom[dom[i].W].Gfz.ie;
    else
      dom[i].Gfz.is = DOM_BUF;
    dom[i].Gfz.isb = dom[i].Gfz.is - DOM_BUF;
    dom[i].Gfz.in = dom[i].xn;
    dom[i].Gfz.inb = dom[i].Gfz.in + 2 * DOM_BUF;
    dom[i].Gfz.ie = dom[i].Gfz.is + dom[i].Gfz.in;
    dom[i].Gfz.ieb = dom[i].Gfz.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfz.js = dom[dom[i].S].Gfz.je;
    else
      dom[i].Gfz.js = DOM_BUF;
    dom[i].Gfz.jsb = dom[i].Gfz.js - DOM_BUF;
    dom[i].Gfz.jn = dom[i].yn;
    dom[i].Gfz.jnb = dom[i].Gfz.jn + 2 * DOM_BUF;
    dom[i].Gfz.je = dom[i].Gfz.js + dom[i].Gfz.jn;
    dom[i].Gfz.jeb = dom[i].Gfz.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfz.ks = dom[dom[i].B].Gfz.ke;
    else
      dom[i].Gfz.ks = DOM_BUF;
    dom[i].Gfz.ksb = dom[i].Gfz.ks - DOM_BUF;
    dom[i].Gfz.kn = dom[i].zn + 1;
    dom[i].Gfz.knb = dom[i].Gfz.kn + 2 * DOM_BUF;
    dom[i].Gfz.ke = dom[i].Gfz.ks + dom[i].Gfz.kn;
    dom[i].Gfz.keb = dom[i].Gfz.ke + DOM_BUF;

    dom[i].Gfz.s1 = dom[i].Gfz.in;
    dom[i].Gfz.s2 = dom[i].Gfz.s1 * dom[i].Gfz.jn;
    dom[i].Gfz.s3 = dom[i].Gfz.s2 * dom[i].Gfz.kn;
    dom[i].Gfz.s1b = dom[i].Gfz.inb;
    dom[i].Gfz.s2b = dom[i].Gfz.s1b * dom[i].Gfz.jnb;
    dom[i].Gfz.s3b = dom[i].Gfz.s2b * dom[i].Gfz.knb;

    dom[i].Gfz._is = DOM_BUF;
    dom[i].Gfz._isb = dom[i].Gfz._is - DOM_BUF;
    dom[i].Gfz._in = dom[i].xn;
    dom[i].Gfz._inb = dom[i].Gfz._in + 2 * DOM_BUF;
    dom[i].Gfz._ie = dom[i].Gfz._is + dom[i].Gfz._in;
    dom[i].Gfz._ieb = dom[i].Gfz._ie + DOM_BUF;

    dom[i].Gfz._js = DOM_BUF;
    dom[i].Gfz._jsb = dom[i].Gfz._js - DOM_BUF;
    dom[i].Gfz._jn = dom[i].yn;
    dom[i].Gfz._jnb = dom[i].Gfz._jn + 2 * DOM_BUF;
    dom[i].Gfz._je = dom[i].Gfz._js + dom[i].Gfz._jn;
    dom[i].Gfz._jeb = dom[i].Gfz._je + DOM_BUF;

    dom[i].Gfz._ks = DOM_BUF;
    dom[i].Gfz._ksb = dom[i].Gfz._ks - DOM_BUF;
    dom[i].Gfz._kn = dom[i].zn + 1;
    dom[i].Gfz._knb = dom[i].Gfz._kn + 2 * DOM_BUF;
    dom[i].Gfz._ke = dom[i].Gfz._ks + dom[i].Gfz._kn;
    dom[i].Gfz._keb = dom[i].Gfz._ke + DOM_BUF;

    dom[i].Gfz._s1 = dom[i].Gfz._in;
    dom[i].Gfz._s2 = dom[i].Gfz._s1 * dom[i].Gfz._jn;
    dom[i].Gfz._s3 = dom[i].Gfz._s2 * dom[i].Gfz._kn;
    dom[i].Gfz._s1b = dom[i].Gfz._inb;
    dom[i].Gfz._s2b = dom[i].Gfz._s1b * dom[i].Gfz._jnb;
    dom[i].Gfz._s3b = dom[i].Gfz._s2b * dom[i].Gfz._knb;
  }

  // set up grid index structs
  // allocate and initialize pressure and velocity vectors
  p0 = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  p = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  //divU = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u0 = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v0 = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w0 = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_star = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v_star = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w_star = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_WE = (real*) malloc(Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real);
  u_SN = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_BT = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  v_WE = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_SN = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real);
  v_BT = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  w_WE = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_SN = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_BT = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real);
  f_x = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  f_y = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  f_z = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
	
	//##########################################################################
	// velocity field initialization(turb or not)
	if(turb_switch == ON) {
		
		// set up the random number generator
		srand(time(NULL));

		for(i = 0; i < Dom.Gcc.s3b; i++) {
			p0[i] = 0.;
			p[i] = 0.;
			//divU[i] = 0.;
		}
		for(i = 0; i < Dom.Gfx.s3b; i++) {
			u[i] = 0.;
			diff0_u[i] = 0.;
			diff_u[i] = 0.;
			u0[i] = 0.;
			conv0_u[i] = 0.;
			conv_u[i] = 0.;
			f_x[i] = 0.;
		}
		for(i = 0; i < Dom.Gfy.s3b; i++) {
			v[i] = 0.;
			diff0_v[i] = 0.;
			diff_v[i] = 0.;
			v0[i] = 0.;
			conv0_v[i] = 0.;
			conv_v[i] = 0.;
			f_y[i] = 0.;
		}
		for(i = 0; i < Dom.Gfz.s3b; i++) {
			w[i] = 0.;
			diff0_w[i] = 0.;
			diff_w[i] = 0.;
			w0[i] = 0.;
			conv0_w[i] = 0.;
			conv_w[i] = 0.;
			f_z[i] = 0.;
		}

		// integral scale
		turbl = (Dom.xl + Dom.yl + Dom.zl) / 3.;
		real urms = 3*turbA*turbl;

		// randomly initialize velocity components
		for(i = 0; i < Dom.Gfx.s3b; i++) {
			tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
			u[i] = tmp;
			u0[i] = tmp;
		}
		for(i = 0; i < Dom.Gfy.s3b; i++) {
			tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
			v[i] = tmp;
			v0[i] = tmp;
		}
		for(i = 0; i < Dom.Gfz.s3b; i++) {
			tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
			w[i] = tmp;
			w0[i] = tmp;
		}

		// calculate the divergence of U
		real vol = (Dom.xn+2*DOM_BUF)*(Dom.yn+2*DOM_BUF)*(Dom.zn+2*DOM_BUF);
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
					W = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
					E = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
					S = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					N = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					B = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
					T = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
					C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
					p[C] = (u[E]-u[W])/Dom.dx + (v[N]-v[S])/Dom.dy + (w[T]-w[B])/Dom.dz;
					p[C] = p[C] / vol;
				}
			}
		}

		real umean = 0.;
		real vmean = 0.;
		real wmean = 0.;
		// subtract off the divergence of U to make U solenoidal
		for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
			for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
				for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
					C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
					W = (i-1) + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
					E = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
					u[C] = u[C] - 0.5*(p[W] + p[E]);
					u0[C] = u0[C] - 0.5*(p[W] + p[E]);
					umean += u[C];
				}
			}
		}
		for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
			for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
				for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
					C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					S = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					N = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					v[C] = v[C] - 0.5*(p[S] + p[N]);
					v0[C] = v0[C] - 0.5*(p[S] + p[N]);
					vmean += v[C];
				}
			}
		}
		for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
			for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
				for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
					C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
					B = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
					T = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
					w[C] = w[C] - 0.5*(p[B] + p[T]);
					w0[C] = w0[C] - 0.5*(p[B] + p[T]);
					wmean += w[C];
				}
			}
		}

		umean /= vol;
		vmean /= vol;
		wmean /= vol;

		// re-scale to give zero mean velocity in each direction
		for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
			for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
				for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
					C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
					u[C] = u[C] - umean;
					u0[C] = u0[C] - umean;
				}
			}
		}
		for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
			for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
				for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
					C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					v[C] = v[C] - vmean;
					v0[C] = v0[C] - vmean;
				}
			}
		}
		for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
			for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
				for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
					C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
					w[C] = w[C] - wmean;
					w0[C] = w0[C] - wmean;
				}
			}
		}
	} else if(turb_switch == OFF) {
		
		// initialize QUIESCENT flow (default)
		for(i = 0; i < Dom.Gcc.s3b; i++) {
			p0[i] = 0.;
			p[i] = 0.;
			//divU[i] = 0.;
		}
		
		for(i = 0; i < Dom.Gfx.s3b; i++) {
			u[i] = 0.;
			diff0_u[i] = 0.;
			diff_u[i] = 0.;
			u0[i] = 0.;
			conv0_u[i] = 0.;
			conv_u[i] = 0.;
			f_x[i] = 0.;
			u_star[i] = 0.;
		}
		for(i = 0; i < Dom.Gfy.s3b; i++) {
			v[i] = 0.;
			diff0_v[i] = 0.;
			diff_v[i] = 0.;
			v0[i] = 0.;
			conv0_v[i] = 0.;
			conv_v[i] = 0.;
			f_y[i] = 0.;
			v_star[i] = 0.;
		}
		for(i = 0; i < Dom.Gfz.s3b; i++) {
			w[i] = 0.;
			diff0_w[i] = 0.;
			diff_w[i] = 0.;
			w0[i] = 0.;
			conv0_w[i] = 0.;
			conv_w[i] = 0.;
			f_z[i] = 0.;
			w_star[i] = 0.;
		}
		
		//#######################
		/*
		for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
			for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
				for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
					real x = ((i-0.5) * Dom.dx) + Dom.xs;
					real y = ((j-0.5) * Dom.dy) + Dom.ys;
					//real r = sqrt(x*x+y*y);
					int C = i+j*Dom.Gfz.s1b+k*Dom.Gfz.s2b;
					w[C] = 1e-3*(4.0-fabs(x)-fabs(y));
					w0[C] = w[C];
				}
			}
		}
		*/
		//##########################
		
		// prevent turbulence linear forcing from being used
		turbl = 0.; // integral scale
	}
	//##########################################################################
	
  for(i = 0; i < Dom.Gfx.jnb*Dom.Gfx.knb; i++) {
    u_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.knb; i++) {
    u_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.jnb; i++) {
    u_BT[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.jnb*Dom.Gfy.knb; i++) {
    v_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.knb; i++) {
    v_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.jnb; i++) {
    v_BT[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.jnb*Dom.Gfz.knb; i++) {
    w_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.knb; i++) {
    w_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.jnb; i++) {
    w_BT[i] = 0.;
  }

  // initialize some variables
  dt = 2 * nu / (Dom.dx * Dom.dx);
  dt += 2 * nu / (Dom.dy * Dom.dy);
  dt += 2 * nu / (Dom.dz * Dom.dz);
  dt = CFL / dt;
  dt0 = -1.;
  stepnum = 0;
  rec_flow_field_stepnum_out = 0;
  rec_paraview_stepnum_out = 0;
  rec_particle_stepnum_out = 0;
  rec_precursor_stepnum_out = 0;
  rec_flow_field_ttime_out = 0;
  rec_paraview_ttime_out = 0;
  rec_particle_ttime_out = 0;
  rec_restart_ttime_out = 0;
  rec_precursor_ttime_out = 0;

  return EXIT_SUCCESS;
}

void cgns_flow_field_Eulerian2(void)
{
	// create the solution file
	char fname[FILE_NAME_SIZE];
	char fnameall[FILE_NAME_SIZE];
	char gname[FILE_NAME_SIZE];
	char gnameall[FILE_NAME_SIZE];
	char snodename[CHAR_BUF_SIZE];
	
	// flow file name
	sprintf(fname, "flow-%f.cgns", ttime);
	// flow file name and path
	sprintf(fnameall, "%s/output/%s", ROOT_DIR, fname);
	
	//sprintf(fnameall2, "%s/output/flow-%s.cgns", ROOT_DIR, format);
	sprintf(snodename, "Solution-");
	//sprintf(snodenameall, "/Base/Zone0/Solution-");
	sprintf(snodename, "%s%f", snodename, ttime);
	//sprintf(snodenameall, "%s%s", snodenameall, format);
	//sprintf(fname, fname2, tout);
	//sprintf(fnameall, fnameall2, tout);
	//sprintf(snodename, snodename, ttime);
	//sprintf(snodenameall, snodenameall, tout);

	// grid file name
	sprintf(gname, "grid.cgns");
	// grid file name and path
	sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
	
	int fn;
	// in our case we have only one base and one zone, need to check in future
	int bn;
	int zn;
	int sn;
	int fnpress;
	int fnu;
	int fnv;
	int fnw;
	//int fnwb;
	int fnnumden;
	//int fnbubmas;
	int fnconcen;
	int fnbubdia;

	// open the output file
	cg_open(fnameall, CG_MODE_WRITE, &fn);
	
	cg_base_write(fn, "Base", 3, 3, &bn);
	cgsize_t size[9];
	size[0] = Dom.xn+1; // cells -> vertices
	size[1] = Dom.yn+1;
	size[2] = Dom.zn+1;
	size[3] = Dom.xn;
	size[4] = Dom.yn;
	size[5] = Dom.zn;
	size[6] = 0;
	size[7] = 0;
	size[8] = 0;
	cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
	cg_goto(fn, bn, "Zone_t", zn, "end");
	cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");
	
	// write the solution node for this time step
	cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
	
	// write the fields
	//==========================================================================
	real *pout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				pout[C] = p[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Pressure", pout, &fnpress);
	
	real *uout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
		for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
			for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = (i-1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC1 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC2 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				int CC3 = (i+2) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
				uout[C] = -0.0625*u[CC0] + 0.5625*u[CC1] + 0.5625*u[CC2] - 0.0625*u[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", uout, &fnu);
	
	real *vout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
		for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
			for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC1 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC2 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				int CC3 = i + (j+2)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
				vout[C] = -0.0625*v[CC0] + 0.5625*v[CC1] + 0.5625*v[CC2] - 0.0625*v[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", vout, &fnv);
	
	real *wout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
		for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
			for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC0 = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
				int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
				int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
				int CC3 = i + j*Dom.Gfz.s1b + (k+2)*Dom.Gfz.s2b;
				wout[C] = -0.0625*w[CC0] + 0.5625*w[CC1] + 0.5625*w[CC2] - 0.0625*w[CC3];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", wout, &fnw);
	/*
	real *wbout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
		for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
			for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
				int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
				wbout[C] = 0.5 * (w_b[CC1] + w_b[CC2]);
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "BubbleVelocityZ", wbout, &fnwb);
	*/
	real *numdenout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				numdenout[C] = numden[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "NumberDensity", numdenout, &fnnumden);
	/*
	real *bubmasout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				bubmasout[C] = bubmas[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "BubbleMass", bubmasout, &fnbubmas);
	*/
	real *concenout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				concenout[C] = concen[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Concentration", concenout, &fnconcen);
	
	real *bubdiaout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				bubdiaout[C] = bubdia[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "Bubble Diameter", bubdiaout, &fnbubdia);
	
	free(pout);
	free(uout);
	free(vout);
	free(wout);
	//free(wbout);
	free(numdenout);
	//free(bubmasout);
	free(concenout);
	free(bubdiaout);
	
	cg_close(fn);
	//==========================================================================
	
	// now add this timestep into time series
	int N = 0;
	char buf[32];
	cgsize_t tsize[1];
	cgsize_t solsize[2];
	real *cgns_tseries;
	char *BaseIterName = malloc(CHAR_BUF_SIZE * sizeof(real));
	
	// open grid file
	int gn;
	bn = 1;
	zn = 1;
	cg_open(gnameall, CG_MODE_MODIFY, &gn);
	
	// get the number of solutions
	cg_goto(gn, bn, "Zone_t", zn, "end");
	cg_link_write(snodename, fname, "Base/Zone0/Solution");
	cg_biter_read(gn, bn, BaseIterName, &N);
	
	cgns_tseries = malloc((N + 1) * sizeof(real));
	cg_goto(gn, bn, "BaseIterativeData_t", 1, "end");
	cg_array_read(1, cgns_tseries);
	cg_biter_write(gn, bn, "BaseIterativeData", N + 1);
	cg_goto(gn, bn, "BaseIterativeData_t", 1, "end");
	cgns_tseries[N] = ttime;
	tsize[0] = N + 1;
	cg_array_write("TimeValues", RealDouble, 1, tsize, cgns_tseries);
	
	const char *solname[N + 1];
	GridLocation_t location;
	for(int i = 0; i < N; i++) {
		cg_sol_info(gn, bn, zn, i + 1, buf, &location);
		//printf("===%s===\n", buf);
		solname[i] = buf;
	}
	solname[N] = snodename;
	cg_goto(gn, bn, "Zone_t", zn, "ZoneIterativeData_t", 1, "end");
	solsize[0] = 32;
	solsize[1] = N + 1;
	cg_array_write("FlowSolutionPointers", Character, 2, solsize, solname);
	
	free(cgns_tseries);
	cg_close(gn);
	free(BaseIterName);
}
