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
	real h; // depth

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
	bubden = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
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
	
	// initialize cell-centered and face-centered bubble density field, 
	// need to define the gravity direction, note that POSITIVE direction have  
	// to be the water surface and NEGATIVE direction have to be the bottom.
	if(gravity_direction == GRAVITY_X) {
		
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
	}
	else if(gravity_direction == GRAVITY_Y) {
		
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
	}
	else if(gravity_direction == GRAVITY_Z) {
		
		bubden_face = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
		cpumem += Dom.Gfz.s3b * sizeof(real);
		
		for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb - 1; k++) {
			for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
				for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
					C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
					h = ((real)(Dom.Gcc.keb - k) - 1.5) * Dom.dz;
					bubden[C] = rho_atm * (1.0 + rho_f * grav_acc * h / pressure_atm);
				}
			}
		}
		for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
			for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
				k = Dom.Gcc.keb - 1;
				C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
				bubden[C] = rho_atm;
			}
		}
		for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb - 1; k++) {
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
				C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
				bubden_face[C] = rho_atm;
			}
		}
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
	printf("gravity_direction ");
	printf("%d\n", gravity_direction);
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
  cg_biter_write(fn, bn, "BaseIterativeData", 1);
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
  
  
  
  // create solution name recorder file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, "solname.rec");
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", "solname.rec");
    exit(EXIT_FAILURE);
  }
  fprintf(rec, "Solution names:\n");
  fclose(rec);
}

void cgns_flow_field_Eulerian(real dtout)
{
	// create the solution file
	char fname[FILE_NAME_SIZE];
	char fname2[FILE_NAME_SIZE];
	char fnameall[FILE_NAME_SIZE];
	char fnameall2[FILE_NAME_SIZE];
	char gname[FILE_NAME_SIZE];
	char gnameall[FILE_NAME_SIZE];
	real tout = ttime; //  = rec_flow_field_stepnum_out * dtout;
	char format[CHAR_BUF_SIZE];
	char snodename[CHAR_BUF_SIZE];
	char snodenameall[CHAR_BUF_SIZE];
	int sigfigs = ceil(log10(1. / dtout));
	if(sigfigs < 1) sigfigs = 1;

	sprintf(format, "%%.%df", sigfigs);

	sprintf(fname2, "flow-%s.cgns", format);
	sprintf(fnameall2, "%s/output/flow-%s.cgns", ROOT_DIR, format);
	sprintf(snodename, "Solution-");
	sprintf(snodenameall, "/Base/Zone0/Solution-");
	sprintf(snodename, "%s%s", snodename, format);
	sprintf(snodenameall, "%s%s", snodenameall, format);
	sprintf(fname, fname2, tout);
	sprintf(fnameall, fnameall2, tout);
	sprintf(snodename, snodename, tout);
	sprintf(snodenameall, snodenameall, tout);

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
	int fnnumden;
	int fnbubmas;
	int fnconcen;

	// check that grid.cgns exists
	if(cg_open(gnameall, CG_MODE_MODIFY, &fn) != 0) {
		fprintf(stderr, "CGNS flow field write failure: no grid.cgns\n");
		exit(EXIT_FAILURE);
	}
	
	cg_sol_write(fn, bn, zn, snodename, CellCenter, &sn);
	
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
	//==========================================================================
	
	free(pout);
	free(uout);
	free(vout);
	free(wout);
	free(numdenout);
	free(bubmasout);
	free(concenout);
	
	cg_close(fn);
	
	// create solution name recorder file
	char path[FILE_NAME_SIZE];
	sprintf(path, "%s/record/%s", ROOT_DIR, "solname.rec");
	FILE *rec = fopen(path, "a");
	fprintf(rec, "%s\n", snodename);
	fprintf(rec, "%lf\n", ttime);
	fclose(rec);
}

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
