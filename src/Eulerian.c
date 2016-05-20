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
		} else if(strcmp(buf, "HYPERBOLICTAN_RANDOM") == 0) {
			BubbleGenerator.BubGen_type = HYPERBOLICTAN_RANDOM;
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
		} else if(strcmp(buf, "HYPERBOLICTAN_RANDOM") == 0) {
			BubbleGenerator.BubGen_type = HYPERBOLICTAN_RANDOM;
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
	
	fret = fscanf(infile, "\n");//==============================================
	
	fret = fscanf(infile, "SIMULATION CONTROL\n");
	fret = fscanf(infile, "quiescent_fluid %s\n", buf);
	if(strcmp(buf, "ON") == 0) {
		quiescent_fluid = ON;
	} else if(strcmp(buf, "OFF") == 0) {
		quiescent_fluid = OFF;
	} else {
		fprintf(stderr, "Eulerian.config read error.\n");
		exit(EXIT_FAILURE);
	}
	
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
	ter_cell = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
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
	
	// cell-centered terminal velocity field
	for(k = 0; k < Dom.Gcc.s3b; k++) {
		ter_cell[k] = 0.;
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
	} else if(BubbleGenerator.BubGen_type == HYPERBOLICTAN_RANDOM) {
		printf("HYPERBOLICTAN_RANDOM");
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
	printf("SIMULATION CONTROL\n");
	printf("quiescent_fluid %d\n", quiescent_fluid);
	
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
	// saturation concentration at different depth
	concen_sat = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	// bubble generator(source term in number density equation)
	BGndot = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	BGmdot = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
	cpumem += Dom.Gcc.s3b * sizeof(real);
	
	// initialize cell-centered and face-centered bubble density field, 
	// need to define the gravity direction, note that POSITIVE direction have 
	// to be the water surface and NEGATIVE direction have to be the bottom.
	if(gravity_direction == GRAVITY_X) {
		
	}
	else if(gravity_direction == GRAVITY_Y) {
		
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
					BGmdot[C] = PI / 6. * pow(BubbleGenerator.BubGen_dia, 3.0) * bubden[C] * BGndot[C];
				}
			}
		}
	} else if(BubbleGenerator.BubGen_type == HYPERBOLICTAN_RANDOM) {
		
		// HYPERBOLICTAN_RANDOM bubble generator used, source terms will be calculated at each time step
		
		// set zeros
		for(C = 0; C < Dom.Gcc.s3b; C++) {
			BGndot[C] = 0.0;
			BGmdot[C] = 0.0;
		}
		
		// initialize arrays
		Ix = (real *)malloc(sizeof(real) * Dom.xn);
		Iy = (real *)malloc(sizeof(real) * Dom.yn);
		Iz = (real *)malloc(sizeof(real) * Dom.zn);
		memset(Ix, 0., Dom.xn);
		memset(Iy, 0., Dom.yn);
		memset(Iz, 0., Dom.zn);
		
		// calculate the inegration
		real xx = Dom.xs + 0.5 * Dom.dx;
		real yy = Dom.ys + 0.5 * Dom.dy;
		real zz = Dom.zs + 0.5 * Dom.dz;
		real XX = 0.5*(1+tanh((xx-BubbleGenerator.BubGen_Lx1)/BubbleGenerator.BubGen_epsx));
		XX = XX * 0.5*(1+tanh((BubbleGenerator.BubGen_Lx2-xx)/BubbleGenerator.BubGen_epsx));
		real YY = 0.5*(1+tanh((yy-BubbleGenerator.BubGen_Ly1)/BubbleGenerator.BubGen_epsy));
		YY = YY * 0.5*(1+tanh((BubbleGenerator.BubGen_Ly2-yy)/BubbleGenerator.BubGen_epsy));
		real ZZ = - pow((zz - BubbleGenerator.BubGen_z0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaZ, 2.0);
		ZZ = exp(ZZ);
		Ix[0] = XX * Dom.dx;
		Iy[0] = YY * Dom.dy;
		Iz[0] = ZZ * Dom.dz;
		
		for(i = 1; i < Dom.xn; i++)	{
			xx = Dom.xs + Dom.dx * (i + 0.5);
			XX =      0.5*(1+tanh((xx-BubbleGenerator.BubGen_Lx1)/BubbleGenerator.BubGen_epsx));
			XX = XX * 0.5*(1+tanh((BubbleGenerator.BubGen_Lx2-xx)/BubbleGenerator.BubGen_epsx));
			Ix[i] = Ix[i - 1] + XX * Dom.dx;
		}
		for(j = 1; j < Dom.yn; j++)	{
			yy = Dom.ys + Dom.dy * (j + 0.5);
			YY =      0.5*(1+tanh((yy-BubbleGenerator.BubGen_Ly1)/BubbleGenerator.BubGen_epsy));
			YY = YY * 0.5*(1+tanh((BubbleGenerator.BubGen_Ly2-yy)/BubbleGenerator.BubGen_epsy));
			Iy[j] = Iy[j - 1] + YY * Dom.dy;
		}
		for(k = 1; k < Dom.zn; k++)	{
			zz = Dom.zs + Dom.dz * (k + 0.5);
			ZZ = - pow((zz - BubbleGenerator.BubGen_z0), 2.0) / 2. / pow(BubbleGenerator.BubGen_sigmaZ, 2.0);
			ZZ = exp(ZZ);
			Iz[k] = Iz[k - 1] + ZZ * Dom.dz;
		}
		
		// normalization
		for(i = 0; i < Dom.xn; i++) {
			Ix[i] = Ix[i] / Ix[Dom.xn - 1];
		}
		for(j = 0; j < Dom.yn; j++) {
			Iy[j] = Iy[j] / Iy[Dom.yn - 1];
		}
		for(k = 0; k < Dom.zn; k++) {
			Iz[k] = Iz[k] / Iz[Dom.zn - 1];
		}
		
		// initialize the random number seed
		srand(1);
	}
	else {
		fprintf(stderr, "Bubble generator type incorrect.\n");
		exit(EXIT_FAILURE);
	}
	
	return EXIT_SUCCESS;
}

void turbulence_init_Eulerian(void)
{
	int i, j, k;    // iterator
	int C, W, E, S, N, B, T;
	real tmp;
	
	// velocity field initialization(turb or not)
	if(turb_switch == ON) {
		
		// set up the random number generator
		srand(time(NULL));
		
		for(i = 0; i < Dom.Gcc.s3b; i++) {
			p0[i] = 0.;
			p[i] = 0.;
			phi[i] = 0.;
			//divU[i] = 0.;
		}
		
		for(i = 0; i < Dom.Gfx.s3b; i++) {
			u[i] = 0.;
			#ifndef IMPLICIT
				diff0_u[i] = 0.;
			#endif
			diff_u[i] = 0.;
			u0[i] = 0.;
			conv0_u[i] = 0.;
			conv_u[i] = 0.;
			f_x[i] = 0.;
			u_star[i] = 0.;
		}
		for(i = 0; i < Dom.Gfy.s3b; i++) {
			v[i] = 0.;
			#ifndef IMPLICIT
				diff0_v[i] = 0.;
			#endif
			diff_v[i] = 0.;
			v0[i] = 0.;
			conv0_v[i] = 0.;
			conv_v[i] = 0.;
			f_y[i] = 0.;
			v_star[i] = 0.;
		}
		for(i = 0; i < Dom.Gfz.s3b; i++) {
			w[i] = 0.;
			#ifndef IMPLICIT
				diff0_w[i] = 0.;
			#endif
			diff_w[i] = 0.;
			w0[i] = 0.;
			conv0_w[i] = 0.;
			conv_w[i] = 0.;
			f_z[i] = 0.;
			w_star[i] = 0.;
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
		real vol = Dom.xn*Dom.yn*Dom.zn;
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
				}
			}
		}
		
		for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
			for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
				for(i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
					C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
					umean += u[C];
				}
			}
		}
		for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
			for(j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
				for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
					C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
					vmean += v[C];
				}
			}
		}
		for(k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
			for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
				for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
					C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
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
	}
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
	int fntervel;

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
	cg_field_write(fn, bn, zn, sn, RealDouble, "BubbleDiameter", bubdiaout, &fnbubdia);
	
	real *tervelout = malloc(Dom.Gcc.s3 * sizeof(real));
	// cpumem += Dom.Gcc.s3 * sizeof(real);
	for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
		for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
			for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
				int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
				int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
				tervelout[C] = ter_cell[CC];
			}
		}
	}
	cg_field_write(fn, bn, zn, sn, RealDouble, "TerminalVelocity", tervelout, &fntervel);
	
	free(pout);
	free(uout);
	free(vout);
	free(wout);
	free(numdenout);
	free(concenout);
	free(bubdiaout);
	free(tervelout);
	
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

void Eulerian_out_restart(void)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/input/restart.config", ROOT_DIR);
  FILE *rest = fopen(path, "w");
  if(rest == NULL) {
    fprintf(stderr, "Could not open file restart.input.\n");
    exit(EXIT_FAILURE);
  }

  fwrite(&ttime, sizeof(real), 1, rest);
  fwrite(&dt0, sizeof(real), 1, rest);
  fwrite(&dt, sizeof(real), 1, rest);
  fwrite(&stepnum, sizeof(int), 1, rest);
  fwrite(&rec_paraview_stepnum_out, sizeof(int), 1, rest);

  fwrite(u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u0, sizeof(real), Dom.Gfx.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_u, sizeof(real), Dom.Gfx.s3b, rest);
#endif
  fwrite(conv0_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(diff_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(conv_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u_star, sizeof(real), Dom.Gfx.s3b, rest);

  fwrite(v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v0, sizeof(real), Dom.Gfy.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_v, sizeof(real), Dom.Gfy.s3b, rest);
#endif
  fwrite(conv0_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(diff_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(conv_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v_star, sizeof(real), Dom.Gfy.s3b, rest);

  fwrite(w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w0, sizeof(real), Dom.Gfz.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_w, sizeof(real), Dom.Gfz.s3b, rest);
#endif
  fwrite(conv0_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(diff_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(conv_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w_star, sizeof(real), Dom.Gfz.s3b, rest);

  fwrite(p, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phi, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(p0, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phase, sizeof(int), Dom.Gcc.s3b, rest);
  fwrite(phase_shell, sizeof(int), Dom.Gcc.s3b, rest);
  //############################################################################
  fwrite(numden, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(bubmas, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(concen, sizeof(real), Dom.Gcc.s3b, rest);
  //############################################################################

  fwrite(flag_u, sizeof(int), Dom.Gfx.s3b, rest);
  fwrite(flag_v, sizeof(int), Dom.Gfy.s3b, rest);
  fwrite(flag_w, sizeof(int), Dom.Gfz.s3b, rest);

  fwrite(bc_plane_pos, sizeof(real), 9, rest);

  fwrite(parts, sizeof(part_struct), nparts, rest);

  fwrite(&coeff_stride, sizeof(int), 1, rest);

  fwrite(pnm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im00, sizeof(real), nparts*coeff_stride, rest);

  fwrite(&rec_flow_field_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_paraview_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_particle_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_restart_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_prec_ttime_out, sizeof(real), 1, rest);
  fwrite(&pid_int, sizeof(real), 1, rest);
  fwrite(&pid_back, sizeof(real), 1, rest);
  fwrite(&gradP.z, sizeof(real), 2, rest);

  // close the file
  fclose(rest);
}

void out_restart_Eulerian(void)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/input/restart.config", ROOT_DIR);
  FILE *rest = fopen(path, "w");
  if(rest == NULL) {
    fprintf(stderr, "Could not open file restart.input.\n");
    exit(EXIT_FAILURE);
  }

  fwrite(&ttime, sizeof(real), 1, rest);
  fwrite(&dt0, sizeof(real), 1, rest);
  fwrite(&dt, sizeof(real), 1, rest);
  fwrite(&stepnum, sizeof(int), 1, rest);
  fwrite(&rec_paraview_stepnum_out, sizeof(int), 1, rest);

  fwrite(u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u0, sizeof(real), Dom.Gfx.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_u, sizeof(real), Dom.Gfx.s3b, rest);
#endif
  fwrite(conv0_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(diff_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(conv_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u_star, sizeof(real), Dom.Gfx.s3b, rest);

  fwrite(v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v0, sizeof(real), Dom.Gfy.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_v, sizeof(real), Dom.Gfy.s3b, rest);
#endif
  fwrite(conv0_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(diff_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(conv_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v_star, sizeof(real), Dom.Gfy.s3b, rest);

  fwrite(w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w0, sizeof(real), Dom.Gfz.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_w, sizeof(real), Dom.Gfz.s3b, rest);
#endif
  fwrite(conv0_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(diff_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(conv_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w_star, sizeof(real), Dom.Gfz.s3b, rest);

  fwrite(p, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phi, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(p0, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phase, sizeof(int), Dom.Gcc.s3b, rest);
  fwrite(phase_shell, sizeof(int), Dom.Gcc.s3b, rest);
  //############################################################################
  fwrite(numden, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(bubmas, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(concen, sizeof(real), Dom.Gcc.s3b, rest);
  //############################################################################

  fwrite(flag_u, sizeof(int), Dom.Gfx.s3b, rest);
  fwrite(flag_v, sizeof(int), Dom.Gfy.s3b, rest);
  fwrite(flag_w, sizeof(int), Dom.Gfz.s3b, rest);

  fwrite(bc_plane_pos, sizeof(real), 9, rest);

  fwrite(parts, sizeof(part_struct), nparts, rest);

  fwrite(&coeff_stride, sizeof(int), 1, rest);

  fwrite(pnm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im00, sizeof(real), nparts*coeff_stride, rest);

  fwrite(&rec_flow_field_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_paraview_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_particle_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_restart_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_prec_ttime_out, sizeof(real), 1, rest);
  fwrite(&pid_int, sizeof(real), 1, rest);
  fwrite(&pid_back, sizeof(real), 1, rest);
  fwrite(&gradP.z, sizeof(real), 2, rest);

  // close the file
  fclose(rest);
}

void in_restart_Eulerian(void)
{
  int fret = 0;
  fret = fret; // prevent compiler warning
  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/restart.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  fret = fread(&ttime, sizeof(real), 1, infile);
  fret = fread(&dt0, sizeof(real), 1, infile);
  fret = fread(&dt, sizeof(real), 1, infile);
  fret = fread(&stepnum, sizeof(int), 1, infile);
  fret = fread(&rec_paraview_stepnum_out, sizeof(int), 1, infile);

  fret = fread(u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u0, sizeof(real), Dom.Gfx.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_u, sizeof(real), Dom.Gfx.s3b, infile);
#endif
  fret = fread(conv0_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(diff_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(conv_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u_star, sizeof(real), Dom.Gfx.s3b, infile);

  fret = fread(v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v0, sizeof(real), Dom.Gfy.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_v, sizeof(real), Dom.Gfy.s3b, infile);
#endif
  fret = fread(conv0_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(diff_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(conv_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v_star, sizeof(real), Dom.Gfy.s3b, infile);

  fret = fread(w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w0, sizeof(real), Dom.Gfz.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_w, sizeof(real), Dom.Gfz.s3b, infile);
#endif
  fret = fread(conv0_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(diff_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(conv_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w_star, sizeof(real), Dom.Gfz.s3b, infile);

  fret = fread(p, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phi, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(p0, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phase, sizeof(int), Dom.Gcc.s3b, infile);
  fret = fread(phase_shell, sizeof(int), Dom.Gcc.s3b, infile);
  //############################################################################
  fret = fread(numden, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(bubmas, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(concen, sizeof(real), Dom.Gcc.s3b, infile);
  //############################################################################

  fret = fread(flag_u, sizeof(int), Dom.Gfx.s3b, infile);
  fret = fread(flag_v, sizeof(int), Dom.Gfy.s3b, infile);
  fret = fread(flag_w, sizeof(int), Dom.Gfz.s3b, infile);

  fret = fread(bc_plane_pos, sizeof(real), 9, infile);

  fret = fread(parts, sizeof(part_struct), nparts, infile);

  fret = fread(&coeff_stride, sizeof(int), 1, infile);

  fret = fread(pnm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im00, sizeof(real), nparts*coeff_stride, infile);

  fret = fread(&rec_flow_field_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_paraview_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_particle_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_restart_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_prec_ttime_out, sizeof(real), 1, infile);
  fret = fread(&pid_int, sizeof(real), 1, infile);
  fret = fread(&pid_back, sizeof(real), 1, infile);
  fret = fread(&gradP.z, sizeof(real), 1, infile);

  // close file
  fclose(infile);
}

void compute_total_number(void)
{
	int i, j, k, C;
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
	totalnumden = totalnumden * Dom.dx * Dom.dy * Dom.dz;
	printf("Total bubble number(on host) = %f\n\n", totalnumden);
}
