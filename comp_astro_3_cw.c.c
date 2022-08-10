#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N 100
#define gamma 1.4
#define x_i 0.0
#define x_f 1.0

double delta_x = (x_f - x_i)/N;
double time = 0.0;
double v[N+2], p[N+2], e[N+2], c_s[N+2], a_pos[N+2];
double x[N+2], E[N+2], rho[N+2], f1, f2, f3, a[N+2];
double delta_t, CFL, a_max, snapshot_time;

////////////////////////////////////////////////////////////////////////////////
/// VECTOR STUFF

// making a struct to define a vector
typedef struct vector{
   double a;
   double b;
   double c;
}vector;

//vector addition/subtraction
//factor is if it is addition or subtraction
vector vect_pm(vector vect_1, double factor, vector vect_2){
   vector add;
   add.a = vect_1.a + factor*vect_2.a;
   add.b = vect_1.b + factor*vect_2.b;
   add.c = vect_1.c + factor*vect_2.c;
   return add;
}

////////////////////////////////////////////////////////////////////////////////
/// EQUATIONS ///
double pressure_calc(double TIE, double rho, double v){
   return (gamma - 1.0)*(TIE - (0.5*rho*pow(v, 2)));
}
double sound_speed(double p, double rho){
   double temp = (gamma*p)/rho;
   return sqrt(temp);
}
double specific_IE(double p, double rho){return p/((gamma - 1.0)*rho);}
////////////////////////////////////////////////////////////////////////////////

//convert state vector into flux vector
vector state_to_flux(vector state_vector){
   vector flux;
   double pres = pressure_calc(state_vector.c, state_vector.a, (state_vector.b/state_vector.a));
   flux.a = state_vector.b;
   flux.b = (pow(state_vector.b, 2)/state_vector.a) + pres;
   flux.c = state_vector.b*((state_vector.c + pres)/state_vector.a);
   return flux;
}
//finds the maximum value in an array
double max(double arr[], int n){
   double maximum = fabs(arr[0]);
   for(int i = 1; i < n; i++){
      if(fabs(arr[i]) > maximum){
         maximum = (arr[i]);
      }
   }
   return maximum;
}

////////////////////////////////////////////////////////////////////////////////

void set_initial_tube(int IC, vector q[], vector f[], FILE* ptr){
   if(IC == 1){
      snapshot_time = 0.2;
      for(int i = 0; i <= N+1; i++){
         x[i] = (double)i/N;
         if(x[i] <= 0.3){rho[i] = 1.0; v[i] = 0.75; p[i] = 1.0;}
         else{rho[i] = 0.125; v[i] = 0.0; p[i] = 0.1;}
      }
   }
   if(IC == 2){
      snapshot_time = 0.012;
      for(int i = 0; i <= N+1; i++){
         x[i] = (double)i/N;
         if(x[i] <= 0.8){rho[i] = 1.0; v[i] = -19.59745; p[i] = 1000.;}
         else{rho[i] = 1.0; v[i] = -19.59745; p[i] = 0.01;}
      }
   }
   for(int i = 0; i<=N+1; i++){
   // total energy density
      E[i] = (p[i]/((gamma - 1.))) + 0.5*rho[i]*v[i]*v[i];
      q[i].a = rho[i];
      q[i].b = rho[i]*v[i];
      q[i].c = E[i];

      f[i] = state_to_flux(q[i]);

      c_s[i] = sound_speed(p[i], rho[i]);
      e[i] = specific_IE(p[i], rho[i]);
      a[i] = fabs(v[i]) + c_s[i];
      fprintf(ptr, "%.3f, %.3f, %.3f, %.3f, %.3f\n", x[i], rho[i], v[i], p[i], e[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////


vector lax_friedrichs(vector f0, vector fp, vector fm, vector q0, vector qp, vector qm){
   /// determine fluxes within cells (ph = plus half, mh = minus half) ///
   // f0 = f[i], fp = f[i+1], fm = f[i-1] (same with q's)
   vector f_ph, f_mh;

   f_ph.a = 0.5*(f0.a + fp.a) + (delta_x/(2.*delta_t))*(q0.a - qp.a);
   f_ph.b = 0.5*(f0.b + fp.b) + (delta_x/(2.*delta_t))*(q0.b - qp.b);
   f_ph.c = 0.5*(f0.c + fp.c) + (delta_x/(2.*delta_t))*(q0.c - qp.c);

   f_mh.a = 0.5*(fm.a + f0.a) + (delta_x/(2.*delta_t))*(qm.a - q0.a);
   f_mh.b = 0.5*(fm.b + f0.b) + (delta_x/(2.*delta_t))*(qm.b - q0.b);
   f_mh.c = 0.5*(fm.c + f0.c) + (delta_x/(2.*delta_t))*(qm.c - q0.c);

   return vect_pm(f_ph, -1, f_mh);
}

vector lax_wendroff(vector q0, vector qp, vector qm, vector f0, vector fp, vector fm){
   /// determine fluxes within cells (ph = plus half, mh = minus half) ///
   vector q_ph, q_mh, f_ph, f_mh;

   q_ph.a = 0.5*(q0.a + qp.a) - (delta_t/(2.*delta_x))*(-f0.a + fp.a);
   q_ph.b = 0.5*(q0.b + qp.b) - (delta_t/(2.*delta_x))*(-f0.b + fp.b);
   q_ph.c = 0.5*(q0.c + qp.c) - (delta_t/(2.*delta_x))*(-f0.c + fp.c);

   q_mh.a = 0.5*(qm.a + q0.a) - (delta_t/(2.*delta_x))*(-fm.a + f0.a);
   q_mh.b = 0.5*(qm.b + q0.b) - (delta_t/(2.*delta_x))*(-fm.b + f0.b);
   q_mh.c = 0.5*(qm.c + q0.c) - (delta_t/(2.*delta_x))*(-fm.c + f0.c);

   f_ph = state_to_flux(q_ph);
   f_mh = state_to_flux(q_mh);

   return vect_pm(f_ph, -1, f_mh);
}


////////////////////////////////////////////////////////////////////////////////


int main(){
	FILE *fp; FILE *fp_orig;
	fp = fopen("out.txt", "w"); fp_orig = fopen("orig.txt", "w");
	if(!fp || !fp_orig){return 1;}

   //allocate memory
	vector *q = malloc(sizeof *q + 3*(N+2)*sizeof(double));
   vector *q_new = malloc(sizeof *q_new + 3*(N+2)*sizeof(double));
   vector *f = malloc(sizeof *f + 3*(N+2)*sizeof(double));
	vector dF;

   // 1 for fig 1, 2 for fig 2
	set_initial_tube(2, q, f, fp_orig);


   int loop_count = 0;
	while(time <= snapshot_time)
   {
      loop_count += 1;
      printf("%i\n", loop_count);

      ///////////////////////////
      /// BOUNDARY CONDITIONS ///
      ///////////////////////////
      q[0].a = q[1].a; q[0].b = q[1].b; q[0].c = q[1].c;
      q[N+1].a = q[N].a; q[N+1].b = q[N].b; q[N+1].c = q[N].c;

      /// redefine flux vector
      for(int i=0; i<=N+1; i++){
         f[i] = state_to_flux(q[i]);
      }

      ///////////////////////
      /// UPDATE TIMESTEP ///
      ///////////////////////
      double a_max = max(a, N);
      delta_t = 0.33*(delta_x/a_max);


      for(int i = 1; i < N+1; i++){
         dF = lax_friedrichs(f[i], f[i+1], f[i-1], q[i], q[i+1], q[i-1]);
         //dF = lax_wendroff(q[i], q[i+1], q[i-1], f[i], f[i+1], f[i-1]);

         q_new[i].a = q[i].a - (delta_t/delta_x)*(dF.a);
         q_new[i].b = q[i].b - (delta_t/delta_x)*(dF.b);
         q_new[i].c = q[i].c - (delta_t/delta_x)*(dF.c);


         rho[i] = q_new[i].a;
         v[i] = q_new[i].b/rho[i];
         E[i] = q_new[i].c;

         p[i] = pressure_calc(E[i], rho[i], v[i]);
         c_s[i] = sound_speed(p[i], rho[i]);
         a[i] = fabs(v[i]) + c_s[i];
         e[i] = specific_IE(p[i], rho[i]);
      }

      // copy calculated state vector into q array of struct
      for(int i = 1; i < N+1; i++){
         q[i].a = q_new[i].a;
         q[i].b = q_new[i].b;
         q[i].c = q_new[i].c;
      }

		///////////////////
		/// UPDATE TIME ///
		///////////////////
      time += delta_t;
	}

	//data dump
	for(int i = 1; i < N+1; i++){
      fprintf(fp, "%.3f, %.3f, %.3f, %.3f, %.3f\n", x[i], q[i].a, v[i], p[i], e[i]);
	}

	fclose(fp); fclose(fp_orig);
	free(q); free(q_new); free(f);

	return 0;
}

