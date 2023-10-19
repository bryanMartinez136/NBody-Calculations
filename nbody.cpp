#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>

#define N 9999     // number of bodies
#define MASS 0     // row in array for mass
#define X_POS 1    // row in array for x position
#define Y_POS 2    // row in array for y position
#define Z_POS 3    // row in array for z position
#define X_VEL 4    // row in array for x velocity
#define Y_VEL 5    // row in array for y velocity
#define Z_VEL 6    // row in array for z velocity
#define G 200      // "gravitational constant" (not really)
#define MU 0.001   // "frictional coefficient"
#define BOXL 100.0 // periodic boundary box length

float dt = 0.05; // time interval

float body[10000][7]; // data array of bodies

void crossProduct(float vect_A[], float vect_B[], float cross_P[]) { 
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}

void norm(float &x, float &y, float &z) {
  float mag = sqrt(x*x+y*y+z*z);
  x/=mag; y/=mag; z/=mag;
}

int main(int argc, char **argv) {

  int tmax = 0;
  float Fx_dir[N];
  float Fy_dir[N];
  float Fz_dir[N];

  if (argc != 2) {
    fprintf(stderr, "Format: %s { number of timesteps }\n", argv[0]);
    exit (-1);
  }

  tmax = atoi(argv[1]);

  // assign each body a random initial positions and velocities
  srand48(time(NULL));
  float vect_A[3];
  float vect_B[3];
  float cross_P[3];

  // black hole at the center
  body[0][MASS] = 4000.0;
  body[0][X_POS] = 0.0;
  body[0][Y_POS] = 0.0;
  body[0][Z_POS] = 0.0;
  body[0][X_VEL] = 0.0;
  body[0][Y_VEL] = 0.0;
  body[0][Z_VEL] = 0.0;
  
  for (int i = 1; i < N; i++) {
    body[i][MASS] = 0.001;

    // TODO: initial coordinates centered on origin, ranging -150.0 to +150.0
    body[i][X_POS] = drand48()*300 -150;
    body[i][Y_POS] = drand48()*300 -150;
    body[i][Z_POS] = drand48()*300 -150;

    // initial velocities directions around z-axis
    vect_A[0]= body[i][X_POS];
    vect_A[1]= body[i][Y_POS];
    vect_A[2]= body[i][Z_POS];
    norm(vect_A[0], vect_A[1], vect_A[2]);
    vect_B[0]= 0.0; vect_B[1]= 0.0; vect_B[2]= 1.0;
    cross_P[0] = 0.0; cross_P[1] = 0.0; cross_P[2] = 0.0; 
    crossProduct(vect_A, vect_B, cross_P);

    // random initial velocities magnitudes
    body[i][X_VEL] = drand48() * 100 * cross_P[0];
    body[i][Y_VEL] = drand48() * 100 * cross_P[1];
    body[i][Z_VEL] = drand48() * 100 * cross_P[2];
  }

  // print out initial positions in PDB format
  printf("MODEL %8d\n", 0);
  for (int i = 0; i < N; i++) {
    printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
           "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
  }
  printf("TER\nENDMDL\n");

  // step through each time step
  for (int t = 0; t < tmax; t++) {
    // force calculation

    // TODO: initialize forces to zero
    for (int i = 0; i < N; i++) {
      // What should I put here?
      Fx_dir[i] = 0.0; 
      Fy_dir[i] = 0.0; 
      Fz_dir[i] = 0.0; 
    }

    for (int x = 0; x < N; x++) {  // force on body x due to
      for (int i = 0; i < N; i++) {   // all other bodies
	// position differences in x-, y-, and z-directions
        float x_diff, y_diff, z_diff;

	      if (i != x) {
	  // TODO: calculate position difference between body i and x in x-,y-, and z-directions
          x_diff = body[i][X_POS] - body[x][X_POS];
          y_diff = body[i][Y_POS] - body[x][Y_POS];
          z_diff = body[i][Z_POS] - body[x][Z_POS];
          
	  // calculate distance (r)
          float rr = (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
          float r = sqrt(rr);

            // force between bodies i and x
          float F = 0;

            // if sufficiently far away, apply gravitation force
          if (r > 50.0) {
              // TODO: compute gravitational force between body i and x
            F = -1.0*(G*body[i][MASS]*body[x][MASS]) / rr; 
            norm(x_diff, y_diff, z_diff); 
            Fx_dir[x] += x_diff*F; 
            Fy_dir[x] += y_diff*F; 
            Fz_dir[x] += z_diff*F; 
          }  
        }
      }
    }

    // update postions and velocity in array
    for (int i = 0; i < N; i++) {

        // TODO: update velocities
      body[i][X_VEL] = body[i][X_VEL] + Fx_dir[i]*dt / body[i][MASS]; 
      body[i][Y_VEL] = body[i][Y_VEL] + Fy_dir[i]*dt / body[i][MASS]; 
      body[i][Z_VEL] = body[i][Z_VEL] + Fz_dir[i]*dt / body[i][MASS]; 

	// TODO: update positions
      body[i][X_POS] = body[i][X_POS] + body[i][X_VEL]*dt; 
      body[i][Y_POS] = body[i][Y_POS] + body[i][Y_VEL]*dt; 
      body[i][Z_POS] = body[i][Z_POS] + body[i][Z_VEL]*dt; 
    }

    // print out positions in PDB format
    printf("MODEL %8d\n", t+1);
    for (int i = 0; i < N; i++) {
	printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
               "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
    }
    printf("TER\nENDMDL\n");
  }  // end of time period loop
}
