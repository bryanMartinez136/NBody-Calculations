#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>

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
#define dt 0.05  // time interval
#define W 7

float body[10000][7]; // data array of bodies


__global__ void init(unsigned int seed, curandState_t* states, int n) {

    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < n) {
  curand_init(seed,
              i,
              0,
              &states[i]);
    }
}
__device__ void crossProduct(float vect_A[], float vect_B[], float cross_P[]) { 
  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]; 
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]; 
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]; 
}
__device__ void norm(float &x, float &y, float &z) {
  float mag = sqrt(x*x+y*y+z*z);
  x/=mag; y/=mag; z/=mag;
}
__global__ void randoms (curandState_t* states, float* boddy, int n) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    float vect_A[3], vect_B[3], cross_P[3];

    if (i < n && i>0) {
      // get x, y, and z positions
      boddy[i* W+MASS] = 0.001; // 0
      boddy[i* W+X_POS] = (curand_uniform(&states[i]))*300 - (150+0.999); // 1
      boddy[i* W+Y_POS] = (curand_uniform(&states[i]))* 300 -(150+0.999); // 2
      boddy[i* W+Z_POS] = (curand_uniform(&states[i]))*300 -(150+0.999); // 3
      // compute norm
      vect_A[0]= boddy[i* W+X_POS];
      vect_A[1]= boddy[i* W+Y_POS];
      vect_A[2]= boddy[i* W+Z_POS];
      norm(vect_A[0], vect_A[1], vect_A[2]);
      // get the cross product
      vect_B[0]= 0.0; vect_B[1]= 0.0; vect_B[2]= 1.0;
      cross_P[0] = 0.0; cross_P[1] = 0.0; cross_P[2] = 0.0; 
      crossProduct(vect_A, vect_B, cross_P);

       // random initial velocities magnitudes
      boddy[i*W+X_VEL] = (curand_uniform(&states[i]))*(100+0.999) *cross_P[0];
      boddy[i*W+Y_VEL] = (curand_uniform(&states[i]))*(100+0.999) *cross_P[1];
      boddy[i*W+Z_VEL] = (curand_uniform(&states[i]))*(100+0.999) *cross_P[2];


    }

}
__global__ void forces(float* dev_body, float* dev_fx, float* dev_fy, float* dev_fz, int n){

  float x_diff, y_diff, z_diff; 
  int i = threadIdx.x + blockIdx.x * blockDim.x; 

  for(int x = 0 ; x < n; x++){
    if(x != i && i < n){
      x_diff = dev_body[i*W+X_POS] - dev_body[x*W+X_POS];
      y_diff = dev_body[i*W+Y_POS] - dev_body[x*W+Y_POS];
      z_diff = dev_body[i*W+Z_POS] - dev_body[x*W+Z_POS];
  
	    // calculate distance (r)
      float rr = (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
      float r = sqrt(rr);

      // force between bodies i and x
      float F = 0;

      // if sufficiently far away, apply gravitation force
      if (r > 50.0) {
        // TODO: compute gravitational force between body i and x
        F = -1.0*(G*dev_body[i*W+MASS]*dev_body[x*W+MASS]) / rr; 
        norm(x_diff, y_diff, z_diff); 
        dev_fx[i] += (x_diff/r)*F; 
        dev_fy[i] += (y_diff/r)*F; 
        dev_fz[i] += (z_diff/r)*F; 
      } 
    }
  }
}

__global__ void update(float* dev_body, float* dev_fx, float* dev_fy, float* dev_fz, int n){

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < n) {

    dev_body[i*W+X_VEL] = dev_body[i*W+X_VEL] + dev_fx[i]*dt / dev_body[i*W+MASS]; 
    dev_body[i*W+Y_VEL] = dev_body[i*W+Y_VEL] + dev_fy[i]*dt / dev_body[i*W+MASS]; 
    dev_body[i*W+Z_VEL] = dev_body[i*W+Z_VEL] + dev_fz[i]*dt / dev_body[i*W+MASS]; 

    // TODO: update positions
    dev_body[i*W+X_POS] = dev_body[i*W+X_POS] + dev_body[i*W+X_VEL]*dt; 
    dev_body[i*W+Y_POS] = dev_body[i*W+Y_POS] + dev_body[i*W+Y_VEL]*dt; 
    dev_body[i*W+Z_POS] = dev_body[i*W+Z_POS] + dev_body[i*W+Z_VEL]*dt;
  } 

}

__global__ void initForce(float* dev_fx,float*  dev_fy,float*  dev_fz, int n){
  int i = threadIdx.x + blockDim.x*blockIdx.x; 
  if(i < n){
    dev_fx[i] = 0; 
    dev_fy[i] = 0; 
    dev_fz[i] = 0; 
  }
}


int main(int argc, char **argv) {

  int tmax = 0;
  float Fx_dir[N], Fy_dir[N], Fz_dir[N];
  float * dev_fx, *dev_fy, *dev_fz;

  if (argc != 2) {
    fprintf(stderr, "Format: %s { number of timesteps }\n", argv[0]);
    exit (-1);
  }

  tmax = atoi(argv[1]);
  dim3 dimBlock(1024);
  dim3 dimGrid((int)ceil((float)N / 1024)); 

  if(tmax<0){
    fprintf(stderr, "No negative values for time allowed\n");

  }

  // assign each body a random initial positions and velocities
  // black hole at the center
  body[0][MASS] = 4000.0; body[0][X_POS] = 0.0;body[0][Y_POS] = 0.0;
  body[0][Z_POS] = 0.0;body[0][X_VEL] = 0.0;body[0][Y_VEL] = 0.0;body[0][Z_VEL] = 0.0;
  
// PARALLELIZED THE RANDOM INITIALIZATION OF MASS AND POSITIONS

  curandState_t* states;
  // allocate space on GPU for random states
  cudaMalloc((void**) &states, N*7*sizeof(curandState_t));
  
  /* invoke the GPU to initialize all of the random states */
  init<<<dimGrid, dimBlock>>>(time(0), states, N*7);
  cudaDeviceSynchronize();

  // allocate array of unsigned ints on CPU and GPU
  float* body_arr;
  cudaMalloc((void**) &body_arr, N*7*sizeof(float));

  // obtain a uniformly random distriubtion of integers, maximum N
  randoms<<<dimGrid, dimBlock>>>(states, body_arr, N*7);
  cudaDeviceSynchronize();
  cudaMemcpy(body, body_arr, N*7*sizeof(float), cudaMemcpyDeviceToHost);
// finished coppying all the initial values to body

  // print out initial positions in PDB format
  printf("MODEL %8d\n", 0);
  for (int i = 0; i < N; i++) {
    printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
           "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
  }
  printf("TER\nENDMDL\n");


// FORCE CALCULATION !!!!

  cudaMalloc((void**) &dev_fx, N*sizeof(float));
  cudaMalloc((void**) &dev_fy, N*sizeof(float));
  cudaMalloc((void**) &dev_fz, N*sizeof(float));

  cudaEvent_t start;
  cudaEventCreate(&start);
  cudaEvent_t stop;
  cudaEventCreate(&stop);

  // start timer
  cudaEventRecord(start,0);

  for (int t = 0; t < tmax; t++) {
    // TODO: initialize forces to zero

    for(int i = 0; i < N; i++) {
      Fx_dir[i] = 0.0; 
      Fy_dir[i] = 0.0; 
      Fz_dir[i] = 0.0; 
      
    }

    // initForce<<<dimGrid, dimBlock>>>(dev_fx,dev_fy,dev_fz, N); 
    cudaMemcpy(dev_fx, Fx_dir, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_fy, Fy_dir, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_fz, Fz_dir, sizeof(float), cudaMemcpyHostToDevice);

    // cudaDeviceSynchronize(); 

    forces<<<dimGrid, dimBlock>>>(body_arr, dev_fx,dev_fy,dev_fz, N); 
    cudaDeviceSynchronize();

    update<<<dimGrid, dimBlock>>>(body_arr, dev_fx, dev_fy, dev_fz, N);
    cudaDeviceSynchronize();
    cudaMemcpy(body, body_arr, N*7*sizeof(float), cudaMemcpyDeviceToHost);

    // print out positions in PDB format
    printf("MODEL %8d\n", t+1);
    for (int i = 0; i < N; i++) {
	printf("%s%7d  %s %s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.3f\n",
               "ATOM", i+1, "CA ", "GLY", "A", i+1, body[i][X_POS], body[i][Y_POS], body[i][Z_POS], 1.00, 0.00);
    }
    printf("TER\nENDMDL\n");

  }
  
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);

  float diff;
  cudaEventElapsedTime(&diff, start, stop);
  printf("time: %f ms\n", diff);

  // deallocate timers
  cudaEventDestroy(start);
  cudaEventDestroy(stop);


  cudaFree(states);
  cudaFree(body_arr);
  cudaFree(dev_fx); 
  cudaFree(dev_fy); 
  cudaFree(dev_fz); 
}
