#include "energy.h"

// pointer to data copied to device
double* x_d; 
double* y_d; 
double* z_d; 
double* rad_d; 
double* eps_d; 
double* chg_d;
double* vol_d; 
double* dis_d;
int* bon_d;
int* clash_d;
double* E_d;

// Constants on device
__constant__ double watRad = 1.4; 
__constant__ double watDia = 4.35; //effective diameter between water molecules
__constant__ double watVol = 107.31; // vol of effective diameter
__constant__ double watPol = 1.47; // water polarization (Murphy WF. J. Chem. Phys. 1977;67:5877–5882)
__constant__ double watEps = 0.15200; 
__constant__ double pi = 3.1415926535;
__constant__ double kc = 332.0636;
__constant__ double kb = 0.0019872041;
__constant__ double t = 300.0;
__constant__ double v = 4.188; // 4/3*pi

// Error checking
inline void _check(cudaError_t code, const char *file, int line)
{
  if (code != cudaSuccess) {
    fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
    exit(code);}
}

///////////Kernels//////////////////////////////////////////////////////////////////

// Kernel function for calculating distances and atom environment volumes
__global__ void calcDistance(double* x, double* y, double* z, double* rad, double* vol, double* dist, int N)
{
  // Compute the global thread index
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Compute the indices of the two points for this thread
  int i = idx / N;
  int j = idx % N;
  int index = i*(N-1)-(i-1)*i/2+j-i-1;

  if (j > i){
    // Calculate distance between i and j
    double distance = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));

    // Estimate volume around atoms occupied by other atoms
    // exposure to solvent can be estimated without a surface calc by estimating cavity volume around atom
    double volJ = v*pow(rad[j],3); double volI = v*pow(rad[i],3);
    if (distance < rad[i]+watDia) {atomicAdd(&vol[i], volJ);}
    if (distance < rad[j]+watDia) {atomicAdd(&vol[j], volI);}

    // save dist to array 
    dist[index] = distance;
  }
}

// Kernel function for calculating distances and atom environment volumes
__global__ void calcClash(double* x, double* y, double* z, double* rad, int* bon, int* clash, int N)
{
  // Compute the global thread index
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Compute the indices of the two points for this thread
  int i = idx / N;
  int j = idx % N;
  int index = i*(N-1)-(i-1)*i/2+j-i-1;

  if (j > i){
    if (bon[index] == 0){
      // Calculate distanceSquared between i and j
      double distanceSq = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]);

      //count clashes between i and j
      if (distanceSq < (rad[i]+rad[j])*(rad[i]+rad[j])) {
        //printf("%d\n",clash[i]);
        atomicAdd(&clash[i],1); atomicAdd(&clash[j],1);
      }
    }
  }
}

// Kernel function for calculating the energy between two atoms
__global__ void calcEnergy(double* rad, double* eps, double* chg, double* vol, double* dist, int* bon, double *E, int N)
{
  // Compute the global thread index
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  // Compute the indices of the two points for this thread
  int i = idx / N;
  int j = idx % N;
  int index = i*(N-1)-(i-1)*i/2+j-i-1;
  

  if (j > i){
    if (bon[index] == 0){
      //Calculate VDW Energies
      double vdw = (sqrt(eps[i]*eps[j])) * (pow(((rad[i]+rad[j])/dist[index]),12) - (2*pow(((rad[i]+rad[j])/dist[index]),6)));
      atomicAdd(E, vdw);

      // Estimate the number of waters occupying remaining shell volume of i to be used for local dielectric
      double polI=0.0; double watersI = 0.0;
      double shellVolI = v*pow((rad[i]+watDia),3);
      double envVolI = (vol[i]+(v*pow(rad[i],3)))/2;
      double waterVolI = shellVolI-envVolI;
      if (waterVolI > 0){watersI = int(waterVolI/watVol); polI = watersI*watPol;}

      // Estimate the number of waters occupying remaining shell volume of j to be used for local dielectric
      double polJ=0.0; double watersJ = 0.0;
      double shellVolJ = v*pow((rad[j]+watDia),3);
      double envVolJ = (vol[j]+(v*pow(rad[j],3)))/2;
      double waterVolJ = shellVolJ-envVolJ;
      if (waterVolJ > 0){watersJ = int(waterVolJ/watVol); polJ = watersJ*watPol;}

      //Calculate the effective dielectric with the Lorentz local field correction
      double dielectricI =2+(8*pi/3)*(polI)/1-(4*pi/3)*(polI);
      double dielectricJ =2+(8*pi/3)*(polJ)/1-(4*pi/3)*(polJ);
      double dielectric = (dielectricI+dielectricJ)/2;
      
      //Calculate Electrostatic Energies
      double ele = (kc * (chg[i] * chg[j]) / dist[index]) / dielectric;
      atomicAdd(E, ele);
    }
  }
  // Calculate solvation energy of each atom (i)
  if (i == j){

    // Estimate the number of waters occupying remaining shell volume of i to be used for solvation
    double shellVolI = v*pow((rad[i]+watDia),3);
    double envVolI = (vol[i]+(v*pow(rad[i],3)))/2;
    double waterVolI = shellVolI-envVolI;
    if (waterVolI > 0){
      double watersI = int(waterVolI/watVol); double polI = watersI*watPol;
    
      //Calculate the effective dielectric with the Lorentz local field correction
      double dielectricI =2+(8*pi/3)*(polI)/1-(4*pi/3)*(polI);
      
      // Calculate electrostatic solvation using born approximation
      double eleSolv = -(kc/2)*(chg[i]*chg[i])/((rad[i]+watRad)*dielectricI)*watersI;
      atomicAdd(E, eleSolv);

      // Calculate ideal vdw interactions with waters
      double vdwSolv = (sqrt(eps[i]*watEps) * -1)*watersI;
      atomicAdd(E, vdwSolv);

      // Calculate hydrophobic effect (water entropy cost)
      double entSolv = -(kb*t)*log(pow(0.5,watersI));
      atomicAdd(E, entSolv);
    }
  }
}
///////////Kernels End//////////////////////////////////////////////////////////////////


//////////Memory allocation////////////////////////////////////////////////////////////

void loadEnergyDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, double* eps_h, double* chg_h, double* vol_h, int* bon_h, double *E_h, int N)
{
  // Allocate memory on the GPU for the arrays
  int bondingSize = N * (N - 1) / 2;
  check(cudaMalloc(&x_d,   N * sizeof(double))); check(cudaMalloc(&y_d,   N * sizeof(double))); check(cudaMalloc(&z_d,   N * sizeof(double)));
  check(cudaMalloc(&rad_d, N * sizeof(double))); check(cudaMalloc(&eps_d, N * sizeof(double))); check(cudaMalloc(&chg_d, N * sizeof(double)));
  check(cudaMalloc(&vol_d, N * sizeof(double))); check(cudaMalloc(&bon_d, bondingSize * sizeof(int))); 
  check(cudaMalloc(&dis_d, bondingSize * sizeof(double))); check(cudaMalloc(&E_d, sizeof(double)));
  
  // Copy the coordinates, radius, epsilon, charge and bonding from the host (CPU) to the device (GPU)
  check(cudaMemcpy(x_d, x_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(y_d, y_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(z_d, z_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(rad_d, rad_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(eps_d, eps_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(chg_d, chg_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(vol_d, vol_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(bon_d, bon_h, bondingSize * sizeof(int), cudaMemcpyHostToDevice));
  check(cudaMemcpy(E_d, E_h, sizeof(double), cudaMemcpyHostToDevice));
}

void loadClashDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, int* bon_h, int* clash_h, int N)
{
  // Allocate memory on the GPU for the arrays
  int bondingSize = N * (N - 1) / 2;
  check(cudaMalloc(&x_d,   N * sizeof(double))); check(cudaMalloc(&y_d,   N * sizeof(double))); check(cudaMalloc(&z_d,   N * sizeof(double)));
  check(cudaMalloc(&rad_d, N * sizeof(double))); check(cudaMalloc(&clash_d, N * sizeof(int))); check(cudaMalloc(&bon_d, bondingSize * sizeof(int))); 
  
  // Copy the coordinates, radius, epsilon, charge and bonding from the host (CPU) to the device (GPU)
  check(cudaMemcpy(x_d, x_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(y_d, y_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(z_d, z_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(rad_d, rad_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(clash_d, clash_h, N * sizeof(int), cudaMemcpyHostToDevice));
  check(cudaMemcpy(bon_d, bon_h, bondingSize * sizeof(int), cudaMemcpyHostToDevice));
}

void loadAllDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, double* eps_h, double* chg_h, double* vol_h, int* clash_h, int* bon_h, double *E_h, int N)
{
  // Allocate memory on the GPU for the arrays
  int bondingSize = N * (N - 1) / 2;
  check(cudaMalloc(&x_d,   N * sizeof(double))); check(cudaMalloc(&y_d,   N * sizeof(double))); check(cudaMalloc(&z_d,   N * sizeof(double)));
  check(cudaMalloc(&rad_d, N * sizeof(double))); check(cudaMalloc(&eps_d, N * sizeof(double))); check(cudaMalloc(&chg_d, N * sizeof(double)));
  check(cudaMalloc(&vol_d, N * sizeof(double))); check(cudaMalloc(&clash_d, N * sizeof(int))); check(cudaMalloc(&bon_d, bondingSize * sizeof(int))); 
  check(cudaMalloc(&dis_d, bondingSize * sizeof(double))); check(cudaMalloc(&E_d, sizeof(double)));
  
  // Copy the coordinates, radius, epsilon, charge and bonding from the host (CPU) to the device (GPU)
  check(cudaMemcpy(x_d, x_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(y_d, y_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(z_d, z_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(rad_d, rad_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(eps_d, eps_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(chg_d, chg_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(vol_d, vol_h, N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(clash_d, clash_h, N * sizeof(int), cudaMemcpyHostToDevice));
  check(cudaMemcpy(bon_d, bon_h, bondingSize * sizeof(int), cudaMemcpyHostToDevice));
  check(cudaMemcpy(E_d, E_h, sizeof(double), cudaMemcpyHostToDevice));
}

void freeEnergyDeviceMem()
{
  // Free the GPU memory
  check(cudaFree(x_d)); check(cudaFree(y_d)); check(cudaFree(z_d));
  check(cudaFree(rad_d)); check(cudaFree(eps_d)); check(cudaFree(chg_d));
  check(cudaFree(bon_d)); check(cudaFree(vol_d)); check(cudaFree(dis_d)); check(cudaFree(E_d));
}

void freeClashDeviceMem()
{
  // Free the GPU memory
  check(cudaFree(x_d)); check(cudaFree(y_d)); check(cudaFree(z_d));
  check(cudaFree(rad_d)); check(cudaFree(clash_d)); check(cudaFree(bon_d)); 
}

void freeAllDeviceMem()
{
  // Free the GPU memory
  check(cudaFree(x_d)); check(cudaFree(y_d)); check(cudaFree(z_d));
  check(cudaFree(rad_d)); check(cudaFree(eps_d)); check(cudaFree(chg_d)); check(cudaFree(clash_d));
  check(cudaFree(bon_d)); check(cudaFree(vol_d)); check(cudaFree(dis_d)); check(cudaFree(E_d));
}
///////////Memory allocation End///////////////////////////////////////////////////////


//////////Functions///////////////////////////////////////////////////////////////////

void calcEnergies(double* x_h, double* y_h, double* z_h, double *E_h, int N)
{
  // Update the coordinates and starting energy from the host to the GPU
  check(cudaMemcpy(x_d, x_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(y_d, y_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(z_d, z_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(E_d, E_h, sizeof(double), cudaMemcpyHostToDevice));
  
  // Bound and invoke the distance kernel
  int tC = N * (N-1) / 2;
  int blocks = (tC+threads_per_block-1)/threads_per_block;
  calcDistance <<< blocks, threads_per_block >>> (x_d,y_d,z_d,rad_d,vol_d,dis_d,N);
  calcEnergy <<< blocks, threads_per_block >>> (rad_d,eps_d,chg_d,vol_d,dis_d,bon_d,E_d,N);
  //check(cudaPeekAtLastError());
  check(cudaDeviceSynchronize());
  
  // Copy the final energy back to the host
  check(cudaMemcpy(E_h, E_d, sizeof(double), cudaMemcpyDeviceToHost));
}

void calcClashes(double* x_h, double* y_h, double* z_h, int* clash_h, int N)
{
  // Update the coordinates and starting energy from the host to the GPU
  check(cudaMemcpy(x_d, x_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(y_d, y_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(z_d, z_h,     N * sizeof(double), cudaMemcpyHostToDevice));
  check(cudaMemcpy(clash_d, clash_h, N * sizeof(int), cudaMemcpyHostToDevice));
  
  // Bound and invoke the distance kernel
  int tC = N * (N-1) / 2;
  int blocks = (tC+threads_per_block-1)/threads_per_block;
  calcClash <<< blocks, threads_per_block >>> (x_d,y_d,z_d,rad_d,bon_d,clash_d,N);
  //check(cudaPeekAtLastError());
  check(cudaDeviceSynchronize());
  
  // Copy the final energy back to the host
  check(cudaMemcpy(clash_h, clash_d, N * sizeof(int), cudaMemcpyDeviceToHost));
}

//////////Functions end///////////////////////////////////////////////////////////////

