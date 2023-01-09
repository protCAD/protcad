#ifdef __CUDA__
#include <cuda_runtime.h>
#endif
#include <cstdio>
#include <iostream>
#include <cmath>
#ifndef ENERGY_H
#define ENERGY_H

const int threads_per_block(256);

//////////Memory allocation////////////////////////////////////////////////////////////
void loadEnergyDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, double* eps_h, double* chg_h, double* vol_h, int* bon_h, double *E_h, int N);
void loadClashDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, int* bon_h, int* clash_h, int N);
void loadAllDeviceMem(double* x_h, double* y_h, double* z_h, double* rad_h, double* eps_h, double* chg_h, double* vol_h, int* clash_h, int* bon_h, double *E_h, int N);
void freeEnergyDeviceMem();
void freeClashDeviceMem();
void freeAllDeviceMem();

//////////Functions///////////////////////////////////////////////////////////////////
void calcEnergies(double* x_h, double* y_h, double* z_h, double* E_h, int N);
void calcClashes(double* x_h, double* y_h, double* z_h, int* clash_h, int N);

//////////Utils//////////////////////////////////////////////////////////////////////
#define check(ans) { _check((ans), __FILE__, __LINE__); }

#endif
