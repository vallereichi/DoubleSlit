#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <fstream>

typedef std::complex<double> complex;

void writeToFile(double* array, int arraySize, std::string filepath, int mode);

complex* createArray (int L, double dx);

complex* gaussWavePacket2D (complex* xarray, complex* yarray, int Nx, int Ny, double x0, double y0, double kx, double sigma);

double* getProbDensities (complex* array, int arraySize);

double* setPotential (int Nx, int Ny, double V0, int i0, int i1, int j0, int j1, int j2, int j3);

complex* Hamiltonian (double* V, int Nx, int Ny, double dx, complex r, double dt);

complex* MatVecMul (complex* Mat, complex* Vec, int vecSize);

complex* VecAdd (complex* vecA, complex* vecB, int vecSize);

complex* VecScale (complex* Vec, double scale, int vecSize);

std::vector<complex*> RungeKutta (complex* H, complex* psi0, int Nx, int Ny, int Nt, double dt);