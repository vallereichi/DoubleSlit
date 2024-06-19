#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <fstream>

#include <Eigen/Sparse>

typedef std::complex<double> complex;
typedef Eigen::SparseMatrix<complex> SpMat;
typedef Eigen::Triplet<complex> T;



void writeToFile(double* array, int arraySize, std::string filepath, int mode);

complex* createArray (int L, double dx);

complex* gaussWavePacket2D (complex* xarray, complex* yarray, int Nx, int Ny, double x0, double y0, double kx, double sigma);

double* getProbDensities (complex* array, int arraySize);

SpMat buildmatrixA (double* V, int Nx, int Ny, double dx, complex r,double dt);

double* setPotential (int Nx, int Ny, double V0, int i0, int i1, int j0, int j1, int j2, int j3);

complex* doubleSlitInteraction (complex* psi, int Nx, int Ny, int i0, int i1, int j0, int j1, int j2, int j3);

std::vector<Eigen::VectorXcd> CrankNicolson2D (SpMat A, SpMat M, complex* psi_init, int Nx, int Ny, int Nt, double dt, int i0, int i1, int j0, int j1, int j2, int j3);


