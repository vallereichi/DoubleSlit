#include "wavePacket.h"

void writeToFile(double* array, int arraySize, std::string filepath, int mode)
{
    std::ofstream file;
   
    if (mode == 0)     
        {
            file.open(filepath);
            for (int i = 0; i < arraySize; i++) 
            {
                file << array[i] << ", ";
            }
        }
    else if (mode == 1)
        {
            file.open(filepath, std::ios::out | std::ios::app);
            for (int i = 0; i < arraySize; i++) 
            {
                file << array[i] << ", ";
            }
        }    
}



complex* createArray (int L, double dx)
{
    int N = floor(L/dx) + 1;
    complex* array = new complex[N];

    for (int i = 0; i < N; i++) {array[i] = i * dx;}

    return array;
}



complex* gaussWavePacket2D (complex* xarray, complex* yarray, int Nx, int Ny, double x0, double y0, double kx, double sigma)
{
    complex* wavePacket = new complex[Nx*Ny];
    complex I(0.0, 1.0);
    int idx = 0;

    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {wavePacket[idx] = 0;}
            else 
            {
                wavePacket[idx] = \
                    exp(I * kx * (xarray[i] - x0)) * \
                    exp(-1.0 / (2.0 * pow(sigma, 2)) * (pow(xarray[i] - x0, 2) + pow(yarray[j] - y0, 2)));
            }
        
            idx++;
        }
    }

    return wavePacket;
}



double* getProbDensities (complex* array, int arraySize)
{
    double* pd = new double[arraySize];
    for (int i = 0; i < arraySize; i++) {pd[i] = std::norm(array[i]);}

    return pd;
}



SpMat buildmatrixA (double* V, int Nx, int Ny, double dx, complex r,double dt)
{
    std::vector<T> TripletList;
    TripletList.reserve(6*Nx*Ny);
    int N = Nx*Ny;
    complex I(0.0, 1.0);
    SpMat A(N, N);

    for (int k = 0; k < N; k++)
    {
        //k = (i-1)*(Nx) + (j-1)
        int i = 1 + floor(k / Nx);
        int j = 1 + k % Nx; 

        //main diagonal
        TripletList.push_back(T(k, k, 1.0 + 4.0 * r + I * dt / 2.0 * V[j * Nx + i]));

        //lower lone diagonal
        if (i != 1) {TripletList.push_back(T(k, (i - 2)*Nx + j -1, -r));}

        //upper lone diagonal
        if (i != Nx) {TripletList.push_back(T(k, i * Nx + j - 1, -r));}

        //lower main diagonal
        if (j != 1) {TripletList.push_back(T(k, k-1, -r));}

        //upper main diagonal
        if (j != Nx) {TripletList.push_back(T(k, k + 1, -r));}
    }

    A.setFromTriplets(TripletList.begin(), TripletList.end());

    return A;
}



double* setPotential (int Nx, int Ny, double V0, int i0, int i1, int j0, int j1, int j2, int j3)
{
    double* V = new double[Nx*Ny];

    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            if ((i >= i0 && i <= i1) && ((j <= j3) || (j > j2 && j < j1) || (j >= j0))) {V[j * Nx + i] = V0;}
            else {V[j * Nx + i] = 0 ;}
        }
    }

    return V;
}


complex* doubleSlitInteraction (complex* psi, int Nx, int Ny, int i0, int i1, int j0, int j1, int j2, int j3)
{
    complex* psiNew = new complex[Nx*Ny];

    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            if ((i >= i0 && i <= i1) && ((j <= j3) || (j > j2 && j < j1) || (j >= j0))) {psiNew[j * Nx + i] = 0;}
            else {psiNew[j * Nx +i] = psi[j * Nx +i];}
        }
    }

    return psiNew;
}


std::vector<Eigen::VectorXcd> CrankNicolson2D (SpMat A, SpMat M, complex* psi_init, int Nx, int Ny, int Nt, double dt, int i0, int i1, int j0, int j1, int j2, int j3)
{
    A.makeCompressed();
    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    std::vector<Eigen::VectorXcd> psi;
    Eigen::VectorXcd psi0 = Eigen::Map<Eigen::VectorXcd>(psi_init, Nx*Ny);
    psi.push_back(psi0);

    for (int t = 0; t < Nt; t++)
    {
        Eigen::VectorXcd b = M * psi[t];
        Eigen::VectorXcd x = solver.solve(b);

        complex* tmp = doubleSlitInteraction(x.data(), Nx, Ny, i0, i1, j0, j1, j2, j3);
        psi.push_back(Eigen::Map<Eigen::VectorXcd>(tmp, Nx*Ny));
        
    }

    return psi;
}

