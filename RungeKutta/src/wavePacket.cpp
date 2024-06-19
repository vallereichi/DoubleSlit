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



complex* Hamiltonian (double* V, int Nx, int Ny, double dx, complex r, double dt)
{
    int N = Nx * Ny;
    complex* H = new complex[N*N];
    complex I(0.0, 1.0);

    for (int k = 0; k < N; k++)
    {
        int j = 1  + floor(k/Nx);
        int i = 1  + k % Nx;

        //main diagonal
        H[k * N + k] = 1.0 + 4.0 * r + I * dt / 2.0 * V[j * Nx + i];

        //upper main diagonal
        if (k + 1 <= N) {H[k * N + k + 1] = -r;}
        //lower main diagonal
        if (k - 1 >= 0) {H[k * N + k - 1] = -r;}
        //upper lone diagonal
        if (j * Nx + i - 1 <= N) {H[k * N + (j * Nx + i - 1)] = -r;}
        //lower lone diagonal
        if ((j - 2) * Nx + i - 1 >= 0) {H[k * N + ((j - 2) * Nx + i - 1)] = -r;}

    }

    return H;
}


complex* MatVecMul (complex* Mat, complex* Vec, int vecSize)
{
    complex* result = new complex[vecSize];

    for (int j = 0; j < vecSize; j++)
    {
        for (int i = 0; i < vecSize; i++)
        {
            result[j] = result[j] + Mat[j * vecSize + i] * Vec[i];
        }
    }

    return result;
}


complex* VecAdd (complex* vecA, complex* vecB, int vecSize)
{
    complex* result = new complex[vecSize];

    for (int i = 0; i < vecSize; i++)
    {
        result[i] = vecA[i] + vecB[i];
    }

    return result;
}


complex* VecScale (complex* Vec, double scale, int vecSize)
{
    complex* result = new complex[vecSize];

    for (int i = 0; i < vecSize; i++) {result[i] = scale * Vec[i];}

    return result;
}


std::vector<complex*> RungeKutta (complex* H, complex* psi0, int Nx, int Ny, int Nt, double dt)
{
    std::vector<complex*> psi_t;
    int N = Nx * Ny;

    psi_t.push_back(psi0);

    for (int t = 0; t < Nt; t++)
    {
        //complex* k1 = MatVecMul(H, psi_t[t], N);

        complex* k1 = new complex[N];
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < N; i++)
            {
                k1[j] = k1[j] + H[j * N + i] * psi_t[t][i];
            }
        }

        //complex* k2 = MatVecMul(H, VecAdd(psi_t[t], VecScale(k1, 0.5, N), N), N);

        complex* k2 = new complex[N];
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < N; i++)
            {
                k2[j] = k2[j] + H[j * N + i] * (psi_t[t][i] + 0.5 * k1[i]);
            }
        }

        //complex* k3 = MatVecMul(H, VecAdd(psi_t[t], VecScale(k2, 0.5, N), N), N);

        complex* k3 = new complex[N];
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < N; i++)
            {
                k3[j] = k3[j] + H[j * N + i] * (psi_t[t][i] + 0.5 * k2[i]);
            }
        }

        //complex* k4 = MatVecMul(H, VecAdd(psi_t[t], k3, N), N);

        complex* k4 = new complex[N];
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < N; i++)
            {
                k4[j] = k4[j] + H[j * N + i] * (psi_t[t][i] + 0.5 * k3[i]);
            }
        }

        complex* tmp = new complex[N];

        for (int i = 0; i < N; i++)
        {
            tmp[i] = psi_t[t][i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        }

        psi_t.push_back(tmp);

        //delete[] k1;
        //delete[] k2;
        //delete[] k3;
        //delete[] k4;
        //delete[] tmp;
    }

    return psi_t;
}
