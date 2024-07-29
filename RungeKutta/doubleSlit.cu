#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <algorithm>

#include <cuda_runtime.h>
#define threadsPerBlock 32

typedef std::complex<double> complex;


void writeToCSV (double *array, int arraySize, std::string name, std::string filepath)
{
    std::ofstream file;
    file.open(filepath, std::ios::out | std::ios::app);
    file << "\n" << name << ",";

    for (int i = 0; i < arraySize; i++)
    {
        if (i != arraySize - 1) file << array[i] << ",";
        else file << array[i];
    }
}



double* getProbDensities (complex* array, int arraySize)
{
    double* pd = new double[arraySize];
    for (int i = 0; i < arraySize; i++) {pd[i] = std::norm(array[i]);}

    return pd;
}


void setEdgesToZero (complex* Psi, int length_x, int length_y)
{
    for (int j = 0; j < length_y; j++)
    {
        for (int i = 0; i < length_x; i++)
        {
            if (i == 0 || i == length_x -1 || j == 0 || j == length_y -1) {Psi[j * length_y + i] = 0;}
        }
    }   
}


cudaError_t cudaCheckError()
{
    cudaError_t code = cudaGetLastError();

    if (code != cudaSuccess)
    {
	    printf("%s: %s\n", cudaGetErrorName(code), cudaGetErrorString(code));
    }

    return code;
}


/*
RK4step performs one step of the Runge Kutta method of 4th order
The steps involved are:
    -> k1 = H * Psi_t
    -> k2 = H * (Psi_t + 0.5 * k1)
    -> k3 = H * (Psi_t + 0.5 * k2)
    -> k4 = H * (Psi_t + k3)

parameters:
    -Mat        -> Matrix representing the Hamiltonian
    -Vec        -> Vector holding the current wave packet
    -length_x   -> x Dimension of the matrix
    -length_y   -> y Dimension of the matrix
    -k_now      -> k vector used to calculate the next k vector
    -k_next     -> k vector that needs to be calculated
    -scale_k    -> weight of the k vector (e.g. in step two and three scale_k = 0.5)

NOTE: The Dimension of the Matrix in the actual implementation will be of size (Nx*Ny)² where Nx and Ny represent 
the number of grid points on the respective axis

*/
__global__ void RKstep(complex *Mat, complex *Vec, size_t length_x, size_t length_y, complex *k_now, complex *k_next, double scale_k)
{
    int j = blockIdx.x * gridDim.x + threadIdx.x;
    double sum = 0;
    
    if (j < length_y*length_x)
    {
        for (int i = 0; i < length_y*length_x; i++)
        {
            sum = sum + Mat[j * length_x*length_y + i] * (Vec[i] + scale_k * k_now[i]);
        }
        k_next[j] = sum;
        sum = 0;
    }
}


int main()
{
    //imaginary unit "i"
    complex I(0.0, 1.0);

    //Parameters
    float L = 10.0;                                    //well of width L
    double dx = 0.1;                                   //spatial step size
    double dt = pow(dx, 2) / 4;                        //temporal step size
    int Nx = floor(L/dx) + 1;                          //Number of points on the x axis
    int Ny = floor(L/dx) + 1;                          //Number of points on the y axis
    int N = Nx * Ny;   
    int Nt = 100;                                      //Number of time steps
    complex r = I / (pow(dx, 2.0));                    //Constant to simplify expressions

    //initial position of the center of the gaussian wave packet
    double x0 = L/5;
    double y0 = L/2;

    //momentum and width of the initial wave packet
    double kx = 15.0 * M_PI;
    double sigma = 0.75;

    //parameters of the double slit
    double w = 0.2;                                   //width of the double slit
    double s = 0.8;                                   //seperation between the slits
    double a = 0.4;                                   //aperture of the slits
    double V0 = 200.0;                                //constant value of the potential


    //indices that parameterize the double slit
    //horizontal axis
    int i0 = int(1 / (2 * dx) * (L - w));
    int i1 = int(1 / (2 * dx) * (L + w));
    //vertical axis
    int j0 = int(1 / (2 * dx) * (L + s) + a / dx);
    int j1 = int(1 / (2 * dx) * (L + s));
    int j2 = int(1 / (2 * dx) * (L - s));
    int j3 = int(1 / (2 * dx) * (L - s) - a / dx);


    std::cout << "simulation parameters:" << std::endl;
    std::cout << "L=" << L << std::endl;
    std::cout << "Nx=" << Nx << std::endl;
    std::cout << "Ny=" << Ny << std::endl;
    std::cout << "dx=" << dx << std::endl;
    std::cout << "dt=" << dt << std::endl;

    //create initial wave packet
    double *xarray = (double*)malloc(Nx * sizeof(double));
    double *yarray = (double*)malloc(Ny * sizeof(double));
    for (int i = 0; i < Nx; i++) {xarray[i] = i * dx;}
    for (int j = 0; j < Nx; j++) {yarray[j] = j * dx;}

    complex *psi0 = (complex*)malloc(Nx * Ny * sizeof(complex));
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            if (i == 0 || i == Nx -1 || j == 0 || j == Ny -1) {psi0[j * Ny + i] = 0;}
            else{
                psi0[j * Ny + i] = \
                exp(I * kx * (xarray[i] - x0)) * \
                exp(-1.0 / (2.0 * pow(sigma, 2)) * (pow(xarray[i] - x0, 2) + pow(yarray[j] - y0, 2)));
            }
        }
    }
    double* probDensities = getProbDensities(psi0, Nx*Ny);
    std::cout << "initial wavepacket set with parameters:" << std::endl;
    std::cout << "x_0=" << x0 << ", y_0=" << y0 << std::endl;
    std::cout << "sigma=" << sigma << std::endl;
    std::cout << "k_x=" << kx << std::endl;



    //set the potential
    double *V = (double*)malloc(Nx * Ny * sizeof(double));
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i< Nx; i++)
        {
            if ((i >= i0 && i <= i1) && ((j <= j3) || (j > j2 && j <j1) || (j >= j0))) {V[j * Nx + i] = V0;}
            else {V[j * Nx + i] = 0;}
        }
    }
    std::cout << "potential set..." << std::endl;



    //create the Hamiltonian
    complex *H = (complex*)malloc(N * N * sizeof(complex));
    for (int k = 0; k < N; k++)
    {
        int j = floor(k/Nx);
        int i = k % Nx;

        //main diagonal
        H[k * N + k] = -4.0 * r - I * V[j * Nx + i];

        //upper main diagonal
        if (k + 1 <= N) {H[k * N + k + 1] = r;}
        //lower main diagonal
        if (k - 1 >= 0) {H[k * N + k - 1] = r;}
        //upper lone diagonal
        if (k + Nx + 1 <= N) {H[k * N + k + Nx + 1] = r;}
        //lower lone diagonal
        if (k - Nx - 1 >= 0) {H[k * N + k - Nx - 1] = r;}
    }
    complex *d_H = NULL;
    cudaMalloc(&d_H, N * N * sizeof(complex));
    cudaMemcpy(d_H, H, N * N * sizeof(complex), cudaMemcpyHostToDevice);
    cudaCheckError();
    std::cout << "Hamiltonian set..." << std::endl;


    //write data to output file
    std::ofstream file;
    file.open("./output/output.csv");
    file.clear();
    writeToCSV(V, Nx * Ny, "Potential", "./output/output.csv");
    writeToCSV(probDensities, Nx * Ny, "Psi-0", "./output/output.csv");
    delete[] probDensities;
    free(xarray);
    free(yarray);
    std::cout << "data has been written to output file..." << std::endl;


    //run Runge Kutta for all time steps
    std::cout << "Running Runge Kutta simulation..." << std::endl;
    complex *psiT = (complex*)malloc(Nx * Ny * sizeof(complex));
    std::copy(psi0, psi0 + Nx * Ny, psiT);
    free(psi0);
    complex *d_psiT = NULL;
    cudaMalloc(&d_psiT, N * sizeof(complex));
    cudaMemcpy(d_psiT, psiT, N * sizeof(complex), cudaMemcpyHostToDevice);
    cudaCheckError();

    for (int t = 0; t < Nt; t++)
    {
        std::ostringstream timeStep;
        timeStep << "Psi-" << t + 1;

        complex *k0 = (complex*)malloc(Nx * Ny * sizeof(complex));
        complex *k1 = (complex*)malloc(Nx * Ny * sizeof(complex));
        complex *k2 = (complex*)malloc(Nx * Ny * sizeof(complex));
        complex *k3 = (complex*)malloc(Nx * Ny * sizeof(complex));
        complex *k4 = (complex*)malloc(Nx * Ny * sizeof(complex));


        complex *d_k0 = NULL, *d_k1 = NULL, *d_k2 = NULL, *d_k3 = NULL, *d_k4 = NULL;

        cudaMalloc(&d_k0, N * sizeof(complex));
        cudaMalloc(&d_k1, N * sizeof(complex));
        cudaMalloc(&d_k2, N * sizeof(complex));
        cudaMalloc(&d_k3, N * sizeof(complex));
        cudaMalloc(&d_k4, N * sizeof(complex));

        cudaMemcpy(d_k0, k0, N * sizeof(complex), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k1, k1, N * sizeof(complex), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k2, k1, N * sizeof(complex), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k3, k1, N * sizeof(complex), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k4, k1, N * sizeof(complex), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        cudaCheckError();

        int threads = threadsPerBlock;
        int blocks = ceil((Ny + threads - 1) / threads);

        RKstep<<<blocks, threads>>>(d_H, d_psiT, Nx, Ny, d_k0, d_k1, 0.0);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_H, d_psiT, Nx, Ny, d_k1, d_k2, 0.5 * dt);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_H, d_psiT, Nx, Ny, d_k2, d_k3, 0.5 * dt);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_H, d_psiT, Nx, Ny, d_k3, d_k4, 1.0 * dt);
        cudaDeviceSynchronize();
        cudaCheckError();

        cudaMemcpy(k1, d_k1, N * sizeof(complex), cudaMemcpyDeviceToHost);
        cudaMemcpy(k2, d_k2, N * sizeof(complex), cudaMemcpyDeviceToHost);
        cudaMemcpy(k3, d_k3, N * sizeof(complex), cudaMemcpyDeviceToHost);
        cudaMemcpy(k4, d_k4, N * sizeof(complex), cudaMemcpyDeviceToHost);
        cudaCheckError();


        for (int i = 0; i < Ny * Nx; i++) 
        {
            psiT[i] = psiT[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        setEdgesToZero(psiT, Nx, Ny);

        double *pd = getProbDensities(psiT, Nx * Ny);
        writeToCSV(pd, Nx * Ny, timeStep.str(), "./output/output.csv");
        delete[] pd;

        cudaFree(d_k0);
        cudaFree(d_k1);
        cudaFree(d_k2);
        cudaFree(d_k3);
        cudaFree(d_k4);

        free(k0);
        free(k1);
        free(k2);
        free(k3);
        free(k4);
    }

    file.close();
    cudaFree(d_H);
    cudaFree(d_psiT);
    free(H);
    free(V);
    free(psiT);
}