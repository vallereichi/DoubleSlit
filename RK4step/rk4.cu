#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <chrono>

#include <cuda_runtime.h>

#define threadsPerBlock 32

cudaError_t cudaCheckError()
{
    cudaError_t code = cudaGetLastError();

    if (code != cudaSuccess)
    {
	    printf("%s: %s\n", cudaGetErrorName(code), cudaGetErrorString(code));
    }

    return code;
}



__global__ void RKstep(double *Mat, double *Vec, size_t length_x, size_t length_y, double *k_now, double *k_next, double scale_k)
{
    int j = blockIdx.x * gridDim.x + threadIdx.x;
    double sum = 0;
    
    if (j < length_y)
    {
        for (int i = 0; i < length_y; i++)
        {
            sum = sum + Mat[j * length_x + i] * (Vec[i] + scale_k * k_now[i]);
        }
        k_next[j] = sum;
        sum = 0;
    }
}



int main ()
{
    std::cout << "Time for calculation on GPU:" << std::endl;

    std::ofstream file;
    file.open("./timecuda.txt", std::ios::out | std::ios::app);

    int Nt = 100;
    for (int t = 0; t < Nt; t++)
    {
        size_t Nx = 100*t;
        size_t Ny = 100*t;

        double *Mat = (double*)malloc(Nx * Ny * sizeof(double));
        double *Vec = (double*)malloc(Ny *sizeof(double));

        for (int i = 0; i < Nx*Ny; i++)     {Mat[i] = 1;}
        for (int i = 0; i < Ny; i++)        {Vec[i] = i + 1;}


        double *k0 = (double*)malloc(Ny * sizeof(double));
        double *k1 = (double*)malloc(Ny * sizeof(double));
        double *k2 = (double*)malloc(Ny * sizeof(double));
        double *k3 = (double*)malloc(Ny * sizeof(double));
        double *k4 = (double*)malloc(Ny * sizeof(double));




        double *d_Mat = NULL, *d_Vec = NULL, *d_k0 = NULL, *d_k1 = NULL, *d_k2 = NULL, *d_k3 = NULL, *d_k4 = NULL;

        cudaMalloc(&d_Mat, Nx*Ny * sizeof(double));
        cudaMalloc(&d_Vec, Ny * sizeof(double));
        cudaMalloc(&d_k0, Ny * sizeof(double));
        cudaMalloc(&d_k1, Ny * sizeof(double));
        cudaMalloc(&d_k2, Ny * sizeof(double));
        cudaMalloc(&d_k3, Ny * sizeof(double));
        cudaMalloc(&d_k4, Ny * sizeof(double));

        cudaMemcpy(d_Mat, Mat, Nx*Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Vec, Vec, Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k0, k0, Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k1, k1, Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k2, k1, Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k3, k1, Ny * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_k4, k1, Ny * sizeof(double), cudaMemcpyHostToDevice);

        cudaCheckError();


        int threads = threadsPerBlock;
        int blocks = ceil((Ny + threads - 1) / threads);



        ////////////////////////////////////////////////
        // TIMER START
        ////////////////////////////////////////////////
        auto t_start = std::chrono::high_resolution_clock::now();


        RKstep<<<blocks, threads>>>(d_Mat, d_Vec, Nx, Ny, d_k0, d_k1, 1);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_Mat, d_Vec, Nx, Ny, d_k1, d_k2, 0.5);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_Mat, d_Vec, Nx, Ny, d_k2, d_k3, 0.5);
        cudaDeviceSynchronize();
        RKstep<<<blocks, threads>>>(d_Mat, d_Vec, Nx, Ny, d_k3, d_k4, 1);
        cudaDeviceSynchronize();
        cudaCheckError();


        
        cudaMemcpy(k1, d_k1, Ny * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(k2, d_k2, Ny * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(k3, d_k3, Ny * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(k4, d_k4, Ny * sizeof(double), cudaMemcpyDeviceToHost);
        cudaCheckError();



        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout << std::to_string(std::chrono::duration<double, std::milli>(t_end - t_start).count()) << std::endl;
        file << std::to_string(std::chrono::duration<double, std::milli>(t_end-t_start).count()) << ", ";
        ////////////////////////////////////////////////
        // TIMER END
        ////////////////////////////////////////////////


        cudaFree(d_Mat);
        cudaFree(d_Vec);
        cudaFree(d_k0);
        cudaFree(d_k1);
        cudaFree(d_k2);
        cudaFree(d_k3);
        cudaFree(d_k4);

        free(Mat);
        free(Vec);
        free(k0);
        free(k1);
        free(k2);
        free(k3);
        free(k4);
    }
    file.close();
    return 0;
}