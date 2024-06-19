#include <iostream>
#include <cuda_runtime.h>

#define threadsPerBlock 32

__global__ void RKstep(double* Mat, double* Vec, size_t length_x, size_t length_y, double* k_now, double* k_next, double scale_k)
{
    int i = blockIdx.x + blockDim.x + threadIdx.x;
    int j = blockIdx.y + blockDim.y + threadIdx.y;
    
    k_next[j] = k_next[j] + Mat[j * length_x + i] * (Vec[i] + scale_k * k_now[i]);
}


int main ()
{
    size_t length_x = 10;
    size_t length_y = 10;

    double* Mat = (double*)calloc(length_x * length_y, sizeof(double));
    double* Mat_d = NULL;
    cudaMalloc(&Mat_d, length_x * length_y, sizeof(double));

    double* Vec = (double*)calloc(length_y, sizeof(double));
    double* Vec_d = NULL;
    cudaMalloc(&Vec_d, length_y, sizeof(double));

    for (size_t i = 0; i < length_x * length_y; i++)        {Mat[i] = 1.0;  }
    for (size_t i = 0; i < length_y; i++)                   {Vec[i] = i + 1;}

    cudaMemcpy(Mat_d, Mat, length_x * length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Vec_d, Vec, length_x * length_y * sizeof(double), cudaMemcpyHostToDevice);

    double* k0 = (double*)calloc(length_y, sizeof(double));
    double* k1 = (double*)calloc(length_y, sizeof(double));
    double* k2 = (double*)calloc(length_y, sizeof(double));
    double* k3 = (double*)calloc(length_y, sizeof(double));
    double* k4 = (double*)calloc(length_y, sizeof(double));

    double *k0_d = NULL, *k1_d = NULL, k2_d = NULL, k3_d = NULL, k4_d = NULL;

    cudaMalloc(&k0_d, length_y, sizeof(double));
    cudaMalloc(&k1_d, length_y, sizeof(double));
    cudaMalloc(&k2_d, length_y, sizeof(double));
    cudaMalloc(&k3_d, length_y, sizeof(double));
    cudaMalloc(&k4_d, length_y, sizeof(double));

    cudaMemcpy(k0_d, k0, length_y * sizeof(double), cudaMemcpyHostToDevcice);
    cudaMemcpy(k1_d, k1, length_y * sizeof(double), cudaMemcpyHostToDevcice);
    cudaMemcpy(k2_d, k2, length_y * sizeof(double), cudaMemcpyHostToDevcice);
    cudaMemcpy(k3_d, k3, length_y * sizeof(double), cudaMemcpyHostToDevcice);
    cudaMemcpy(k4_d, k4, length_y * sizeof(double), cudaMemcpyHostToDevcice);

    double* result = (double*)calloc(length_y, sizeof(double));

    int threads = threadsPerBlock;
    int blocks_x = (length_x + threads - 1) / threads;
    int blocks_y = (length_y + threads - 1) / threads;

    dim3 blockDim(threads, threads);
    dim3 grifDim(blocks_x, blocks_y);


    
    /*
    //1st step of Runge Kutta method -> Matrix Vector Multiplication
    for (int j = 0; j < length_y; j++)
    {
        for (int i = 0; i < length_x; i++)
        {
            k1[j] = k1[j] + Mat[j * length_x + i] * Vec[i];
        }
    }

    //2nd step of Runge Kutta method -> 1. Vector Assition 2. Matrix Vector Multiplication
    for (int j = 0; j < length_y; j++)
    {
        for (int i = 0; i < length_x; i++)
        {
            k2[j] = k2[j] + Mat[j * length_x + i] * (Vec[i] + 0.5 * k1[i]);
        }
    }
    */
    
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k0, k1, 1.0);
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k1, k2, 0.5);
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k2, k3, 0.5);
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k3, k4, 1.0);
    cudaDeviceSynchronize();

    cudaMemcpy(k0, k0_d, length_y * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(k1, k1_d, length_y * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(k2, k2_d, length_y * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(k3, k3_d, length_y * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(k4, k4_d, length_y * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(Mat_d);
    cudaFree(Vec_d);
    cudaFree(k0_d);
    cudaFree(k1_d);
    cudaFree(k2_d);
    cudaFree(k3_d);
    cudaFree(k4_d);


    //averaging the different steps with appropriate weights
    for (int i = 0; i < length_y; i++) 
    {
        result[i] = Vec[i] + (1.0 / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        std::cout << result[i] << std::endl;
    }

    free(Mat);
    free(Vec);
    free(k0);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(result);
}

