#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>




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
    int i = blockIdx.x + blockDim.x + threadIdx.x;
    int j = blockIdx.y + blockDim.y + threadIdx.y;
    
    k_next[j] = k_next[j] + Mat[j * length_x + i] * (Vec[i] + scale_k * k_now[i]);
}


int main ()
{
    size_t length_x = 10;
    size_t length_y = 10;

    double *Mat = (double*)calloc(length_x * length_y, sizeof(double));
    double *Mat_d = NULL;
    cudaMalloc(&Mat_d, length_x * length_y * sizeof(double));
    cudaError_t code_mat = cudaCheckError();

    double *Vec = (double*)calloc(length_y, sizeof(double));
    double *Vec_d = NULL;
    cudaMalloc(&Vec_d, length_y * sizeof(double));
    cudaError_t code_vec = cudaCheckError();

    if (code_mat != cudaSuccess || code_vec != cudaSuccess)
    {
        if (Mat_d != NULL) cudaFree(Mat_d);
        if (Vec_d != NULL) cudaFree(Vec_d);
        free(Mat);
        free(Vec);
        return code_mat;
    }

    for (size_t i = 0; i < length_x * length_y; i++)        {Mat[i] = 1.0;  }
    for (size_t i = 0; i < length_y; i++)                   {Vec[i] = i + 1;}

    cudaMemcpy(Mat_d, Mat, length_x * length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(Vec_d, Vec, length_x * length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaCheckError();

    double *k0 = (double*)calloc(length_y, sizeof(double));
    double *k1 = (double*)calloc(length_y, sizeof(double));
    double *k2 = (double*)calloc(length_y, sizeof(double));
    double *k3 = (double*)calloc(length_y, sizeof(double));
    double *k4 = (double*)calloc(length_y, sizeof(double));

    double *k0_d = NULL, *k1_d = NULL, *k2_d = NULL, *k3_d = NULL, *k4_d = NULL;

    cudaMalloc(&k0_d, length_y * sizeof(double));
    cudaMalloc(&k1_d, length_y * sizeof(double));
    cudaMalloc(&k2_d, length_y * sizeof(double));
    cudaMalloc(&k3_d, length_y * sizeof(double));
    cudaMalloc(&k4_d, length_y * sizeof(double));

    cudaError_t code = cudaCheckError();
    cudaError_t code_1 = cudaCheckError();
    cudaError_t code_2 = cudaCheckError();
    cudaError_t code_3 = cudaCheckError();
    cudaError_t code_4 = cudaCheckError();

    if (code != cudaSuccess || code_1 != cudaSuccess || code_2 != cudaSuccess || code_3 != cudaSuccess || code_4 != cudaSuccess) 
    {
        if (k0_d != NULL) cudaFree(k0_d);
        if (k1_d != NULL) cudaFree(k1_d);
        if (k2_d != NULL) cudaFree(k2_d);
        if (k3_d != NULL) cudaFree(k3_d);
        if (k4_d != NULL) cudaFree(k4_d);
        free(k0);
        free(k1);
        free(k2);
        free(k3);
        free(k4);
        return code;
    }

    cudaMemcpy(k0_d, k0, length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(k1_d, k1, length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(k2_d, k2, length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(k3_d, k3, length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(k4_d, k4, length_y * sizeof(double), cudaMemcpyHostToDevice);
    cudaCheckError();

    double *result = (double*)calloc(length_y, sizeof(double));

    int threads = threadsPerBlock;
    int blocks_x = (length_x + threads - 1) / threads;
    int blocks_y = (length_y + threads - 1) / threads;

    dim3 blockDim(threads, threads);
    dim3 grifDim(blocks_x, blocks_y);


    
    
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k0, k1, 1.0);
    cudaCheckError();
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k1, k2, 0.5);
    cudaCheckError();
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k2, k3, 0.5);
    cudaCheckError();
    cudaDeviceSynchronize();
    RKstep<<<gridDim, blockDim>>>(Mat, Vec, length_x, length_y, k3, k4, 1.0);
    cudaCheckError();
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

