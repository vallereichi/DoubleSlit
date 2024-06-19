#include <iostream>
#include <string>
#include <fstream>
#include <chrono>

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
void RKstep(double* Mat, double* Vec, size_t length_x, size_t length_y, double* k_now, double* k_next, double scale_k)
{
    for (int j = 0; j < length_y; j++)
    {
        for (int i = 0; i < length_x; i++)
        {
            k_next[j] = k_next[j] + Mat[j * length_x + i] * (Vec[i] + scale_k * k_now[i]);
        }
    }
}


int main ()
{
    std::ofstream file;
    file.open("../timecpp.txt", std::ios::out | std::ios::app);

    for (int t = 0; t < 100; t++)
    {
        size_t length_x = 100*t;
        size_t length_y = 100*t;

        double* Mat = (double*)calloc(length_x * length_y, sizeof(double));

        double* Vec = (double*)calloc(length_y, sizeof(double));

        double* k0 = (double*)calloc(length_y, sizeof(double));
        double* k1 = (double*)calloc(length_y, sizeof(double));
        double* k2 = (double*)calloc(length_y, sizeof(double));
        double* k3 = (double*)calloc(length_y, sizeof(double));
        double* k4 = (double*)calloc(length_y, sizeof(double));

        double* result = (double*)calloc(length_y, sizeof(double));


        for (size_t i = 0; i < length_x * length_y; i++)        {Mat[i] = 1.0;  }
        for (size_t i = 0; i < length_y; i++)                   {Vec[i] = i + 1;}

        auto t_start = std::chrono::high_resolution_clock::now();

        RKstep(Mat, Vec, length_x, length_y, k0, k1, 1.0);
        RKstep(Mat, Vec, length_x, length_y, k1, k2, 0.5);
        RKstep(Mat, Vec, length_x, length_y, k2, k3, 0.5);
        RKstep(Mat, Vec, length_x, length_y, k3, k4, 1.0);

        auto t_end = std::chrono::high_resolution_clock::now();

        //averaging the different steps with appropriate weights
        for (int i = 0; i < length_y; i++) 
        {
            result[i] = Vec[i] + (1.0 / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            std::cout << result[i] << std::endl;
        }


        
        file << std::to_string(std::chrono::duration<double, std::milli>(t_end-t_start).count()) << ", ";
        std::cout << std::to_string(std::chrono::duration<double, std::milli>(t_end-t_start).count()) << std::endl;
        


        free(Mat);
        free(Vec);
        free(k0);
        free(k1);
        free(k2);
        free(k3);
        free(k4);
        free(result);
    }
    file.close();
}   

