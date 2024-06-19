#include <iostream>

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
    size_t length_x = 10;
    size_t length_y = 10;

    double* Mat = (double*)calloc(length_x * length_y, sizeof(double));
    if (Mat == NULL)
    {
        std::cout << "Not enough RAM for an array of size " << length_x << " x " << length_y << std::endl;
        return 1;
    }

    double* Vec = (double*)calloc(length_y, sizeof(double));

    double* k0 = (double*)calloc(length_y, sizeof(double));
    double* k1 = (double*)calloc(length_y, sizeof(double));
    double* k2 = (double*)calloc(length_y, sizeof(double));
    double* k3 = (double*)calloc(length_y, sizeof(double));
    double* k4 = (double*)calloc(length_y, sizeof(double));

    double* result = (double*)calloc(length_y, sizeof(double));


    for (size_t i = 0; i < length_x * length_y; i++)        {Mat[i] = 1.0;  }
    for (size_t i = 0; i < length_y; i++)                   {Vec[i] = i + 1;}

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
    
    RKstep(Mat, Vec, length_x, length_y, k0, k1, 1.0);
    RKstep(Mat, Vec, length_x, length_y, k1, k2, 0.5);
    RKstep(Mat, Vec, length_x, length_y, k2, k3, 0.5);
    RKstep(Mat, Vec, length_x, length_y, k3, k4, 1.0);

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

