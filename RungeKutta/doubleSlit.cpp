#include <iostream>
#include <fstream>
#include <string>
#include <complex>

typedef std::complex<double> complex;


void writeToCSV (double *array, int arraySize, std::string name, std::string filepath)
{
    std::ofstream file;
    file.open(filepath, std::ios::out | std::ios::app);
    file << "\n" << name << ",";

    for (int i = 0; i < arraySize; i++)
    {
        file << abs(array[i]) << ",";
    }
}



double* getProbDensities (complex* array, int arraySize)
{
    double* pd = new double[arraySize];
    for (int i = 0; i < arraySize; i++) {pd[i] = std::norm(array[i]);}

    return pd;
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
void RKstep(complex* Mat, complex* Vec, size_t length_x, size_t length_y, complex* k_now, complex* k_next, double scale_k)
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
    //imaginary unit "i"
    complex I(0.0, 1.0);

    //Parameters
    float L = 10;                                    //well of width L
    double dx = 1.0;                                  //spatial step size
    double dt = pow(dx, 2) / 4;                     //temporal step size
    int Nx = floor(L/dx) + 1;                       //Number of points on the x axis
    int Ny = floor(L/dx) + 1;                       //Number of points on the y axis
    int Nt = 500;                                    //Number of time steps
    complex r = - dt / (2.0 * I * pow(dx, 2));      //Constant to simplify expressions

    //initial position of the center of the gaussian wave packet
    double x0 = L/5;
    double y0 = L/2;

    //momentum and width of the initial wave packet
    double k = 15.0 * M_PI;
    double sigma = 0.5;

    //parameters of the double slit
    double w = 0.2;                                 //width of the double slit
    double s = 0.8;                                 //seperation between the slits
    double a = 0.4;                                 //aperture of the slits
    double V0 = 200.0;                             //constant value of the potential


    //indices that parameterize the double slit
    //horizontal axis
    int i0 = int(1 / (2 * dx) * (L - w));
    int i1 = int(1 / (2 * dx) * (L + w));
    //vertical axis
    int j0 = int(1 / (2 * dx) * (L + s) + a / dx);
    int j1 = int(1 / (2 * dx) * (L + s));
    int j2 = int(1 / (2 * dx) * (L - s));
    int j3 = int(1 / (2 * dx) * (L - s) - a / dx);

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
                exp(I * k * (xarray[i] - x0)) * \
                exp(-1.0 / (2.0 * pow(sigma, 2)) * (pow(xarray[i] - x0, 2) + pow(yarray[j] - y0, 2)));
            }
        }
    }
    double* probDensities = getProbDensities(psi0, Nx*Ny);
    

    //set the potential
    double *V = (double*)malloc(Nx * Ny * sizeof(double));
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i< Nx; i++)
        {
            if ((i <= i0 && i <= i1) && ((j <= j3) || (j > j2 && j <j1) || (j >= j0))) {V[j * Nx + i] = V0;}
            else {V[j * Nx + i] = 0;}
        }
    }


    //write data to output file
    std::ofstream file;
    file.open("../output/output.csv");
    file.clear();
    writeToCSV(xarray, Nx, "X-Array", "../output/output.csv");
    writeToCSV(yarray, Ny, "Y-Array", "../output/output.csv");
    writeToCSV(probDensities, Nx * Ny, "Psi 0", "../output/output.csv");

    /*
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            if (i == Nx - 1)        {std::cout << V[j * Nx + i] << "\n";}
            else                    {std::cout << V[j * Nx + i] << " "; }
        }
    }
    */


    //create the Hamiltonian
    complex *H = (complex*)malloc(Nx * Nx * Ny * Ny * sizeof(complex));
    for (int k = 0; k < Nx * Ny; k++)
    {
        int j = 1 + floor(k/Nx);
        int i = 1 + k % Nx;

        H[k * Nx * Ny + k] = 1.0 + 4.0 * r * I * dt / 2.0 * V[j * Nx + i];

        if (k + 1 <= Nx * Ny) {H[k * Nx * Ny + k + 1] = -r;}
        if (k - 1 >= Nx * Ny) {H[k * Nx * Ny + k - 1] = -r;}
        if (j * Nx + i - 1 <= Nx * Ny) {H[k * Nx * Ny + (j * Nx + i - 1)] = -r;}
        if ((j - 2) * Nx + i - 1 >= 0) {H[k * Nx * Ny + ((j - 2) * Nx + i - 1)] = -r;}
    }





    /*
    for (int j = 0; j < Ny*Nx; j++)
    {
        for (int i = 0; i < Nx*Ny; i++)
        {
            if (i == Nx*Ny - 1)        {std::cout << H[j * Nx*Ny + i] << "\n";}
            else                       {std::cout << H[j * Nx*Ny + i] << " "; }
        }
    }
    */

    //run runge kutta for all time steps
    //std::vector<complex*> psi_t = RungeKutta(H, psi0, Nx, Ny, Nt, dt);
    //delete[] H;
    //delete[] psi0;





}