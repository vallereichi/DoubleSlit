#include <iostream>

#include "wavePacket.h"

int main ()
{
    //imaginary unit "i"
    complex I(0.0, 1.0);

    //Parameters
    float L = 3;                                    //well of width L
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
    complex* xarray = createArray(L, dx);
    complex* yarray = createArray(L, dx);

    complex* psi0 = gaussWavePacket2D(xarray, yarray, Nx, Ny, x0, y0, k, sigma);

    delete[] xarray;
    delete[] yarray;

    //for (int i = 0; i < Nx*Ny; i++) std::cout << psi0[i] << std::endl;
    
    double* probDensities = getProbDensities(psi0, Nx*Ny);
    writeToFile(probDensities, Nx*Ny, "../output/psi0.txt", 0);
    delete[] probDensities;

    //set the potential
    double* V = setPotential(Nx, Ny, V0, i0, i1, j0, j1, j2, j3);

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
    complex* H = Hamiltonian(V, Nx, Ny, dx, r, dt);
    //delete[] V;

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
    std::vector<complex*> psi_t = RungeKutta(H, psi0, Nx, Ny, Nt, dt);
    //delete[] H;
    //delete[] psi0;





}