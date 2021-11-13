#ifndef SIMULATIONPARAMETER_H
#define SIMULATIONPARAMETER_H

#include <vector>

class SimulationParameter
{
private:


public:

    SimulationParameter() {}
    // Method for setting Dirichlet BC on left side of the domain.
    void setLeftDiricBC(const double &dy ,const double &h, double **&newTemp, double **&oldTemp, const size_t &gridsize)
    {
        for (size_t y = 0; y < gridsize ; y++)
        {
            newTemp[y][0] = ( y*dy * (h - y*dy) ) / 4;
            oldTemp[y][0] = ( y*dy * (h - y*dy) ) / 4;
        }
    }
    // Method for setting Dirichlet BC on right side of the domain.
    void setRightDiricBC(double **&newTemp, double **&oldTemp, const size_t &gridsize)
    {
        for (size_t y = 0; y < gridsize; y++)
        {
            newTemp[y][(int)gridsize-1] = 0.0;
            oldTemp[y][(int)gridsize-1] = 0.0;
        }
    }
    // Method for setting Dirichlet BC on top side of the domain.
    void setTopDiricBC(double **&newTemp, double **&oldTemp, const size_t &gridsize)
    {
        for (size_t x = 0; x < gridsize; x++)
        {
            newTemp[0][x] = 0.0;
            oldTemp[0][x] = 0.0;
        }
    }
    // Method for setting Dirichlet BC on bottom side of the domain.
    void setBottomDiricBC(double **&newTemp, double **&oldTemp, const size_t &gridsize)
    {
        for (size_t x = 0; x < gridsize; x++)
        {
            newTemp[gridsize-1][x] = 0.0;
            oldTemp[gridsize-1][x] = 0.0;
        }
    }
    // Method for setting Dirichlet BC on left side of the domain for test case 2.
    void setLeftNeumannBC(const double &dy , double **&newTemp, double **&oldTemp, const size_t &gridsize)
    {
        for (size_t i = (size_t)(0.1*gridsize); i < (int)(0.5*gridsize); i++)
        {
            newTemp[i][0] = 1.0;
            oldTemp[i][0] = 1.0;
        }
    }

};

#endif