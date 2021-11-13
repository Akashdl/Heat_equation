#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "SimulationParameter.h"
#include "SolveGaussSeidel.h"
#include "VtkData.h"

int main()
{
    SimulationParameter parameter; // Initialization of Simulation parameters includes Number of grid points and Boundary Conditions.

    size_t gridSizes[3] {65,129,257};
    double centreTemp[3] {0.0, 0.0, 0.0}; // Variable to store Temperature at the centre of the domain for different grid numbers
    double convergence = 0.0;             // Store the convergence value
    double domainLength;                  // Length of the Square domain
    double tol;                           // Store tolerance value for the error.

    std::cout<<"Poisson's equation for Heat Conduction"<<std::endl;               // Project Details
    std::cout<<"Task1: Inhomogeneous Diriclet boundary condition: "<<std::endl;

    std::cout<<"Enter the side of a Square domain: ";
    std::cin>>domainLength;

    while(domainLength < 0)  // Error checking for length.
    { 
        // Throw error if length of Square domain is less than or equal to zero.
        std::cout<<"Invalid side dimensions please provide positive, non-zero value: "<<std::endl;
        std::cin >> domainLength;
    }

    /* Input of tolerance value */
    std::cout << "Enter the tolerance limit in range of 1e-1 to 1e-15: " ; 
    std::cin >> tol; 
    while(tol <= 10e-15 || tol >1) // Error checking for Tolerance variable
    { 
        if(tol <= 10e-15){
            std::cout << "Given tolerance is very small please provide tolerance in range of 1e-1 to 1e-15: " <<std::endl;
            std::cin >> tol;
        }
        else if(tol > 1){
            std::cout << "Given tolerance is very large please provide tolerance in range of 1e-1 to 1e-15: " <<std::endl;
            std::cin >> tol;
        }
    }

    for (size_t i = 0; i<3; i++) // Run the simulation for different number of grid points.
    {
        std::cout<<"Solving for Grid size: "<<gridSizes[i]<<" x "<<gridSizes[i]<<std::endl;
        double dx = domainLength/((double)gridSizes[i] - 1.0);  // step size along x-axis
        double dy = domainLength/((double)gridSizes[i] - 1.0);  // step size along y-axis
        
        double** oldTemp = new double*[gridSizes[i]];   // Declare and Initialize buffer Temperature matrix
        for(size_t j = 0; j < gridSizes[i]; j++)
        {
            oldTemp[j] = new double[gridSizes[i]];
            for(size_t k = 0; k < gridSizes[i]; k++)
            {
                oldTemp[j][k] = 0.0;
            }
        }

        double** newTemp = new double*[gridSizes[i]];   // Declare and Initialize  Temperature matrix
        for(size_t j = 0; j < gridSizes[i]; j++)
        {
            newTemp[j] = new double[gridSizes[i]];
            for(size_t k = 0; k < gridSizes[i]; k++)
            {
                newTemp[j][k] = 0.0;
            }
        }

        /* Assign Boundary Values */
        parameter.setLeftDiricBC(dy, domainLength, newTemp, oldTemp, gridSizes[i]);   // Left BC
        parameter.setRightDiricBC(newTemp, oldTemp, gridSizes[i]);                    // Right BC
        parameter.setTopDiricBC(newTemp, oldTemp, gridSizes[i]);                      // Top Bc
        parameter.setBottomDiricBC(newTemp, oldTemp, gridSizes[i]);                   // Bottom Bc
       
        /* Solving the FD equations*/
        SolveGaussSeidel solve;     // Initialization of Solver class
        double error = 9e9;         // Initialize error to a large number.
        size_t iter = solve.solveIterDirichlet(oldTemp, newTemp, tol, dx, dy, gridSizes[i], error);// Solve and assign number iterations

        centreTemp[i] = newTemp[(gridSizes[i]/2) +1][(gridSizes[i]/2) +1];

        std::cout<<"Finished solving for grid size: "<<gridSizes[i]<<" x "<<gridSizes[i]<<std::endl;
        std::cout<<"Number of iterations to reach tolerance of "<<tol<<" is : "<<iter<<std::endl;

        /*To write VTK file*/
        VtkData vtk;
        
        std::vector<Vec3> points;                   // vector for Points
        std::vector<std::vector<int>> cells;        // vector for cells
        std::vector<double> temperatureField;       // Scalr Temperature field  with zero initialising to contain the Temperatures at nodes

        //Points generation
        for (size_t x = 0; x < gridSizes[i]; x++)
        {
            for (size_t y = 0; y < gridSizes[i]; y++)
            {
                points.push_back(Vec3(x*dx, y*dy, 0.0));  //  Passing the 2D coordinates
            }
        }

        // Cells Generation
        size_t noCells = (gridSizes[i] -1) * (gridSizes[i] -1);       // Calulating number of Elements for Cells generation
        for (size_t ele = 1; ele <= noCells+(gridSizes[i]-2); ele++)  //  Numbering of grid points and forming the Cells for each
        {
            if (ele % gridSizes[i] == 0)        // Depending upon the numbering of Nodes, starting from 0, the Nodes repeat in a pattern
                continue;                       // in every multiple of grid size
            else 
            {
                size_t nodeZero = ele;                              //Corresponds to node zero of the quad Element definition
                size_t nodeOne = nodeZero + gridSizes[i];           //Corresponds to node one of the quad Element definition
                size_t nodeTwo = nodeOne - 1;                       //Corresponds to node two of the quad Element definition
                size_t nodeThree = nodeZero - 1;                    //Corresponds to node three of the quad Element definition
                cells.push_back({(int)nodeZero, (int)nodeOne, (int)nodeTwo, (int)nodeThree});    //Creating the Quadelemnt
            }
        }

        // Scalar Temperature field to hold the Temperatures at nodes
        for(size_t k = 0; k < gridSizes[i]; k++)
        {
            for(size_t j = 0; j < gridSizes[i]; j++)
            {
                temperatureField.push_back(newTemp[k][j]);       //  Passing the Temperature data.
            }
        }


        vtk.setPoints(points);                                  //filling the point data for VTK object
        vtk.setCells(cells);                                    //filling the cell data for VTK object
        vtk.setScalarPointData(temperatureField);               //filling the Scalar value (Temperature) data for VTK object

        std::string size{std::to_string(gridSizes[i])};
        const std::string filename{"Heat_Conduction_DirichletBC_" + size + "x" + size+".vtk"};
        vtk.writeVtkFile(filename);

        /*Delete the temperature matrix*/
        for(size_t k = 0; k < gridSizes[i]; k++)
            delete[] newTemp[k];
        delete[] newTemp;

        /*Delete the buffer temperature matrix*/
        for(size_t k = 0; k < gridSizes[i]; k++)
            delete[] oldTemp[k];
        delete[] oldTemp;
    }

    convergence = - log( (centreTemp[2] - centreTemp[1]) / (centreTemp[1] - centreTemp[0])) / log(2);
    std::cout << "The empirical rate of convergence is: " << convergence << std::endl;
    return 0;
}