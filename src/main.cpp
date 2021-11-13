#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "Mapper2D.h"
#include "GaussSeidel.h"
#include "SimulationParameter.h"
#include "VtkWriter.h"


int main()
{

    SimulationParameter parameter;

    const Mapper2D actualGrid(parameter.getGridNx(), parameter.getGridNy());

    std::cout << "Grid size: " << actualGrid.nxValue() << "x" << actualGrid.nyValue() << std::endl;

    // The matrix A of linear system of equation which is banded.
    const Mapper2D matrixGrid((actualGrid.nxValue() * actualGrid.nyValue()) , (actualGrid.nxValue() * actualGrid.nyValue()));

     std::cout << "Matrix size: " << matrixGrid.nxValue() << "x" << matrixGrid.nyValue() << std::endl;

    // Create a 1D Vector of size x*y
    std::vector<double> matrixA;

    /* initialization including boundary conditions*/
    for (size_t i = 0; i < matrixGrid.sizeValue(); i++)
    {
        matrixA.push_back(0.0);
    }

    ///////////////////////////////////////////////////////////////////////
    //To display Matrix.
    
    // for(size_t y = 0;  y < matrixGrid.nyValue(); y++)
    // {
    //     std::cout<<y<<"  ";
    //     for(size_t x = 0; x < matrixGrid.nxValue(); x++)
    //     {
    //         std::cout<<matrixA[matrixGrid.pos(x,y)]<<"  ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;

    /////////////////////////////////////////////////////////////////////////
    // Construct Matrix A for the finite difference scheme

    for (size_t y = 0; y < matrixGrid.nyValue(); y++)
    {
        for (size_t x = 0; x < matrixGrid.nxValue(); x++ )
        {

            // diagonal elements
            if( x == y)
            {
                matrixA[matrixGrid.pos(x,y)] = -4.0;
            }
            else {
                
                if( x == ( y + 1 ) ) // upper diagonal
                {   
                    matrixA[matrixGrid.pos(x,y)] = 1.0;
                    if ( ( ( y + 1 ) % parameter.getGridNy() ) == 0)
                    {
                        matrixA[matrixGrid.pos(x,y)] = 0.0;
                    }
                }
                // lower diagonal
                if( x == ( y - 1 ) )
                {
                    matrixA[matrixGrid.pos(x,y)] = 1.0;
                    if ( ( ( x + 1 ) % parameter.getGridNx() ) == 0)
                    {
                        matrixA[matrixGrid.pos(x,y)] = 0.0;
                    }
                }
                // upper ones
                if ( x == ( y + parameter.getGridNx() ) )
                {
                    matrixA[matrixGrid.pos(x,y)] = 1.0;
                }
                // lower ones
                if ( x == ( y - parameter.getGridNy() ) )
                {
                    matrixA[matrixGrid.pos(x,y)] = 1.0;
                }

            }
        }
    }

    //////////////////////////////////////////////////////////////////////
    // To display Matrix after constructing it with finite difference scheme
    
    // for(size_t y = 0;  y < matrixGrid.nyValue(); y++)
    // {
    //     std::cout<<y<<"  ";
    //     for(size_t x = 0; x < matrixGrid.nxValue(); x++)
    //     {
    //         std::cout<<matrixA[matrixGrid.pos(x,y)]<<"  ";
    //     }
    //     std::cout<<std::endl;
    // }

    /////////////////////////////////////////////////////////////////////////
    // To display BC

    // for (size_t x = 0; x < parameter.getGridNx(); x++)
    // {
    //     std::cout<<parameter.getBCTop(x)<<"     "<<parameter.getBCBottom(x)<<std::endl;
    // }    
    // for(size_t y = 0; y < parameter.getGridNy(); y++)
    // {
    //     std::cout<<std::endl<<parameter.getBCRight(y)<<"     "<<parameter.getBCLeft(y)<<std::endl;
    // }

    //////////////////////////////////////////////////////////////////////////
    // Calculating RHS

    std::vector<double> rhs;

    for(size_t x = 0; x < parameter.getGridNx()  ; x++)
    {
        for(size_t y = 0; y < parameter.getGridNy() ; y++)
        {
            // Top Left Corrner
            if(x == 0 && y == 0)
                rhs.push_back(-parameter.getBCTop(0) - parameter.getBCLeft(0));

            // Bottom Left Corner
            if ( x == 0 && y == parameter.getGridNy() - 1 )
                rhs.push_back(-parameter.getBCBottom(0) - parameter.getBCLeft(y));

            // Top Right Corner
            if (x == parameter.getGridNx() - 1 && y == 0 )
                rhs.push_back(- parameter.getBCTop(x) - parameter.getBCRight(0));

            // Bottom Right Corner
            if (x == parameter.getGridNx() - 1 && y == parameter.getGridNy() - 1)
                rhs.push_back(- parameter.getBCBottom(x) - parameter.getBCRight(y));

            // Left Intermediate
            if( x == 0 && y > 0 && y < parameter.getGridNy() - 1 )
                rhs.push_back( - parameter.getBCLeft(y) );
                
            // Right Intermediate
            if( x == parameter.getGridNx() - 1 && y > 0 && y < parameter.getGridNy() - 1 )
                rhs.push_back( - parameter.getBCRight(y) );

            // Top Intermediate
            if (x > 0 && x < parameter.getGridNx() - 1 && y == 0 )
                rhs.push_back( - parameter.getBCTop(x) );

            // Bottom Intermediate
            if ( x > 0 && x < parameter.getGridNx() - 1 && y == parameter.getGridNy() - 1 )
                rhs.push_back( - parameter.getBCBottom(x) );

            // Inner
            if( x > 0 && x < parameter.getGridNx() - 1 &&  y > 0 && y < parameter.getGridNy() - 1 )
                rhs.push_back( 0.0 );
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // To display RHS

    // std::cout<<std::endl;
    // for(size_t i = 0; i < parameter.getGridSize(); i++)\
    // {
    //     std::cout<<rhs[i]<<std::endl;
    // }

    ////////////////////////////////////////////////////////////////////////////
    // Solving the Equations

    GaussSeidel solve(matrixGrid.nxValue(), matrixGrid.nyValue());

    // Create and Intialize Unknown Temperature grid points with Initial guess 
    std::vector<double> unknownTemperature((int)parameter.getGridSize(), 0.0);

    solve.solveIter(matrixA, rhs, unknownTemperature);

    ///////////////////////////////////////////////////////////////////////////////
    // To display Calculated Temperatures

    // std::cout<<std::endl;
    // for(size_t i = 0; i < parameter.getGridSize(); i++)
    // {
    //     std::cout<<unknownTemperature[i]<<std::endl;
    // }

    /////////////////////////////////////////////////////////////////////////////////
    // To display temperature on a grid

    // for( size_t x = 0; x < parameter.getGridNx(); x++ )
    // {
    //     for( size_t y = 0; y < parameter.getGridNy(); y++ )
    //     {
    //         std::cout<< unknownTemperature[actualGrid.pos(x,y)]<<std::setw(10);
    //     }
    //     std::cout<<std::endl;
    // }

    /////////////////////////////////////////////////////////////////////////////////
    // To write VTK file

    VtkWriter vtkfile(parameter.getGridNx(), parameter.getGridNy());
    vtkfile.createPoints();
    vtkfile.createScalarPointData(unknownTemperature);
    std::string calcName{"Heat_Conduction_DirichletBC_"};
    std::string size{std::to_string(parameter.getGridNx())};
    std::string filename{calcName + size + "x" + size};
    vtkfile.createVtkFile(filename);


    return 0;
}