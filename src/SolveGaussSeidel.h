#ifndef SOLVEGAUSSSEIDEL_H
#define SOLVEGAUSSSEIDEL_H

#include <cmath> 
#include <vector>

class SolveGaussSeidel
{
private:    

public:
    
    SolveGaussSeidel() {}

    size_t solveIterDirichlet(double **&oldTemp, 
                                double **&newTemp, 
                                const double &tol, 
                                const double &dx, 
                                const double &dy, 
                                const size_t &gridsize, 
                                double &error)
    {
        double a = pow(dx,2);   
        double b = pow(dy,2);
        double k = 2*((1/a)+(1/b)); //Truncation of node in FD Formula in 2D
        size_t iter = 0;

        while( error > tol )
        {
            /* FD Discritization scheme */
            for (size_t i = 1; i < gridsize-1; i++ )
            {
                for (size_t j = 1; j < gridsize-1; j++ )
                {
                    //Gauß Siedel FD fomula
                    newTemp[i][j] = ( ( ( ( newTemp[i-1][j] + oldTemp[i+1][j] ) / a ) + ( ( newTemp[i][j-1] + oldTemp[i][j+1] ) / b )) / k);
                }
            }

            error  = fabs(oldTemp[(gridsize/2) +1][2] - newTemp[(gridsize/2) +1][2]);

            if(iter<100)    // Making sure loop does not end for small iterations.
                error = 10;
            
            /* Updating the new temperature values in buffer matrix */
            for (size_t i = 1; i < gridsize-1; i++ )
            {
                for (size_t j = 1; j < gridsize-1; j++ )
                {
                    oldTemp[i][j] = newTemp[i][j];
                }
            }

            iter++;
        }
        return iter;
    }


    size_t solveIterNeumann(double **&oldTemp, 
                                double **&newTemp, 
                                const double &tol, 
                                const double &dx, 
                                const double &dy, 
                                const size_t &gridsize, 
                                double &error,
                                const double &domainLength)
    {
        double a = pow(dx,2);   
        double b = pow(dy,2);
        double k = 2*((1/a)+(1/b)); //Truncation of node in FD Formula in 2D
        size_t iter = 0;

        double** buf = new double*[gridsize];   
        for(size_t j = 0; j < gridsize; j++)
        {
            buf[j] = new double[gridsize];
            for(size_t k = 0; k < gridsize; k++)
            {
                buf[j][k] = sqrt( pow( k*dy, 2 ) + pow( (domainLength - (j*dx)), 2 ) );  //Distance of first node from the last node of 1st column
            }
        }

        while(error>tol)
        {
            for(size_t i = 1; i < gridsize-1; i++)
            {
                for(size_t j =1; j < gridsize-1; j++)
                {
                    if (buf[i][j] >= 0.5*domainLength && buf[i][j] <= 0.9*domainLength)
                    {
                        if(buf[i+1][j] <= 0.5*domainLength) newTemp[i+1][j] = newTemp[i][j]; //Applying insulation boundary condition
                        if(buf[i-1][j] >= 0.9*domainLength) newTemp[i-1][j] = newTemp[i][j]; //Applying insulation boundary condition
                        if(buf[i][j+1] >= 0.9*domainLength) newTemp[i][j+1] = newTemp[i][j]; //Applying insulation boundary condition
                        if(buf[i][j-1] <= 0.5*domainLength) newTemp[i][j-1] = newTemp[i][j]; //Applying insulation boundary condition

                        //Gauß Siedel FD fomula
                        newTemp[i][j] = ( ( ( ( newTemp[i-1][j] + oldTemp[i+1][j] ) / a ) + ( ( newTemp[i][j-1] + oldTemp[i][j+1] ) / b )) / k);
                    }
                }
            }

            error  = fabs(oldTemp[(gridsize/2) +1][2] - newTemp[(gridsize/2) +1][2]);
            
            if(iter<100)  // Making sure loop does not end for small iterations.
                error = 10;
            
            /* Updating the new temperature values in buffer matrix */
            for (size_t i = 1; i < gridsize-1; i++ )
            {
                for (size_t j = 1; j < gridsize-1; j++ )
                {
                    oldTemp[i][j] = newTemp[i][j];
                }
            }

            iter++;
        }

        for(size_t k = 0; k < gridsize; k++)
            delete[] buf[k];
        delete[] buf;
        
        return iter;

    }

};

#endif