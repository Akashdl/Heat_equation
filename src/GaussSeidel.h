#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include <vector>
#include <iostream>
#include "Mapper2D.h"

class GaussSeidel
{
private:
    std::vector<double> initialGuess;
    size_t gridNx, gridNy, convergenceInterval;
    const double tol = 1.0e-6;
    

public:
    
    GaussSeidel(const size_t gridNx, const size_t gridNy) : gridNx(gridNx) , gridNy(gridNy)
    {
        for(size_t i = 0; i < gridNy ; i++)
            initialGuess.push_back(0.0);
    }

    void solveIter(std::vector<double> &matrix, std::vector<double> &rhs, std::vector<double> &unknownTemperature)
    {
        const Mapper2D matrixGrid(gridNx, gridNy);
        int m = 20;
        int flag = 0;
        std::vector<double> tmp((int)matrixGrid.nyValue(), 0.0);
        while(m>0)
        {
            for (int i = 0; i < matrixGrid.nyValue(); i++)
            {
                tmp[i] = (rhs[i] / matrix[matrixGrid.pos(i,i)]);

                for (int j = 0; j < matrixGrid.nyValue(); j++)
                {
                    if (j == i)
                        continue;
                    tmp[i] = tmp[i] - ((matrix[matrixGrid.pos(j,i)] / matrix[matrixGrid.pos(i,i)] ) * unknownTemperature[j]);
                }
                // if( abs ( unknownTemperature[i] - tmp[i] ) <= tol)
                //     flag++;
                unknownTemperature[i] = tmp[i];
            }

            m--;

        } //while(flag < matrixGrid.nyValue());

    }

};


#endif