#ifndef VTKWRITER_H
#define VTKWRITER_H

#include "VtkData.h"

class VtkWriter
{

private:
    VtkData vtkData;
    std::vector< Vec3 > points;
    std::vector< std::vector<int> >  cells;
    size_t nx, ny;

public:
    VtkWriter(const size_t nx, const size_t ny);
    void createPoints();
    // void createCells();
    void createScalarPointData(std::vector<double> &unknownTemperature);
    void createVtkFile(const std::string filename);
};


#endif