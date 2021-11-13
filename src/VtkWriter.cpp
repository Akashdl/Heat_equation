#include "VtkWriter.h"
#include "VtkData.h"

VtkWriter::VtkWriter(const size_t nx, const size_t ny) : nx(nx), ny(ny) {}

void VtkWriter::createPoints()
{
    for (size_t i = 0; i < this->nx; i++)
    {
        for (size_t j = 0; j < this->ny; j++)
        {
            points.push_back( Vec3(j,i,0) );
        }
    }
    vtkData.setPoints(points);
}

void createScalarData(std::vector<double> &unknownTemperature)
{
    vtkData.setScalarPointData(unknownTemperature);
}

// void VtkWriter::createCells()
// {

// }

void VtkWriter::createVtkFile(const std::string filename)
{
    vtkData.writeVtkFile(filename);
}