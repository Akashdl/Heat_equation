#include "VtkData.h"

int main()
{
    VtkData vtkData;

    std::vector< Vec3 > points;

    points.push_back( Vec3( 0,0,1));
    points.push_back( Vec3( 1,0,1));
    points.push_back( Vec3( 1,1,0));
    points.push_back( Vec3( 0,1,0));

    std::vector< std::vector<int> >  cells;

    cells.push_back( {0,2} );
    cells.push_back( {0,1,3} );
    cells.push_back({0,1,2,3});

    std::vector< double > scalarPointData(4);
    std::vector< Vec3 > vectorPointData(4);

    vtkData.setPoints(points);
    vtkData.setCells(cells);
    vtkData.setScalarPointData(scalarPointData);
    vtkData.setVectorPointData(vectorPointData);

    vtkData.writeVtkFile("myFile.vtk");

    return 0;
}