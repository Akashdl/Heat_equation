#ifndef MAPPER2D_H
#define MAPPER2D_H

#include <cstddef>

class Mapper2D
{
public:
    Mapper2D(size_t x, size_t y) : nx(x), ny(x), size(x * y) {}

    size_t sizeValue() const 
    { 
        return size; 
    }

    size_t pos(size_t x, size_t y) const 
    {
         return nx * y + x; 
    }

    size_t xForPos(size_t pos) const 
    {
         return pos % nx;
    }

    size_t yForPos(size_t pos) const 
    {
         return pos / nx; 
    }

    size_t nxValue() const
    {
        return nx;
    }

    size_t nyValue() const
    {
        return ny;
    }
    
private:
    size_t nx;
    size_t ny;
    size_t size;

};


#endif
