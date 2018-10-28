#include "griddata.h"


GridData::GridData(int size)
{
    this->setSize(size);
    heights = new int[size*size];
    for (int i=0; i< size*size; i++)
        heights[i] = 0;
}


int GridData::getAt(int x, int y)
{
    if (x < size && x >= 0 && y < size && y >= 0)
        return heights[x+y*size];
    else
        return 0;
}

void GridData::setAt(int x, int y, int val)
{
    heights[x + y*size] = val;
}
