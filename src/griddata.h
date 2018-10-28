#ifndef GRIDDATA_H
#define GRIDDATA_H


class GridData
{
public:
    GridData(int size);
    int getSize() { return size; }
    void setSize(int size) { this->size = size; }
    void setAt(int x, int y, int val);
    int getAt(int x, int y);
    int *getHeights() { return heights; }
protected:
    int *heights;
    int size;
};

#endif // GRIDDATA_H

