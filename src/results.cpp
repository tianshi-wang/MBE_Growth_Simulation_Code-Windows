#include "results.h"
#include <iostream>
#include <string>

Results::Results()
{
    temperature = 1300.0;
    time = 0.0;
    timescale = 0.0;
    bondEnergy = 0.5;
    desorptionEnergy = 20.0;
    hoppingEnergy = 2.0;
    arrivalRate = 1;
    reset = false;
    resume = false;
    stop = false;
    resetTimer = false;
    steps = 3;
    speedratio = 0.0;
    printf("Pass");
    avgHeight = 0.0;
}


void Results::initHeight(){
    int i, j;
    for (i=0;i<this->data->getSize();i++) {
        for (j=0; j<this->data->getSize(); j++) {
            this->data->setAt(i,j, i*steps/(this->data->getSize()) ); //steps
        }
    }
}

