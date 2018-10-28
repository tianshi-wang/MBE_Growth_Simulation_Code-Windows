#ifndef RESULTS_H
#define RESULTS_H
#include "griddata.h"


class Results
{
public:
    Results();
    double time;  //the growth time
    double timescale;  //the growth time over code execution time updated for each loop
    double bondEnergy, desorptionEnergy, hoppingEnergy, arrivalRate, temperature, avgHeight;
    bool reset,resetTimer, resume, stop;
    int steps;
    double speedratio;
    GridData *data = new GridData(64);
    void initHeight();

};

//- (void)lock;
//- (void)unlock;
//- (void)lockForRun;
//- (void)unlockForRun;



#endif // RESULTS_H

