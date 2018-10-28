#ifndef SOLVER_H
#define SOLVER_H
#include <QThread>
#include <QString>
#include <QTime>

#include  "results.h"
#include <stdio.h>
#include <math.h>

#define DEPTH 9
#define STEPS 0
#define L (1<<DEPTH)

class Solver : public QThread
{
public:
    Solver(Results * myresult);
    Results * results;
    void run();

//public slots:
//    void setSpeed(double value);

//signals:


};


#endif // SOLVER_H
