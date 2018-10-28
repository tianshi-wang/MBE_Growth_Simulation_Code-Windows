// myglwidget.h

#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QGLWidget>
#include "results.h"
#include "solver.h"
#include "QThread"
#include <QLabel>
#include <Qstring>
#include <string>


class MyGLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit MyGLWidget(QWidget *parent = 0);
    ~MyGLWidget();
    QTimer *timer  ;
    GridData * newdata;
    void emitHeightSignal (double height){
        emit HeightSignal((int)(height/0.01)*0.01);
    }
    void emitCurrentTimeSignal(double currentTime){
        emit currentTimeSignal((int)(currentTime/0.01)*0.01);
    }
    void emitSimulationEfficiencySignal(QString simulationEfficiencyStr){
        emit simulationEfficiencySignal(simulationEfficiencyStr); //in percentage
    }

    void emitUncertaintySignal(double uncertainty){
        emit UncertaintySignal((int)(uncertainty/0.01)*0.01);
    }


protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

public slots:
    // slots for xyz-rotation slider

    void setXRotation(int angle);
//    void setYRotation(int angle);
    void setZRotation(int angle);
    void setSpeed(int val);

    void setDesorptionEnergy(double desorptionEnergy);
    void setDesorptionEnergyfromSlider(int desorptionEnergySlider);


    void setBondingEnergy(double bondingEnergy);
    void setBondingEnergyfromSlider(int bondingEnergySlider);

    void setHoppingEnergy(double hoppingEnergy);
    void setHoppingEnergyfromSlider(int hoppingEnergySlider);

    void setArrivalRate(double arrivalRate);
    void setArrivalRateSlider(int arrivalRateSlider);

    void setTemperature(int temperature);
    void setResumeSignal ();
    void setStopSignal ();
    void setResetSignal();
    void setGridSize(int size);
    void setSteps(int steps);
    void popOutAbout();


signals:
    // signaling rotation from mouse movement
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
    void HeightSignal(double height);
    void currentTimeSignal (double crrentTime);
    void simulationEfficiencySignal (QString simulationEfficiencyStr);

    void UncertaintySignal(double uncertainty);
    void gridSizeSignal(int i);
    void bondingEnergySignal(int i);
    void bondingEnergySignalSlider(double bondingEnergySlider);
    void hoppingEnergySignal(int i);
    void hoppingEnergySignalSlider(double desorptionEnergySlider);
    void desorptionEnergySignal(int i);
    void desorptionEnergySignalSlider(double desorptionEnergySlider);
    void arrivalRateSignal(int arrivalRate);
    void arrivalRateSignalSlider(double arrivalRateSlider);


private:
    void draw();
    int xRot;
    int yRot;
    int zRot;
    double simulationEfficiency;
    QString simulationEfficiencyStr;
    double currentTime, runTime;
    bool reset, resume, stop;
    QPoint lastPos;
};

#endif // MYGLWIDGET_H

