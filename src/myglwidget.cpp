// myglwidget.cpp

#include <QtWidgets>
#include <QtOpenGL>
#include <QThread>
#include <thread>
#include "GL/glu.h"
#include <Qstring>
#include <math.h>
#include "griddata.h"
#include <Qstring>
#include <string>
#include <iostream>




#include "myglwidget.h"
#define max(x,y) (x)>(y)?(x):(y)
#define PI 3.14259

QString str = QString::number(123.4567890123456, 'e', 2);

Results *myresult = new Results();

Solver mysolver(myresult);
int size = myresult->data->getSize();
int newheight[L][L];

MyGLWidget::MyGLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    xRot = 315;
    yRot = 0;
    zRot = 45;
}

MyGLWidget::~MyGLWidget()
{
}

QSize MyGLWidget::minimumSizeHint() const
{
    return QSize(300, 300);
}

QSize MyGLWidget::sizeHint() const
{
    return QSize(520, 520);
}

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360;
    while (angle > 360)
        angle -= 360;
}

void MyGLWidget::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}


void MyGLWidget::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void MyGLWidget::setSpeed(int speed)
{
        myresult->speedratio = log10(speed);
        updateGL();
}

void  MyGLWidget::setDesorptionEnergy(double desorptionEnergy){
    myresult->desorptionEnergy = desorptionEnergy;
    emit desorptionEnergySignal((int)(desorptionEnergy*100));
    updateGL();
}
void MyGLWidget::setDesorptionEnergyfromSlider(int desorptionEnergySlider){
    myresult->desorptionEnergy = desorptionEnergySlider*0.01;
    emit desorptionEnergySignalSlider((double)desorptionEnergySlider*0.01);
    updateGL();
}

void  MyGLWidget::setBondingEnergy(double bondingEnergy){
    myresult->bondEnergy = bondingEnergy;
    emit bondingEnergySignal((int)(bondingEnergy*100));
    updateGL();
}

void MyGLWidget::setBondingEnergyfromSlider(int bondingEnergySlider){
    myresult->bondEnergy = bondingEnergySlider*0.01;
    emit bondingEnergySignalSlider((double)bondingEnergySlider*0.01);
    updateGL();
}
void  MyGLWidget::setHoppingEnergy(double hoppingEnergy){
    myresult->hoppingEnergy = hoppingEnergy;
    emit hoppingEnergySignal((int)(hoppingEnergy*100));
    updateGL();
}

void MyGLWidget::setHoppingEnergyfromSlider(int hoppingEnergySlider){
    myresult->hoppingEnergy = hoppingEnergySlider*0.01;
    emit hoppingEnergySignalSlider((double)hoppingEnergySlider*0.01);
    updateGL();
}

void  MyGLWidget::setTemperature(int temperature){
    myresult->temperature = (double)temperature;
    updateGL();
}
void  MyGLWidget::setArrivalRate(double arrivalRate){
    myresult->arrivalRate = (double)arrivalRate;
    emit arrivalRateSignal((int)arrivalRate*100);
    updateGL();
}
void MyGLWidget::setArrivalRateSlider(int arrivalRateSlider){
    myresult->arrivalRate = arrivalRateSlider*0.01;
    emit arrivalRateSignalSlider((double)arrivalRateSlider*0.01);
    updateGL();
}

void MyGLWidget::setGridSize(int size){
     size=1<<size;
     emit gridSizeSignal(size);
     myresult->data = new GridData(size);
     myresult->reset = true;
}

void MyGLWidget::setResumeSignal(){
    myresult->resume = true;
    myresult->stop = false; 

}
void MyGLWidget::setStopSignal(){
    myresult->stop = true;
    myresult->resume = false;
}
void MyGLWidget::setResetSignal(){
    myresult->reset = true;
}
void MyGLWidget::setSteps(int steps){
    myresult->steps = steps;
}
void MyGLWidget::popOutAbout(){
    QMessageBox aboutBox;
    aboutBox.setWindowTitle("About the MBE growth simulation code");
    aboutBox.setStyleSheet("QLabel{min-width: 700px;}");
    aboutBox.setText("The MBE growth simulation code can simulate growth morphology at different conditions for molecular beam epitaxy (MBE) using the kinetic Monte Carlo method. \n"
                     "The code is based on KMCInterative code written by Michael Grundmann (mgrundmann@ece.ucsb.edu). \n"
                     "Developer: Tianshi Wang (tswang@udel.edu) & Wei Li (verali@udel.edu). \n"
                     "Souce files and releases can be found at https://github.com/tianshi-wang/ for Mac and Windows.");
    aboutBox.exec();
}



void setTemperature(double temperature);
void setFlux(double flux);

void MyGLWidget::initializeGL()
{
    qglClearColor(Qt::white);
    myresult->initHeight();
    QTimer *timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(update()));
    timer->start(100);

//    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_CULL_FACE);
//    glShadeModel(GL_SMOOTH);
//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);

//    static GLfloat lightPosition[4] = { 0, 0, 10, 1.0 };
//    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

void MyGLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);
    draw();

    //glRotatef(-phi, 1.0f, 0.0f, 0.0f);
    //glRotatef(theta, 0.0f, 0.0f, 1.0f);

}

void MyGLWidget::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0);

//#ifdef QT_OPENGL_ES_1
//    glOrthof(-2, +2, -2, +2, 1.0, 15.0);
//#else
//    glOrtho(-2, +2, -2, +2, 1.0, 15.0);
//#endif
    glMatrixMode(GL_MODELVIEW);
}

void MyGLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void MyGLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        if (dy>0) setXRotation( (xRot +  dy)<350?(xRot +  dy):350);
        if (dy<0) setXRotation( (xRot +  dy)>285?(xRot +  dy):285);
        if (dx>0) setZRotation( (zRot +  dx)<90?(zRot +  dx):90);
        if (dx<0) setZRotation( (zRot +  dx)>0?(zRot +  dx):0);
    }
    else if (event->buttons() & Qt::RightButton) {
        if (dy>0) setXRotation( (xRot +  dy)<350?(xRot +  dy):350);
        if (dy<0) setXRotation( (xRot +  dy)>285?(xRot +  dy):285);
        if (dx>0) setZRotation( (zRot +  dx)<90?(zRot +  dx):90);
        if (dx<0) setZRotation( (zRot +  dx)>0?(zRot +  dx):0);
    }

    lastPos = event->pos();
}


void MyGLWidget::draw()
{
    if (myresult->reset==true){
        myresult->initHeight();
        mysolver.start();
    }

    int i,j, size, height, nHeight, minHeight, maxHeight,heightRange;
    double avgHeight=0.0;
    double aspect;
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(20.0, 1, 4.0, 15.0);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);

    glDepthMask(GL_TRUE);
    glLoadIdentity();

    glTranslatef(0.0f, 0.0f, -7.0f);
    glRotatef(xRot, 1.0, 0.0, 0.0);
//    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
//      glRotatef(-45, 1.0, 0.0, 0.0);

    glRotatef(zRot, 0.0, 0.0, 1.0);


//    glMatrixMode(GL_MODELVIEW);
//    glEnable(GL_DEPTH_TEST);
//    glDepthMask(GL_TRUE);
//    glLoadIdentity();

    minHeight = 1<<20;
    maxHeight = -minHeight;



    if (true) {
        size = myresult->data->getSize();

        for (i=0;i < size; i++)
            for (j=0; j < size; j++){
               newheight[i][j] = myresult->data->getAt(i,j);

            }


        for (i=0;i < size; i++) {
            for (j=0; j < size; j++) {
                height = newheight[i][j];

                avgHeight+=height;
                if (height<minHeight)
                    minHeight = height;
                if (height>maxHeight)
                    maxHeight = height;
            }
        }


        size =  myresult->data->getSize();
        currentTime = myresult->time;
        reset = myresult->reset;
        heightRange= max(maxHeight-minHeight,5);
//        printf("HeightRange is %d \n", heightRange);
        //heightRange=10;
        avgHeight/=(size*size);

        double uncertaintySquare;
        for (i=0;i < size; i++) {
            for (j=0; j < size; j++) {
                uncertaintySquare += pow(((double)newheight[i][j]-avgHeight),2.0);
            }
         }
        double uncertainty = sqrt(uncertaintySquare)/(double)size;
//        if(runTime<0.001)
//            simulationEfficiency = 0;
//        else simulationEfficiency = currentTime/runTime;

        emitHeightSignal(avgHeight);
        emitCurrentTimeSignal(currentTime);
        simulationEfficiencyStr = QString::number(myresult->timescale, 'e', 0 );
        emitSimulationEfficiencySignal(simulationEfficiencyStr);
        emitUncertaintySignal(uncertainty);

//        std::string varAsString = std::to_string(simulationEfficiency);
//        QString str1 = QString::number( simulationEfficiency, 'g', 6 );




//        printf("avgHeight is %f \n", avgHeight);

        glTranslatef(0.0f, 0.0f, -avgHeight/(double)size);  //define avgHeight=0
        glBegin(GL_QUADS);


        for (i=0;i < size; i++) {
            for (j=0; j < size; j++) {
                height = newheight[i][j];


                if (heightRange !=0)
                    glColor3f(max(cos(PI*((float)((float)height-avgHeight)/(float)(heightRange)-0.5)),0.0),
                        max(cos(PI*(((float)height-avgHeight)/(float)(heightRange))),0.0),
                        max(cos(PI*(((float)height-avgHeight)/(float)(heightRange)+0.5)),0.0));
                glVertex3f((i-size/2)/(size/2.0), (j-size/2)/(size/2.0), height/(float)size);
                glVertex3f((i-size/2+1)/(size/2.0), (j-size/2)/(size/2.0), height/(float)size);
                glVertex3f((i-size/2+1)/(size/2.0), (j-size/2+1)/(size/2.0), height/(float)size);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2+1)/(size/2.0), height/(float)size);
             }
            }


        for (i=0;i <= size; i++) {
            for (j=0; j <= size; j++) {
                height = newheight[i%size][j%size];
                if (j!=size) {
                nHeight = newheight[(i-1+size)%size][j%size];
                glColor3f(0.25f, 0.25f, 0.25f);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2)/(size/2.0), height/(float)size);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2)/(size/2.0), nHeight/(float)size);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2+1)/(size/2.0), nHeight/(float)size);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2+1)/(size/2.0), height/(float)size);
                }
                if (i!=size) {
                nHeight = newheight[i%size][(j-1+size)%size];
                glColor3f(0.25f, 0.25f, 0.25f);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2)/(size/2.0), height/(float)size);
                glVertex3f((i-size/2)/(size/2.0), (j-size/2)/(size/2.0), nHeight/(float)size);
                glVertex3f((i-size/2+1)/(size/2.0), (j-size/2)/(size/2.0), nHeight/(float)size);
                glVertex3f((i-size/2+1)/(size/2.0), (j-size/2)/(size/2.0), height/(float)size);
                }

            }
        }



        glEnd();
    }

}
