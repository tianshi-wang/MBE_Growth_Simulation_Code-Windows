// window.cpp

#include <QtWidgets>
#include "window.h"
#include "ui_window.h"
#include <Qstring>
#include "myglwidget.h"


Window::Window(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Window)
{


    ui->setupUi(this);

    //Manully add some connections. If the variables in myglwidget are changed, modify the UI display.
    connect(ui->myGLWidget, SIGNAL(xRotationChanged(int)), ui->rotXSlider, SLOT(setValue(int)));
    connect(ui->myGLWidget, SIGNAL(zRotationChanged(int)), ui->rotZSlider, SLOT(setValue(int)));

    connect(ui->myGLWidget, SIGNAL(HeightSignal(double)), ui->label_11, SLOT(setNum(double)));
    connect(ui->myGLWidget, SIGNAL(currentTimeSignal(double)), ui->label_12, SLOT(setNum(double)));
    connect(ui->myGLWidget, SIGNAL(simulationEfficiencySignal(QString)), ui->label_24, SLOT(setText(QString)));

    connect(ui->myGLWidget, SIGNAL(UncertaintySignal(double)), ui->label_22, SLOT(setNum(double)));
    connect(ui->myGLWidget, SIGNAL(gridSizeSignal(int)), ui->label_17, SLOT(setNum(int)));
    //Change both the spinbonx (1st) and slider (2nd).
    connect(ui->myGLWidget, SIGNAL(bondingEnergySignal(int)), ui->horizontalSlider_3, SLOT(setValue(int)));
    connect(ui->myGLWidget, SIGNAL(bondingEnergySignalSlider(double)), ui->doubleSpinBox_2, SLOT(setValue(double)));
    connect(ui->myGLWidget, SIGNAL(hoppingEnergySignal(int)), ui->horizontalSlider_4, SLOT(setValue(int)));
    connect(ui->myGLWidget, SIGNAL(hoppingEnergySignalSlider(double)), ui->doubleSpinBox_3, SLOT(setValue(double)));

    connect(ui->myGLWidget, SIGNAL(desorptionEnergySignal(int)), ui->horizontalSlider_5, SLOT(setValue(int)));
    connect(ui->myGLWidget, SIGNAL(desorptionEnergySignalSlider(double)), ui->doubleSpinBox, SLOT(setValue(double)));
    connect(ui->myGLWidget, SIGNAL(arrivalRateSignal(int)), ui->horizontalSlider, SLOT(setValue(int)));
    connect(ui->myGLWidget, SIGNAL(arrivalRateSignalSlider(double)), ui->doubleSpinBox_4, SLOT(setValue(double)));

    //add tooltips for some labels
//    ui->label_8->setToolTip("ml stands for monolayer");


}

Window::~Window()
{
    delete ui;
}


void Window::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_Escape)
        close();
    else
        QWidget::keyPressEvent(e);
}
