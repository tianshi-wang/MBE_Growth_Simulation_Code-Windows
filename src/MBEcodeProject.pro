#-------------------------------------------------
#
# Project created by Tianshi Wang, Nov. 2017
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MBE_Simulation
TEMPLATE = app
RC_ICONS = Icon.ico

SOURCES += main.cpp\
        window.cpp \
    myglwidget.cpp \
    griddata.cpp \
    results.cpp \
    solver.cpp

HEADERS  += window.h \
    myglwidget.h \
    griddata.h \
    results.h \
    solver.h

FORMS    += window.ui

#For Windows user: change the route to your library position to link OpenGL and glu.
LIBS += -LC:\Qt\5.6.3\mingw49_32\lib\libQt5OpenGL.a -lopengl32
LIBS += -LC:\Qt\5.6.3\mingw49_32\lib -lglu32
#For Mac, comment the two lines above and uncomment the following 4 lines. Also change "#include "GL/glu.h"" to "#include "OpenGL/glu.h"" in myglwidget.cpp.
#mac: LIBS += -framework OpenGL
#mac: LIBS += -framework GLUT

