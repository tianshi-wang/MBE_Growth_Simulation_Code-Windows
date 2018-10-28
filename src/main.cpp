// main.cpp

#include <QApplication>
#include <QDesktopWidget>

#include "window.h"
#include "results.h"
#include "solver.h"


int main(int argc, char *argv[])
{

    QApplication app(argc, argv);
    Window window;
    window.resize(500,935);

    int desktopArea = QApplication::desktop()->width() *
                     QApplication::desktop()->height();
    int widgetArea = window.width() * window.height();

    window.setWindowTitle("MBE Growth Simulation Tool");
    window.show();

//    if (((float)widgetArea / (float)desktopArea) < 0.75f)
//        window.show();
//    else
//        window.showMaximized();
    return app.exec();
}
