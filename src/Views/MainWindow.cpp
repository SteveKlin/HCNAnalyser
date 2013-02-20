#include "MainWindow.h"

#include "WidgetHCN.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    widgetHCN = new WidgetHCN(this);

    setCentralWidget(widgetHCN);
}
