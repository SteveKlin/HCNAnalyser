#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class WidgetHCN;

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    
    WidgetHCN* widgetHCN;

signals:
    
public slots:
    
};

#endif // MAINWINDOW_H
