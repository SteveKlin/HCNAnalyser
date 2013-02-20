QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


HEADERS += \
    src/Views/MainWindow.h \
    ../ms/ms.h \
    src/MSSampleGenerator.h \
    src/ms-wrapper/sample_stats.h \
    src/ms-wrapper/mswrapper.h \
    src/HCNGenerator.h \
    src/Views/WidgetHCN.h \
    src/Defines.h \
    src/stringTools.h

SOURCES += \
    src/Views/MainWindow.cpp \
    src/main.cpp \
    ../ms/tajd.c \
    ../ms/streec.c \
    ../ms/rand1t.c \
    ../ms/ms.c \
    src/MSSampleGenerator.cpp \
    src/HCNGenerator.cpp \
    src/ms-wrapper/mswrapper.cpp \
    src/Views/WidgetHCN.cpp \
    src/stringTools.cpp
