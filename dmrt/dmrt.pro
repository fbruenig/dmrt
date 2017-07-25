#TEMPLATE = app
TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
CONFIG += libc++

QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -pedantic -fPIC -shared -fopenmp

INCLUDEPATH +=  /home/rottee/Python/include/python2.7
#INCLUDEPATH +=  /usr/include/x86_64-linux-gnu/c++/4.9

LIBS += -L/home/rottee/Python/lib/python2.7 -lpython2.7  -lgomp

SOURCES += dmrtreader.cpp \
#    dmrtpy.cpp \
#    main.cpp \
    dmrtmain.cpp \
    libdmrt.cpp \
    dmrtalg2.cpp


include(deployment.pri)
qtcAddDeployment()

HEADERS += Python.h\
    dmrtreader.h \
#    dmrtpy.h \
    dmrtmain.h \
    dmrtalg2.h


