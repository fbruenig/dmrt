TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -pedantic -fopenmp

LIBS += -lgomp

SOURCES += main.cpp \
    ../dmrt/dmrtalg.cpp \
    ../dmrt/dmrtmain.cpp \
    ../dmrt/dmrtreader.cpp \
    ../dmrt/dmrtalg2.cpp

include(deployment.pri)
qtcAddDeployment()

SUBDIRS += \
    ../dmrt/dmrt.pro

HEADERS += \
    ../dmrt/dmrtalg.h \
    ../dmrt/dmrtmain.h \
    ../dmrt/dmrtreader.h \
    ../dmrt/dmrtalg2.h

