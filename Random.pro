#-------------------------------------------------
#
# Project created by QtCreator 2013-04-28T15:46:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Random
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    describing_statistic.cpp \
    statistical_criteria.cpp \
    gauss_generator.cpp

HEADERS  += mainwindow.h \
    describing_statistic.h \
    statistical_criteria.h \
    gauss_generator.h

FORMS    += mainwindow.ui
