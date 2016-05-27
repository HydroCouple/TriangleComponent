#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2016
#License Same as Shewchucks'.
#Triangulation Component

TEMPLATE = lib
TARGET   = TriangleComponent
QT      -= gui

DEFINES += TRIANGLECOMPONENT_LIBRARY


INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include

PRECOMPILED_HEADER = ./include/stdafx.h

HEADERS += ./include/stdafx.h\
           ./include/trianglecomponent_global.h \
           ./include/core/triangle.h \
           ./include/component/trianglecomponent.h

SOURCES += ./src/stdafx.cpp \
           ./src/core/triangle.c \
           ./src/component/trianglecomponent.cpp


CONFIG(debug, debug|release) {

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK.1.0.0
   }
}

CONFIG(release, debug|release) {

    DESTDIR = lib
    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib -lHydroCoupleSDK.1.0.0
   }
}