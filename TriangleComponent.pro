#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2016
#License Same as Shewchucks'.
#Triangulation Component


DEFINES += TRIANGLECOMPONENT_LIBRARY

contains(DEFINES,TRIANGLECOMPONENT_LIBRARY){
  TEMPLATE = lib
  message("Compiling as library")
} else {
  TEMPLATE = app
  CONFIG-=app_bundle
  message("Compiling as application")
}

VERSION = 1.0.0.0
TARGET   = TriangleComponent
QT      += core concurrent

INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include

PRECOMPILED_HEADER = ./include/stdafx.h

HEADERS += ./include/stdafx.h\
           ./include/trianglecomponent_global.h \
           ./include/core/triangle.h \
           ./include/component/trianglecomponent.h \
           ./include/component/trianglecomponentinfo.h

SOURCES += ./src/stdafx.cpp \
           ./src/core/triangle.c \
           ./src/core/tricall.c \
           ./src/component/trianglecomponent.cpp \
           ./src/component/trianglecomponentinfo.cpp \
           ./src/main.cpp


macx{
INCLUDEPATH += /usr/local \
               /usr/local/include

}


CONFIG(debug, debug|release) {

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK
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
    LIBS += -L./../HydroCoupleSDK/lib -lHydroCoupleSDK
   }
}
