# Connect the necessary Qt modules.
QT += core gui datavisualization widgets printsupport

# Choosing the C++ standard.
CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

RC_ICONS = MultilayerPlateSE.ico

HEADERS += \
    MainWindow.h  \
    MatLibForm.h  \
    qcustomplot.h

SOURCES += \
    Main.cpp        \
    MainWindow.cpp  \
    MatLibForm.cpp  \
    qcustomplot.cpp

FORMS += \
    MainWindow.ui \
    MatLibForm.ui

RESOURCES += \
    Resources.qrc

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
