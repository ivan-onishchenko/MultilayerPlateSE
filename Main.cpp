#include "MainWindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    QFile file(":/styleSheets/styleSheets/EasyCode.qss");
    file.open(QFile::ReadOnly);
    QString styleSheet = QLatin1String(file.readAll());
    w.setStyleSheet(styleSheet);
    w.show();
    return a.exec();
}
