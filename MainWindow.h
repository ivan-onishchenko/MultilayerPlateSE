#pragma once

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define M_PI (3.14159265358979323846)
#define M_E (2.7182818284590451)

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <QMainWindow>
#include <QTabBar>
#include <QDebug>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurface3DSeries>

const double c = 299792458;
const double eps0 = (1 / (4 * M_PI * c * c)) * pow(10,7);
const double mu0 = 4 * M_PI * pow(10,-7);
const double Z0 = 376.73;
const std::complex<double> jj(0.0, 1.0);

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    double CalcAEM(double k);
    double CalcREM(double f);
    double CalcBEM(double f, double k);
    double CalcRE (double f, double r);
    double CalcBE (double f, double k, double r);
    double CalcRM (double f, double r);
    double CalcBM (double f, double k, double r);

    std::complex<double> CalcZ(double mu, double eps);
    std::complex<double> CalcZ(double mu, double eps, double f, double sigma);
    std::complex<double> CalcGamma(double mu, double eps, double f);
    std::complex<double> CalcGamma(double mu, double eps, double f, double sigma);
    std::complex<double>* CalcA(double mu, double eps, double l, double f);
    std::complex<double>* CalcA(double mu, double eps, double l, double f, double sigma);
    std::complex<double>* CalcGeneralA(double f);

    double CalcSE(double f);

private slots:
	void on_actionHome_triggered();
	void on_actionOpen_DB_triggered();
	
	void point_selected(QtDataVisualization::QSurface3DSeries* series);
	void on_nearButton_clicked();
	void on_farButton_clicked();

	void on_magneticButton_clicked();
	void on_electricButton_clicked();
	void on_nearBackButton_clicked();

	void on_metalButton_clicked();
	void on_compositeButton_clicked();
	void on_farBackButton_clicked();

	void on_run2DCalcButton_nz_clicked();
	void on_run3DCalcButton_nz_clicked();
	void on_nearSettingsBackButton_clicked();

	void on_addPlateButton_clicked();
	void on_deletePlateButton_clicked();
	void on_run2DCalcButton_fz_clicked();
	void on_farSettingsBackButton_clicked();

	void on_calc2DBackButton_clicked();
	void on_calc3DBackButton_clicked();
	
private:
    int _curMode = 0;
    int _platesNum = 1;
    QtDataVisualization::Q3DSurface* _surface;
    QWidget* _container;
    Ui::MainWindow* _ui;
};

#endif // MAINWINDOW_H
