#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "MatLibForm.h"

using namespace QtDataVisualization;

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), _ui(new Ui::MainWindow)
{
    _ui->setupUi(this);
    QTabBar* tabBar = _ui->tabWidget->findChild<QTabBar*>();
    tabBar->hide();
    _ui->secondPlate->hide();
    _ui->thirdPlate->hide();
    _ui->fourthPlate->hide();
    setWindowTitle("MultilayerPlateSE");
}

MainWindow::~MainWindow()
{
    delete _ui;
}

double MainWindow::CalcAEM(double k)
{
    return (8.69 * k * _ui->t_nz->text().toDouble()) / sqrt(2);
}

double MainWindow::CalcREM(double f)
{
    double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
    double Zm = sqrt((2 * M_PI * f * mu) / _ui->sigma_nz->text().toDouble());
    double Zd = sqrt(mu0 / eps0);
    return 20 * log10(abs(pow(Zd + Zm, 2) / (4 * Zd * Zm)));
}

double MainWindow::CalcBEM(double f, double k)
{
    double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
    double sigma = _ui->sigma_nz->text().toDouble();
    double t = _ui->t_nz->text().toDouble();
    double Zm = sqrt((2 * M_PI * f * mu) / sigma);
    double Zd = sqrt(mu0 / eps0);
    return 20 * log10(abs(1 - (pow((Zd - Zm) / (Zd + Zm), 2) * pow(M_E, -2 * k * t))));
}

double MainWindow::CalcRE(double f, double r)
{
	double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
	double sigma = _ui->sigma_nz->text().toDouble();
	double Zm = sqrt((2 * M_PI * f * mu) / sigma);
    double Zd = 1 / (2 * M_PI * f * eps0 * r);
    return 20 * log10(abs(pow(Zd + Zm, 2) / (4 * Zd * Zm)));
}

double MainWindow::CalcBE(double f, double k, double r)
{
    double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
    double sigma = _ui->sigma_nz->text().toDouble();
    double t = _ui->t_nz->text().toDouble();
	double Zm = sqrt((2 * M_PI * f * mu) / sigma);
    double Zd = 1 / (2 * M_PI * f * eps0 * r);
    return 20 * log10(abs(1 - (pow((Zd - Zm) / (Zd + Zm), 2) * pow(M_E, -2 * k * t))));
}

double MainWindow::CalcRM(double f, double r)
{
    double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
    double sigma = _ui->sigma_nz->text().toDouble();
    double Zm = sqrt((2 * M_PI * f * mu) / sigma);
    double Zd = 2 * M_PI * f * mu0 * r;
    return 20 * log10(abs(pow(Zd + Zm, 2) / (4 * Zd * Zm)));
}

double MainWindow::CalcBM(double f, double k, double r)
{
	double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
    double sigma = _ui->sigma_nz->text().toDouble();
    double t = _ui->t_nz->text().toDouble();
    double Zm = sqrt((2 * M_PI * f * mu) / sigma);
    double Zd = 2 * M_PI * f * mu0 * r;
    return 20 * log10(abs(1 - (pow((Zd - Zm) / (Zd + Zm), 2) * pow(M_E, -2 * k * t))));
}

// Composite material.
std::complex<double> MainWindow::CalcZ(double mu, double eps)
{
    return Z0 * sqrt(mu / eps);
}

// Metal.
std::complex<double> MainWindow::CalcZ(double mu, double eps, double f, double sigma)
{
    return sqrt((jj * 2.0 * M_PI * f * mu * mu0) / (sigma + jj * 2.0 * M_PI * f * eps * eps0));
}

// Composite material.
std::complex<double> MainWindow::CalcGamma(double mu, double eps, double f)
{
    return jj * 2.0 * M_PI * f * sqrt(mu * eps) * pow(c, -1);
}

// Metal.
std::complex<double> MainWindow::CalcGamma(double mu, double eps, double f, double sigma)
{
    return sqrt(jj * 2.0 * M_PI * f * mu * mu0 * (sigma + jj * 2.0 * M_PI * f * eps * eps0));
}

// Composite material.
std::complex<double>* MainWindow::CalcA(double mu, double eps, double l, double f)
{
    std::complex<double> Z = CalcZ(mu, eps);
    std::complex<double> k = CalcGamma(mu, eps, f);
    std::complex<double>* A = new std::complex<double> [4];

    A[0] = cosh(k * l);
    A[1] = Z * sinh(k * l);
    A[2] = sinh(k * l) / Z;
    A[3] = cosh(k * l);

    return A;
}

// Metal.
std::complex<double>* MainWindow::CalcA(double mu, double eps, double l, double f, double sigma)
{
    std::complex<double> Z = CalcZ(mu, eps, f, sigma);
    std::complex<double> k = CalcGamma(mu, eps, f, sigma);
    std::complex<double>* A = new std::complex<double> [4];

    A[0] = cosh(k * l);
    A[1] = Z * sinh(k * l);
    A[2] = sinh(k * l) / Z;
    A[3] = cosh(k * l);

    return A;
}

std::complex<double>* MainWindow::CalcGeneralA(double f)
{
    std::complex<double>* mainMatrix = new std::complex<double> [4];
    std::complex<double>* tempA = new std::complex<double> [4];
    std::complex<double>* tempB = new std::complex<double> [4];

    double mu = _ui->mu_fz_1->text().toDouble();
    double eps = _ui->eps_fz_1->text().toDouble();
    double l = _ui->l_fz_1->text().toDouble();

    if (_platesNum >= 1)
    {
        if (_curMode == 3)
        {
            mainMatrix = CalcA(mu, eps, l, f, _ui->sigma_fz_1->text().toDouble());
        }
        else if (_curMode == 4)
        {
            mainMatrix = CalcA(mu, eps, l, f);
        }
    }

    if (_platesNum >= 2)
    {
        tempA = mainMatrix;
        mu = _ui->mu_fz_2->text().toDouble();
        eps = _ui->eps_fz_2->text().toDouble();
        l = _ui->l_fz_2->text().toDouble();

        if (_curMode == 3)
        {
            tempB = CalcA(mu, eps, l, f, _ui->sigma_fz_2->text().toDouble());
        }
        else if (_curMode == 4)
        {
            tempB = CalcA(mu, eps, l, f);            
        }

        mainMatrix[0] = tempA[0] * tempB[0] + tempA[1] * tempB[2];
        mainMatrix[1] = tempA[0] * tempB[1] + tempA[1] * tempB[3];
        mainMatrix[2] = tempA[2] * tempB[0] + tempA[3] * tempB[2];
        mainMatrix[3] = tempA[2] * tempB[1] + tempA[3] * tempB[3];
    }

    if (_platesNum >= 3)
    {
        tempA = mainMatrix;
        mu = _ui->mu_fz_3->text().toDouble();
        eps = _ui->eps_fz_3->text().toDouble();
        l = _ui->l_fz_3->text().toDouble();

        if (_curMode == 3)
        {
            tempB = CalcA(mu, eps, l, f, _ui->sigma_fz_3->text().toDouble());
        }
        else if (_curMode == 4)
        {
            tempB = CalcA(mu, eps, l, f);
        }

        mainMatrix[0] = tempA[0] * tempB[0] + tempA[1] * tempB[2];
        mainMatrix[1] = tempA[0] * tempB[1] + tempA[1] * tempB[3];
        mainMatrix[2] = tempA[2] * tempB[0] + tempA[3] * tempB[2];
        mainMatrix[3] = tempA[2] * tempB[1] + tempA[3] * tempB[3];
    }

    if (_platesNum == 4)
    {
        tempA = mainMatrix;
        mu = _ui->mu_fz_4->text().toDouble();
        eps = _ui->eps_fz_4->text().toDouble();
        l = _ui->l_fz_4->text().toDouble();

        if (_curMode == 3)
        {
            tempB = CalcA(mu, eps, l, f, _ui->sigma_fz_4->text().toDouble());
        }
        else if (_curMode == 4)
        {
            tempB = CalcA(mu, eps, l, f);
        }

        mainMatrix[0] = tempA[0] * tempB[0] + tempA[1] * tempB[2];
        mainMatrix[1] = tempA[0] * tempB[1] + tempA[1] * tempB[3];
        mainMatrix[2] = tempA[2] * tempB[0] + tempA[3] * tempB[2];
        mainMatrix[3] = tempA[2] * tempB[1] + tempA[3] * tempB[3];
    }

    return mainMatrix;
}

double MainWindow::CalcSE(double f)
{
    std::complex<double>* A = new std::complex<double> [4];
    A = CalcGeneralA(f);
    return 20 * log10(abs((A[0] * Z0 + A[1] + A[2] * Z0 * Z0 + A[3] * Z0) / (2.0 * Z0)));
}

void MainWindow::on_actionHome_triggered()
{
	_ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_actionOpen_DB_triggered()
{
	MatLibForm* matLib = new MatLibForm;
	matLib->setWindowTitle("Material Library");
	matLib->show();
}

void MainWindow::point_selected(QSurface3DSeries* series)
{
    qDebug() << "there";
}

void MainWindow::on_nearButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(1);
}

void MainWindow::on_farButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(3);
}

void MainWindow::on_magneticButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(2);
    _curMode = 1;
}

void MainWindow::on_electricButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(2);
    _curMode = 2;
}

void MainWindow::on_nearBackButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_metalButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(4);
	_ui->sigma_fz_1->show();
	_ui->label_14->show();
	_ui->sigma_fz_2->show();
	_ui->label_19->show();
	_ui->sigma_fz_3->show();
	_ui->label_24->show();
	_ui->sigma_fz_4->show();
	_ui->label_29->show();
    _curMode = 3;
}

void MainWindow::on_compositeButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(4);
	_ui->sigma_fz_1->hide();
	_ui->label_14->hide();
	_ui->sigma_fz_2->hide();
	_ui->label_19->hide();
	_ui->sigma_fz_3->hide();
	_ui->label_24->hide();
	_ui->sigma_fz_4->hide();
	_ui->label_29->hide();
    _curMode = 4;
}

void MainWindow::on_farBackButton_clicked()
{
	_ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_run2DCalcButton_nz_clicked()
{
	double mu = _ui->mu_nz->text().toDouble() * mu0;
	double sigma = _ui->sigma_nz->text().toDouble();
	double fmin = _ui->fMin_nz->text().toDouble();
	double fmax = _ui->fMax_nz->text().toDouble();
	double fstep = _ui->fStep_nz->text().toDouble();
	double r = _ui->r_nz->text().toDouble();
	QVector<double> A, B, R, SE, fs;

	// Magnetic calculation.
	if (_curMode == 1)
	{
		for (double f = fmin; f <= fmax; f += fstep)
		{
			auto k = abs(sqrt(2 * M_PI * f * mu * sigma));
			A.push_back(CalcAEM(k));
			R.push_back(CalcRM(f, r));
			B.push_back(CalcBM(f, k, r));
			SE.push_back(A.last() + R.last() + B.last());
			fs.push_back(f);
		}
	}

	// Electric calculation.
	else if (_curMode == 2)
	{
		for (double f = fmin; f <= fmax; f += fstep)
		{
			auto k = abs(sqrt(2 * M_PI * f * mu * sigma));
			A.push_back(CalcAEM(k));
			R.push_back(CalcRE(f, r));
			B.push_back(CalcBE(f, k, r));
			SE.push_back(A.last() + R.last() + B.last());
			fs.push_back(f);
		}
	}

	_ui->plot->addGraph();
	_ui->plot->graph(0)->setData(fs, A);
	_ui->plot->graph(0)->setPen(QPen(Qt::blue));
	_ui->plot->graph(0)->rescaleAxes(true);
	_ui->plot->addGraph();
	_ui->plot->graph(1)->setData(fs, B);
	_ui->plot->graph(1)->setPen(QPen(Qt::green));
	_ui->plot->graph(1)->rescaleAxes(true);
	_ui->plot->addGraph();
	_ui->plot->graph(2)->setData(fs, R);
	_ui->plot->graph(2)->setPen(QPen(Qt::black));
	_ui->plot->graph(2)->rescaleAxes(true);
	_ui->plot->addGraph();
	_ui->plot->graph(3)->setData(fs, SE);
	_ui->plot->graph(3)->setPen(QPen(Qt::red));
	_ui->plot->graph(3)->rescaleAxes(true);
	_ui->plot->replot();
	_ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	_ui->tabWidget->setCurrentIndex(5);
}

void MainWindow::on_run3DCalcButton_nz_clicked()
{
	_surface = new Q3DSurface();
	_container = QWidget::createWindowContainer(_surface);

	_surface->setHorizontalAspectRatio(1.0);
	_surface->axisX()->setTitle("Distance to shield, m");
	_surface->axisY()->setTitle("SE, dB");
	_surface->axisZ()->setTitle("Frequency, Hz");
	_surface->axisZ()->setLabelFormat("%.0f");
	_surface->axisX()->setLabelAutoRotation(30.0f);
	_surface->axisY()->setLabelAutoRotation(30.0f);
	_surface->axisZ()->setLabelAutoRotation(30.0f);
	_surface->axisX()->setTitleVisible(true);
	_surface->axisY()->setTitleVisible(true);
	_surface->axisZ()->setTitleVisible(true);
	_ui->surfLayout->addWidget(_container, 1);
	_ui->tabWidget->setCurrentIndex(6);

	double r = _ui->r_nz->text().toDouble();
	double mu = _ui->mu_nz->text().toDouble() * 4 * M_PI * pow(10, -7);
	double sigma = _ui->sigma_nz->text().toDouble();
	double fmin = _ui->fMin_nz->text().toDouble();
	double fmax = _ui->fMax_nz->text().toDouble();
	int numPoints = _ui->fStep_nz->text().toInt();
	int step = (fmax - fmin) / numPoints;

	QSurfaceDataArray* data = new QSurfaceDataArray;
	QSurface3DSeries* series = new QSurface3DSeries;
	double x, y, z;

	// Magnetic calculation.
	if (_curMode == 1)
	{
		for (int i = 0; i < numPoints; i++)
		{
			QSurfaceDataRow *dataRow = new QSurfaceDataRow(numPoints);

			x = r;
			for (int j = 0; j < numPoints; j++)
			{
				auto k = abs(sqrt(2 * M_PI * step * (j + 1) * mu * sigma));
				y = CalcAEM(k) + CalcRM(step * (j + 1), r) + CalcBM(step * (j + 1), k, r);
				z = (j + 1) * step;
				(*dataRow)[j].setPosition(QVector3D(x, y, z));
			}
			*data << dataRow;
			r += 1e-3;
		}
	}

	// Electric calculation.
	else if (_curMode == 2)
	{
		for (int i = 0; i < numPoints; i++)
		{
			QSurfaceDataRow* dataRow = new QSurfaceDataRow(numPoints);

			x = r;
			for (int j = 0; j < numPoints; j++)
			{
				auto k = abs(sqrt(2 * M_PI * step * (j + 1) * mu * sigma));
				y = CalcAEM(k) + CalcRE(step * (j + 1), r) + CalcBE(step*(j + 1), k, r);
				z = (j + 1) * step;
				(*dataRow)[j].setPosition(QVector3D(x, y, z));
			}
			*data << dataRow;
			r += 1e-3;
		}
	}

	_surface->axisX()->setRange(0, r);
	_surface->axisZ()->setRange(0, step * numPoints);

	series->dataProxy()->resetArray(data);
	series->setDrawMode(QSurface3DSeries::DrawSurface);
	_surface->addSeries(series);

	Q3DTheme* theme = new Q3DTheme(Q3DTheme::ThemeStoneMoss);
	_surface->setActiveTheme(theme);
	QLinearGradient gr;
	gr.setColorAt(0.0, QColor(11, 30, 32));
	gr.setColorAt(0.4, QColor(50, 105, 115));
	gr.setColorAt(0.8, QColor(69, 213, 233));
	_surface->seriesList().at(0)->setBaseGradient(gr);
	_surface->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

	_ui->surfLayout->update();
}

void MainWindow::on_nearSettingsBackButton_clicked()
{
	_ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_addPlateButton_clicked()
{
    if (_ui->secondPlate->isHidden())
    {
        _ui->secondPlate->show();
        _platesNum++;
    }
    else if (_ui->thirdPlate->isHidden())
    {
        _ui->thirdPlate->show();
        _platesNum++;
    }
    else
    {
        _ui->fourthPlate->show();
        _platesNum++;
    }
}

void MainWindow::on_deletePlateButton_clicked()
{
    if (!_ui->fourthPlate->isHidden())
    {
        _ui->fourthPlate->hide();
        _platesNum--;
    }
    else if (!_ui->thirdPlate->isHidden())
    {
        _ui->thirdPlate->hide();
        _platesNum--;
    }
    else
    {
        _ui->secondPlate->hide();
        _platesNum--;
    }
}

void MainWindow::on_run2DCalcButton_fz_clicked()
{
    _ui->plot->clearGraphs();
    double fmin = _ui->fMin_fz->text().toDouble();
    double fmax = _ui->fMax_fz->text().toDouble();
    double fstep = _ui->fStep_fz->text().toDouble();
    QVector<double> SE, fs;

    for (double f = fmin; f <= fmax; f += fstep)
    {
        fs.push_back(f);
        SE.push_back(CalcSE(f));
    }

    _ui->plot->addGraph();
    _ui->plot->graph(0)->setData(fs, SE);
    _ui->plot->graph(0)->setPen(QPen(Qt::red));
    _ui->plot->graph(0)->rescaleAxes(true);
    _ui->plot->replot();
    _ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    _ui->tabWidget->setCurrentIndex(5);
}

void MainWindow::on_farSettingsBackButton_clicked()
{
	_ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_calc2DBackButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(0);
}

void MainWindow::on_calc3DBackButton_clicked()
{
    _ui->tabWidget->setCurrentIndex(0);
    _ui->surfLayout->removeWidget(_container);
}
