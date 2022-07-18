#include "MatLibForm.h"
#include "ui_MatLibForm.h"

MatLibForm::MatLibForm(QWidget* parent) : QWidget(parent), _ui(new Ui::MatLibForm)
{
    _ui->setupUi(this);
}

MatLibForm::~MatLibForm()
{
    delete _ui;
}

void MatLibForm::on_tableWidget_cellClicked(int row, int column)
{
    _ui->tableWidget->selectRow(row);
    QString name  = _ui->tableWidget->item(row, 0)->text();
    QString descr = _ui->tableWidget->item(row, 2)->text();
    _ui->mat_name->setText(name);
    _ui->descript_browser->setText(descr);
    if(_ui->tableWidget->item(row, 1)->text() == "Metal")
		_ui->comboBox->setCurrentIndex(0);
    else
		_ui->comboBox->setCurrentIndex(1);
}
