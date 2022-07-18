#ifndef MATLIBFORM_H
#define MATLIBFORM_H

#include <QWidget>

namespace Ui { class MatLibForm; }

class MatLibForm : public QWidget
{
    Q_OBJECT

public:
    explicit MatLibForm(QWidget* parent = nullptr);
    ~MatLibForm();

private slots:
    void on_tableWidget_cellClicked(int row, int column);

private:
    Ui::MatLibForm* _ui;
};

#endif // MATLIBFORM_H
