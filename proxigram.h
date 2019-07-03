#ifndef PROXIGRAM_H
#define PROXIGRAM_H

#include <QWidget>
#include "IPlugin.h"
#include "IStateCreator.h"
#include <QObject>
#include "CurvePlotter.h"
#include "geometry.h"
#include "proxigramstools.h"
#include <QLabel>
#include <QDialog>

///TODO На будущее - реализовать буфер результатов проксиграммы

namespace Ui {
class Proxigram;
}

class ErrorDialog : public QDialog
{
    Q_OBJECT

public:

    ErrorDialog(QWidget* pwgt = 0);

    void settest(const QString& message);

private:

    QLabel* ptr_error;
   // QLabel* ptr_pix;
};

class Proxigram : public QWidget
{
    Q_OBJECT

public:
    explicit Proxigram(QWidget *parent = 0);
    ~Proxigram();
    void Plotter(const ElementSymbolInfo elt, QVector<double>* X, QVector<double>* Y, double Max);
  //  void Plot_kraken(std::vector<std::vector<double>> data);
private slots:
    void on_Kraken_button_clicked();

    void on_NumBins_Edit_textEdited(const QString &arg1);

    void on_Xmax_Edit_textChanged(const QString &arg1);

    void on_Xmin_Edit_textChanged(const QString &arg1);

    void on_elt_widget_SelectionChanged(const QVector<bool> &Elements);

    void on_Separation_button_clicked();

    void on_pushButton_clicked();

    void on_Show_button_clicked();

    void on_Hide_button_clicked();

    void on_lineEdit_returnPressed();

    void on_horizontalSlider_valueChanged(int value);

    void on_horizontalSlider_sliderReleased();

private:
    std::vector<size_t> showID;
    bool isAtomInside(const Point&, const std::vector<std::size_t>&);
    double calculateLength(const std::vector<std::size_t>&);
    void setupIsosurface();
    void computeIsoAtoms();
    void calculatedata();
    void plotdata();
    double getMinimumDistance(const Point& atopPoint);
    std::vector<std::vector<double>> _plotdata;
    QVector<double> XData;
    double xmax, xmin;
    double step;
    QVector<bool> ChosenAtomsA;
    std::vector<bool> ForClusters;
    Ui::Proxigram *ui;
    void AElemntSelectorSlot(QVector<bool> selectedElements);
    proxiTriangles proxi;
    GraphSeparator gs;
    std::vector<Isosurface> isos;
};


class ProxigramPlugin: public IPlugin {
 Q_OBJECT
public:
    ProxigramPlugin();
    void connectToGui(IConnectableWindow* windows);
private:
     Proxigram* wid;
};

#endif // PROXIGRAM_H
