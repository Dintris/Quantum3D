#include "proxigram.h"
#include "ui_proxigram.h"
#include "naivepluginManager.h"
#include "VoxelDraw.h"
#include "Overlord3D.h"
#include <iostream>
#include <QMessageBox>
#include <ClusteIdentificatorPlugin.h>

ErrorDialog::ErrorDialog(QWidget *pwgt):
    QDialog(pwgt, Qt::WindowTitleHint | Qt::WindowSystemMenuHint)
{
    ptr_error = new QLabel;

    QPushButton* pcmdOk     = new QPushButton("&Больше так не буду");
    QPushButton* pcmdCancel = new QPushButton("&Я все сделал правильно");

    connect(pcmdOk, SIGNAL(clicked()), SLOT(accept()));
    connect(pcmdCancel, SIGNAL(clicked()), SLOT(reject()));

    QGridLayout* layout = new QGridLayout;

    layout->addWidget(ptr_error, 1, 0);
    layout->addWidget(pcmdOk, 2, 0);
    layout->addWidget(pcmdCancel, 3, 0);
    setLayout(layout);

}

void ErrorDialog::settest(const QString &message)
{
    ptr_error->setText(message);
}

ProxigramPlugin::ProxigramPlugin()
: wid(new Proxigram)
{

}

void ProxigramPlugin::connectToGui(IConnectableWindow* windows)
{
    auto csw = windows->CreateDockableWidget(wid, "Proxigram");
    QIcon icon;
    QAction* toolbarAction = windows->CreateToolbarButton(icon, tr("Proxigram"), tr("Proxigram"));
    QAction* menubarAction = windows->CreateMenubarSubmenu(icon, "Proxigram");
    QObject::connect(toolbarAction, &QAction::triggered, [csw](bool) {csw ->toggle(); });
    QObject::connect(menubarAction, &QAction::triggered, [csw](bool) {csw ->toggle(); });
}

Proxigram::Proxigram(QWidget *parent) :
    QWidget(parent),
    gs(0),
    ui(new Ui::Proxigram)
{
    ui->setupUi(this);
    xmax = 10.0;
    xmin = -10.0;
    step = 1;
    ui->Xmax_Edit->setText(QString::number(xmax));
    ui->Xmin_Edit->setText(QString::number(xmin));
    ui->NumBins_Edit->setText(QString::number(step));
    ui->graph_widget->plot()->setAxisTitle(2,"Distance");
    ui->graph_widget->plot()->setAxisTitle(0,"Percentage");

}

Proxigram::~Proxigram()
{
    delete ui;
}

void Proxigram::Plotter(const ElementSymbolInfo elt, QVector<double> *X, QVector<double> *Y, double Max)
{
    double Diapason = xmax - xmin;
    ui->graph_widget->plot()->setAxisScale(QwtPlot::yLeft,0.f,Max,Max/5.0);
    ui->graph_widget->plot()->setAxisScale(QwtPlot::xBottom,-Diapason/2.f,Diapason/2.f,Diapason/10.f);
    ui->graph_widget->GetCurve(elt.symbol)->setPen(QColor(elt.fRed * 255, elt.fGreen * 255, elt.fBlue * 255));
    ui->graph_widget->UpdateCurveData(elt.symbol, *X, *Y);
    ui->graph_widget->replot();
}

void Proxigram::on_Kraken_button_clicked()
{
    calculatedata();
}

void Proxigram::on_NumBins_Edit_textEdited(const QString &arg1)
{
    step = arg1.toDouble();
}

void Proxigram::on_Xmax_Edit_textChanged(const QString &arg1)
{
    xmax = arg1.toDouble();
}

void Proxigram::on_Xmin_Edit_textChanged(const QString &arg1)
{
    xmin = arg1.toDouble();
}

void Proxigram::on_elt_widget_SelectionChanged(const QVector<bool> &Elements)
{
    AElemntSelectorSlot(Elements);
}

void Proxigram::calculatedata()
{
    _plotdata.clear();
    XData.clear();

    if(xmin >= xmax)
    {
        ErrorDialog* ptr_error = new ErrorDialog;
        ptr_error->settest("Некорректные значения границ диапазона");
        if(ptr_error->exec() == QDialog::Accepted)
        {
            QMessageBox::information(0, "Информация", "Вот и молодец");
        }
        else
        {
            QMessageBox::information(0, "Информация", "А вот и нет");
            delete ptr_error;
            return;
        }
        delete ptr_error;
        return;
    }

    StateMachineReader smr("Voxel Draw");
    const State& state = smr.getState();
    if(!state.EnvironmentLoaded())
        return;

    const std::vector<ElementSymbolInfo>& elts = state.GetEnvironment().GetElements();
    const EIntervalsVector& eltsIntervals =  state.GetEnvironment().GetElementIntervals();
    const auto& atoms = state.GetAtoms();
    Overlord3D* p1 = dynamic_cast<Overlord3D*>(NaivePluginLoader::instance().get("Overlord3D"));
    VoxelDraw* p = dynamic_cast<VoxelDraw*>(p1->GetDrawWidget()->GetVoxelDraw());

    proxi.buildTriangle(p->_elts[0].draw->pointSet, p->_elts[0].draw->normals);

    if(p->_elts[0].draw->pointSet.size() == 0) return;

    const auto& cL = gs.getCompsList();
    if(cL.size() == 0) return;

    //Выбираем треугольник только в выбранных изоповерхностях
    std::size_t selectSize = showID.size();
    std::vector<std::pair<Triangle, QVector3D>> selectedTriangles;
    for(std::size_t i = 0; i < selectSize; ++i)
    {
        std::size_t currentIndex = showID[i];
        const auto& currentIsosurface = cL[currentIndex];
        for(std::size_t i = 0; i < currentIsosurface.size(); ++i)
        {
            //Заносим треугольники выбранной изоповерхности
            std::size_t currentTriangleIndex = currentIsosurface[i];
            auto insertedItem = std::make_pair(proxi[currentTriangleIndex].first, proxi[currentTriangleIndex].second);
            selectedTriangles.emplace_back(insertedItem);
        }
    }

    std::cout << "Number of Selected IsoSurface = " << selectSize << std::endl;
    if(selectedTriangles.size() == 0) return;

    std::size_t NumBins = (xmax - xmin) / step;

    _plotdata.resize(elts.size(), std::vector<double>(NumBins, 0.0));

    std::vector<proxiDatagram> datagramStore;
    datagramStore.resize(elts.size());
    for(size_t i = 0; i < datagramStore.size(); ++i)
    {
        datagramStore[i].update(NumBins, xmin, xmax);
    }

    for(const auto& atom : atoms)
    {
        //=================================================================
        Point atomPoint(atom.fPx, atom.fPy, atom.fPz);
        const std::vector<size_t>& intervalsIds = state.GetEnvironment().GetIntervalsIdsForMass(atom.GetMass());
        double distance = selectedTriangles[0].first.ComputeDistance(atomPoint);
        QVector3D normal = selectedTriangles[0].second;
        Triangle tri = selectedTriangles[0].first;
        std::size_t SIZE = selectedTriangles.size();
        for(size_t i = 1; i < SIZE; ++i)
        {
            double d = selectedTriangles[i].first.ComputeDistance(atomPoint);
            if(d < distance)
            {
                distance = d;
                normal = selectedTriangles[i].second;
                tri = selectedTriangles[i].first;
            }
        }
        QVector3D ptmp(atom.fPx - tri.vertex[0]->X(), atom.fPy - tri.vertex[0]->Y() , atom.fPz - tri.vertex[0]->Z());
        if(QVector3D::dotProduct(ptmp, normal) < 0)
        {
            distance = -distance;
        }
        //====================================================================
        for(size_t i = 0; i < intervalsIds.size(); ++i)
        {
            const size_t iNumElt = eltsIntervals[intervalsIds[i]].iNumElt;
            const float increment = eltsIntervals[intervalsIds[i]].fIncrementForConcentration;
            datagramStore[iNumElt].fillValue(distance, increment);
        }
    }

    for(size_t i = 0; i < datagramStore.size(); ++i)
    {
        _plotdata[i] = datagramStore[i].getPlotData();
    }

    std::vector<double> sumConcVector = getSumConcentration(_plotdata);

    for(size_t i = 0; i < _plotdata[0].size(); ++i)
    {
        for(size_t k = 0; k < _plotdata.size(); ++k)
        {
            if(sumConcVector[i] != 0.0)
                _plotdata[k][i] /= sumConcVector[i];
        }
    }

    XData.clear();
    XData.resize(NumBins);
    XData[0] = xmin + step/2;
    for(size_t i = 1; i < NumBins; ++i){
        XData[i] = XData[i-1] + step;
    }

    plotdata();
}

void Proxigram::plotdata()
{
    StateMachineReader smr("Voxel Draw");
    const State& state = smr.getState();
    if(!state.EnvironmentLoaded())
        return;
    const std::vector<ElementSymbolInfo> elts = state.GetEnvironment().GetElements();
    auto AElements = ChosenAtomsA ;
    QVector<double> Y;

    for(int i =0; i< elts.size(); ++i){
        Y.clear();
       if(AElements[i]){
        for(int j =0; j<_plotdata[i].size(); ++j){
            Y.push_back(_plotdata[i][j]*100.0);
        }

       Plotter(elts[i],&XData, &Y, 100);
    }
    }
}

void Proxigram::AElemntSelectorSlot(QVector<bool> selectedElements) {
    ChosenAtomsA = selectedElements;
    if(!_plotdata.empty()){
    ui->graph_widget->RemoveCurves();
    plotdata();
    }
    for(size_t i=0;i<selectedElements.size();++i)
    {
        if(selectedElements[i])
        {
            ForClusters.push_back(true);
        }
        else
        {
            ForClusters.push_back(false);
        }
    }
}

void Proxigram::on_Separation_button_clicked()
{
    StateMachineReader smr("Voxel Draw");
    const State& state = smr.getState();
    if(!state.EnvironmentLoaded())
        return;

    const std::vector<ElementSymbolInfo> elts = state.GetEnvironment().GetElements();
    Overlord3D* p1 = dynamic_cast<Overlord3D*>(NaivePluginLoader::instance().get("Overlord3D"));
    VoxelDraw* p = dynamic_cast<VoxelDraw*>(p1->GetDrawWidget()->GetVoxelDraw());

    if(p->_elts[0].draw->pointSet.size() == 0) return;

    proxi.buildTriangle(p->_elts[0].draw->pointSet, p->_elts[0].draw->normals);

    gs.update(proxi.getSize());
    isos.clear();

    gs.PrintState();

    for(std::size_t i = 0; i < proxi.getSize(); ++i)
    {
        const auto& tri1 = proxi[i].first;
        for(std::size_t k = 0; k < proxi.getSize(); ++k)
        {
            const auto& tri2 = proxi[k].first;
            if(isTriangleConnect(tri1, tri2))
            {
                gs.inputEdge(i, k);
            }
        }
    }

    gs.findComps(1);

    const auto& cL = gs.getCompsList();

    setupIsosurface();
    computeIsoAtoms();

    const auto& atoms = state.GetAtoms();

    std::vector<double> length;
    length.reserve(cL.size());

    for(std::size_t i = 0; i < cL.size(); ++i)
    {
        double value = calculateLength(cL[i]);
        length.push_back(value);
        isos[i].setLength(value);
    }

    auto minMax = std::minmax_element(std::begin(length), std::end(length));

    int minValue = *minMax.first * 10;
    int maxValue = *minMax.second * 10;

    ui->label_3->setText("Min Length = " + QString::number(*minMax.first));
    ui->label_4->setText("Max Length = " + QString::number(*minMax.second));
    ui->horizontalSlider->setRange(minValue, maxValue);
    ui->IsosurfaceIDList->disconnect();
    ui->IsosurfaceIDList->clearContents();
    ui->IsosurfaceIDList->setRowCount(cL.size());
    for(size_t i = 0; i < cL.size(); ++i)
    {
        ui->IsosurfaceIDList->setItem(i, 0, new QTableWidgetItem(QString::number(i)));
        QTableWidgetItem* pit = new QTableWidgetItem();
        pit->setFlags( pit->flags() | Qt::ItemIsUserCheckable);
        pit->setCheckState(Qt::Unchecked);
        ui->IsosurfaceIDList->setItem(i, 1, pit);
        ui->IsosurfaceIDList->setItem(i, 2, new QTableWidgetItem(QString::number(isos[i].getAtomsNumber())));
        ui->IsosurfaceIDList->setItem(i, 3, new QTableWidgetItem(QString::number(isos[i].getLength())));
    }
    QObject::connect(ui->IsosurfaceIDList, &QTableWidget::itemChanged, [this](){on_pushButton_clicked();});
    on_pushButton_clicked();
}

void Proxigram::on_pushButton_clicked()
{
    showID.clear();

    for(size_t i=0; i<ui->IsosurfaceIDList->rowCount();++i)
    {
        if(ui->IsosurfaceIDList->item(i,1)->checkState() == Qt::Checked){
            showID.push_back(ui->IsosurfaceIDList->item(i,0)->text().toUInt());
        }
    }
}

void Proxigram::on_Show_button_clicked()
{
    for(size_t i=0; i<ui->IsosurfaceIDList->rowCount();++i){
            ui->IsosurfaceIDList->item(i,1)->setCheckState(Qt::Checked);
    }
}

void Proxigram::on_Hide_button_clicked()
{
    for(size_t i=0; i<ui->IsosurfaceIDList->rowCount();++i){
            ui->IsosurfaceIDList->item(i,1)->setCheckState(Qt::Unchecked);
    }
}

bool Proxigram::isAtomInside(const Point& atom, const std::vector<size_t>& isosurface)
{
    std::size_t SIZE = isosurface.size();
    for(std::size_t i = 0; i < SIZE; ++i)
    {
        QVector3D& normal = proxi[isosurface[i]].second;
        Triangle& tri = proxi[isosurface[i]].first;
        QVector3D ptmp(atom.X() - tri.vertex[0]->X(),
                       atom.Y() - tri.vertex[0]->Y() ,
                       atom.Z() - tri.vertex[0]->Z());
        if(QVector3D::dotProduct(normal, ptmp) > 0)
            return false;
    }
    return true;
}

double Proxigram::calculateLength(const std::vector<std::size_t>& triangles)
{
    std::size_t SIZE = triangles.size();
    if(SIZE == 0) return 0.0;

    std::size_t index = triangles[0];
    auto& triangle = proxi[index].first;

    double xmax, xmin, ymax, ymin, zmax, zmin;
    double xLength, yLength, zLength;

    xmax = triangle.vertex[0]->X();
    xmin = triangle.vertex[0]->X();
    ymax = triangle.vertex[0]->Y();
    ymin = triangle.vertex[0]->Y();
    zmax = triangle.vertex[0]->Z();
    zmin = triangle.vertex[0]->Z();

    for(std::size_t i = 1; i < 3; ++i)
    {
        if(xmax < triangle.vertex[i]->X())
            xmax = triangle.vertex[i]->X();
        if(xmin > triangle.vertex[i]->X())
            xmin = triangle.vertex[i]->X();

        if(ymax < triangle.vertex[i]->Y())
            ymax = triangle.vertex[i]->Y();
        if(ymin > triangle.vertex[i]->Y())
            ymin = triangle.vertex[i]->Y();

        if(zmax < triangle.vertex[i]->Z())
            zmax = triangle.vertex[i]->Z();
        if(zmin > triangle.vertex[i]->Z())
            zmin = triangle.vertex[i]->Z();
    }

    for(std::size_t i = 1; i < SIZE; ++i)
    {
        std::size_t index = triangles[i];
        triangle = proxi[index].first;
        for(std::size_t k = 0; k < 3; ++k)
        {
            if(xmax < triangle.vertex[k]->X())
                xmax = triangle.vertex[k]->X();
            if(xmin > triangle.vertex[k]->X())
                xmin = triangle.vertex[k]->X();

            if(ymax < triangle.vertex[k]->Y())
                ymax = triangle.vertex[k]->Y();
            if(ymin > triangle.vertex[k]->Y())
                ymin = triangle.vertex[k]->Y();

            if(zmax < triangle.vertex[k]->Z())
                zmax = triangle.vertex[k]->Z();
            if(zmin > triangle.vertex[k]->Z())
                zmin = triangle.vertex[k]->Z();
        }
    }

    xLength = xmax - xmin;
    yLength = ymax - ymin;
    zLength = zmax - zmin;

    return std::max({xLength, yLength, zLength});
}

void Proxigram::setupIsosurface()
{
    std::size_t isosurfaceNumber = gs.getCompsList().size();
    const auto& targetIso = gs.getCompsList();
    isos.resize(isosurfaceNumber);
    for(std::size_t i = 0; i < isosurfaceNumber; ++i)
    {
        isos[i].setUp(targetIso[i]);
    }
}

void Proxigram::computeIsoAtoms()
{
    StateMachineReader smr("Voxel Draw");
    const State& state = smr.getState();
    if(!state.EnvironmentLoaded())
        return;
    //// Опасный момент на грани темной магии
    State& st = smr.UnsafeGetState();
    AtomsVector& atoms = st.GetAtoms();

    const auto& cL = gs.getCompsList();

    std::size_t SIZE = cL.size();
    size_t at_ind = 0;
    st.CreateClustersData();
    st.clusters().SetSearchParams(ForClusters);
    st.clusters().resize(isos.size()+1);
    st.CalculateSubobjects();
    for(const auto& atom: atoms)
    {

        Point p(atom.fPx, atom.fPy, atom.fPz);
        double distance = proxi[cL[0][0]].first.ComputeDistance(p);
        QVector3D normal = proxi[cL[0][0]].second;
        Triangle tri = proxi[cL[0][0]].first;
        std::size_t index = 0;
        //Для всех изоповерхностей находим самый близкий к атому треугольник
        for(std::size_t i = 0; i < SIZE; ++i)
        {
            std::size_t isoSize = cL[i].size();
            //Для всех треугольников изоповерхности
            for(std::size_t k = 0; k < isoSize; ++k)
            {
                double d = proxi[cL[i][k]].first.ComputeDistance(p);
                if(d < distance)
                {
                    std::size_t idx = cL[i][k];
                    distance = d;
                    index = i;
                    normal = proxi[cL[i][k]].second;
                    tri = proxi[idx].first;
                    //break;
                }
            }
         }
        //Если атом внутри, то увеличиваем счетчик атомов изоповерхности нужно индекса
        QVector3D ptmp(atom.fPx - tri.vertex[0]->X(), atom.fPy - tri.vertex[0]->Y() , atom.fPz - tri.vertex[0]->Z());
//        Atom* at = const_cast<Atom*>(&atom);
//        State* st = const_cast<State*>(&state);



        if(QVector3D::dotProduct(ptmp, normal) <= 0)
        {
            isos[index].incrementAtomsNumber();
            atoms[at_ind].ClusterId = index+1;
            atoms[at_ind].AtomType =  Atom::Base;
        }
        else
        {
            atoms[at_ind].ClusterId = 0;
            atoms[at_ind].AtomType =  Atom::Matrix;
        }
        //std::cout << atom.ClusterId << std::endl;
        at_ind++;
    }

}
void Proxigram::on_lineEdit_returnPressed()
{
    for(size_t i=0; i<ui->IsosurfaceIDList->rowCount();++i)
    {
        ui->IsosurfaceIDList->item(i,1)->setCheckState(Qt::Unchecked);
    }
    showID.clear();
    int value = ui->lineEdit->text().toDouble();
    ui->horizontalSlider->setValue(value*10);
    for(std::size_t i = 0; i < isos.size(); ++i)
    {
        if(isos[i].getLength() >= value)
        {
            showID.push_back(i);
            ui->IsosurfaceIDList->item(i,1)->setCheckState(Qt::Checked);
        }
    }
    for(std::size_t i = 0; i < showID.size(); ++i)
    {
        std::cout << "ID selected = " << showID[i] << std::endl;

    }

}

void Proxigram::on_horizontalSlider_valueChanged(int value)
{
    ui->lineEdit->setText(QString::number(value/10));
}

void Proxigram::on_horizontalSlider_sliderReleased()
{
    on_lineEdit_returnPressed();
}
