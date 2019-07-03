#ifndef PROXIGRAMSTOOLS_H
#define PROXIGRAMSTOOLS_H

#include "geometry.h"
#include "deloneytriangulation.h"
#include "VoxelDraw.h"
#include <cassert>

///TODO: пишем функцию сепаратора в проксиграмме
///Тестируем сепаратор на тестовых данных
///Пишем класс изоповерхности, в которой должны быть ее треугольники и атомы попавшие
///внутрь

std::vector<double> getSumConcentration(const std::vector<std::vector<double>>& plotData);
bool isTriangleConnect(const Triangle& first, const Triangle& second);

class proxiDatagram
{
private:

    using dataIterator = std::vector<double>::iterator;
    using dataConstIterator =std::vector<double>::const_iterator;

    std::size_t NumBins;
    double xmin, xmax;
    std::vector<double> plotData;
    bool needUpdate = false;

public:

    proxiDatagram(std::size_t NumBins = 20, double xmin = -10.0, double xmax = 10.0);

    void update(std::size_t NumBins, double xmin, double xmax);

    //Операторы доступа к данным ==============================================

    std::size_t getSize() const {return plotData.size();}
    const double& operator[](std::size_t index) const;
    double& operator[](std::size_t index);
    std::vector<double> getPlotData() const;
    void fillValue(double value, float increment);

    //Итераторы ===============================================================

    dataIterator dataBegin();
    dataConstIterator dataCBegin() const;
    dataIterator dataEnd();
    dataConstIterator dataCEnd() const;

};

class proxiTriangles
{
public:

    using Normal = QVector3D; //нормаль к треугольнику

    using TrianglesVector = std::vector<std::pair<Triangle, Normal>>;
    using TrianglesVectorIterator =
          std::vector<std::pair<Triangle, Normal>>::iterator;
    using TrianglesVectorConstIterator =
          std::vector<std::pair<Triangle, Normal>>::const_iterator;

    explicit proxiTriangles();

    void buildTriangle(const std::vector<Point>& pointSet,
                       const std::vector<Normal>& normals);

    //Операторы доступа к данным ==============================================
    const std::pair<Triangle, Normal> operator[] (std::size_t index) const;

    std::pair<Triangle, Normal> operator[] (std::size_t index);

    std::size_t getSize() const {return triangles.size();}

    //Итераторы ===============================================================
    TrianglesVectorIterator begin();

    TrianglesVectorConstIterator cbegin() const;

    TrianglesVectorIterator end();

    TrianglesVectorConstIterator cend() const;

private:

    bool needUpdate;

    TrianglesVector triangles;

};

class GraphSeparator
{
public:

    using Normal = QVector3D;
    using connectivityList = std::vector<std::vector<std::size_t>>;

    GraphSeparator();

    explicit GraphSeparator(std::size_t numVertex);

    void findComps(std::size_t index);
    void inputEdge(std::size_t iVertex1, std::size_t iVertex2);
    void update(std::size_t newNumVertex);

    void PrintState() const
    {
        std::cout << "NumVertex = " << numVertex << std::endl;
        std::cout << "GraphData size = " << graphData.size() << std::endl;
        std::cout << "CompStore size = " << compStore.size() << std::endl;
        std::cout << "Components size = " << components.size() << std::endl;
    }

    const connectivityList& getCompsList() const;
    connectivityList& getCompsList();

private:

    void dfs(std::size_t v);
    void nonRecursDfs(std::size_t v);

    std::size_t numVertex;
    QVector<bool> used;
    connectivityList graphData;
    connectivityList compStore;
    std::vector<std::size_t> components;
};

class Isosurface
{
public:

    Isosurface();

    void setUp(const std::vector<std::size_t>& triangles);
    void setLength(double);
    void incrementAtomsNumber();
    std::size_t getAtomsNumber() const;
    double getLength() const;
    std::size_t getSize() const;
    std::size_t operator [](std::size_t index) const;

    std::vector<std::size_t> triangles;

private:

    /// Треугольники изоповерхности
    //std::vector<std::size_t> triangles;
    /// Кол-во атомов попавших в изоповерхность
    std::size_t atomsNumber = 0;
    /// Характерный размер изоповерхности
    double length = 0.0;
};

//Шаблон вычисляющий значение Value для объекта Target по отношению к каждому объекту
//из диапазона last - first
template<class InputIt, class Target, typename Value>
auto computeData(InputIt first, InputIt last, const Target& target,
                 Value result)->decltype(Value)
{
    while(first != last)
    {
        Value tmp = first->computeDistance(target);
        if(tmp < result)
            result = tmp;
        first++;
    }
    return result;
}


#endif // PROXIGRAMSTOOLS_H
