#include "proxigramstools.h"
#include <stack>

std::vector<double> getSumConcentration(const std::vector<std::vector<double>>& plotData)
{
    std::vector<double> sumConcetrationVector;
    sumConcetrationVector.resize(plotData[0].size());

    double sumConc = 0.0;

    for(size_t i = 0; i < plotData[0].size(); ++i)
    {
        for(size_t k = 0; k < plotData.size(); ++k)
        {
            sumConc += plotData[k][i];
        }
        sumConcetrationVector[i] = sumConc;
        sumConc = 0.0;
    }
    return sumConcetrationVector;
}

bool isTriangleConnect(const Triangle &first, const Triangle &second)
{
    int count = 0;
    if(*first.vertex[0] == *second.vertex[0] ||
       *first.vertex[0] == *second.vertex[1] ||
       *first.vertex[0] == *second.vertex[2]) count++;

    if(*first.vertex[1] == *second.vertex[0] ||
       *first.vertex[1] == *second.vertex[1] ||
       *first.vertex[1] == *second.vertex[2]) count++;

    if(*first.vertex[2] == *second.vertex[0] ||
       *first.vertex[2] == *second.vertex[1] ||
       *first.vertex[2] == *second.vertex[2]) count++;

    if(count == 2) return true;
    if(count == 3) return false;
    return false;
}

//proxTriangles implementation =========================================================

proxiTriangles::proxiTriangles():
    needUpdate(false)
{
    triangles.resize(0);
}

std::pair<Triangle, proxiTriangles::Normal> proxiTriangles::operator[](size_t index)
{
    //checkIndex
    assert(index == triangles.size());
    return triangles[index];
}

const std::pair<Triangle, proxiTriangles::Normal> proxiTriangles::operator[](size_t index) const
{
    //checkIndex
    assert(index == triangles.size());
    return triangles[index];
}

proxiTriangles::TrianglesVectorIterator proxiTriangles::begin()
{
    return std::begin(triangles);
}

proxiTriangles::TrianglesVectorConstIterator proxiTriangles::cbegin() const
{
    return std::cbegin(triangles);
}

proxiTriangles::TrianglesVectorIterator proxiTriangles::end()
{
    return std::end(triangles);
}

proxiTriangles::TrianglesVectorConstIterator proxiTriangles::cend() const
{
    return std::cend(triangles);
}

void proxiTriangles::buildTriangle(const std::vector<Point> &pointSet,
                                   const std::vector<proxiTriangles::Normal> &normals)
{
    if(needUpdate) triangles.clear();
    bool cheak = pointSet.size() % 3 != 0;
    if(cheak) return;

    std::size_t SIZE = pointSet.size();
    for(std::size_t i = 0; i < SIZE; i += 3)
    {
        Triangle t(pointSet[i], pointSet[i+1], pointSet[i+2]);
        triangles.emplace_back(std::make_pair(std::ref(t), normals[i]));
    }

    needUpdate = true;
}

//proxDatagram implementation =========================================================

proxiDatagram::proxiDatagram(size_t NumBins, double xmin, double xmax):
    NumBins(NumBins), xmin(xmin), xmax(xmax)
{
    update(NumBins, xmin, xmax);
}

void proxiDatagram::update(size_t NumBins, double xmin, double xmax)
{
    if(needUpdate)
        plotData.clear();
    this->NumBins = NumBins;
    this->xmin = xmin;
    this->xmax = xmax;
    plotData.resize(NumBins, 0.0);
    needUpdate = true;
}

const double& proxiDatagram::operator[](size_t index) const
{
    //check index
    return plotData[index];
}

double& proxiDatagram::operator[](size_t index)
{
    //check index
    return plotData[index];
}

std::vector<double> proxiDatagram::getPlotData() const
{
    return plotData;
}

void proxiDatagram::fillValue(double value, float increment)
{
    if(value < xmin || value > xmax) return;
    std::size_t position = abs((xmin - value)*NumBins / (xmax - xmin));
    if(position == plotData.size())
        position -= 1;
    plotData[position] += increment;
}

proxiDatagram::dataIterator proxiDatagram::dataBegin()
{
    return std::begin(plotData);
}

proxiDatagram::dataConstIterator proxiDatagram::dataCBegin() const
{
    return std::cbegin(plotData);
}

proxiDatagram::dataIterator proxiDatagram::dataEnd()
{
    return std::end(plotData);
}

proxiDatagram::dataConstIterator proxiDatagram::dataCEnd() const
{
    return std::cend(plotData);
}

//GraphSeparator implementation=============================================
GraphSeparator::GraphSeparator()
{
    update(0);
}

GraphSeparator::GraphSeparator(std::size_t numVertex): numVertex(numVertex)
{
    update(numVertex);
}

void GraphSeparator::findComps(std::size_t index = 0)
{
    for(auto& element: used)
        element = false;
    for(std::size_t i = 0; i < numVertex; ++i)
    {
        if(!used[i])
        {
            components.clear();
            if(index == 0)
                dfs(i);
            else
                nonRecursDfs(i);
            compStore.push_back(components);
        }
    }
}

void GraphSeparator::inputEdge(size_t iVertex1, size_t iVertex2)
{
    if(iVertex1 >= numVertex || iVertex2 >= numVertex)
        return;

    auto first = std::begin(graphData[iVertex1]);
    auto last = std::end(graphData[iVertex1]);
    if(std::find(first, last, iVertex2) != last)
        return;

    graphData[iVertex1].push_back(iVertex2);
    graphData[iVertex2].push_back(iVertex1);
}

const GraphSeparator::connectivityList& GraphSeparator::getCompsList() const
{
    return compStore;
}

GraphSeparator::connectivityList& GraphSeparator::getCompsList()
{
     return compStore;
}

void GraphSeparator::dfs(size_t v)
{
    used[v] = true;
    components.push_back(v);
    for (size_t i = 0; i < graphData[v].size(); ++i)
    {
        int to = graphData[v][i];
        if (!used[to])
            dfs(to);
    }
}

void GraphSeparator::nonRecursDfs(std::size_t v)
{
    std::stack<std::size_t> st;
    st.push(v);
    used[v] = true;
    components.push_back(v);
    while(!st.empty())
    {
        std::size_t ver = st.top();
        st.pop();
        for(std::size_t i = 0; i < graphData[ver].size(); ++i)
        {
            if(!used[graphData[ver][i]])
            {
                st.push(graphData[ver][i]);
                used[graphData[ver][i]] = true;
                components.push_back(graphData[ver][i]);
            }
        }
    }
}

void GraphSeparator::update(size_t newNumVertex)
{
    numVertex = newNumVertex;
    used.clear();
    used.resize(numVertex);
    for(auto& element: used)
    {
        element = false;
    }
    graphData.clear();
    graphData.resize(numVertex);
    components.clear();
    compStore.clear();
}


Isosurface::Isosurface()
{

}

void Isosurface::setUp(const std::vector<size_t> &triangles)
{
    this->triangles = triangles;
}

void Isosurface::setLength(double length)
{
    this->length = length;
}

void Isosurface::incrementAtomsNumber()
{
    atomsNumber += 1;
}

size_t Isosurface::getAtomsNumber() const
{
    return atomsNumber;
}

double Isosurface::getLength() const
{
    return length;
}

size_t Isosurface::getSize() const
{
    return triangles.size();
}

size_t Isosurface::operator [](size_t index) const
{
    return triangles[index];
}
