#ifndef DATAHELPER_HH
#define DATAHELPER_HH

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

enum InterpolationType
{
  NearestNeighbor
  , Bilinear
//  , Bicubic
};

class DataHelper
{

private:
  std::vector<std::vector<double>> mData;
  int dimY;
  int dimX;
  double minX;
  double maxX;
  double minY;
  double maxY;
  bool valid;
  InterpolationType mType;

public:

  DataHelper()
    : mData(0)
    , dimY(0)
    , dimX(0)
    , valid(false)
    , mType(NearestNeighbor)
  {
  }

  static DataHelper* forFile(const std::string& fileName, InterpolationType type = Bilinear)
  {
    DataHelper* result = new DataHelper(fileName, type);
    if (!result->isValid())
      {
        delete result;
        result = NULL;
      }
    return result;
  }
  
  DataHelper(const std::string& fileName, InterpolationType type = Bilinear)
    : mData(0)
    , dimY(0)
    , dimX(0)
    , minX(0)
    , maxX(0)
    , minY(0)
    , maxY(0)
    , valid(false)
    , mType(type)
  {
    initializeFromFile(fileName);
  }

  DataHelper(int numRows, int numCols, const double &value=double() )
    : mData(numRows)
    , dimY(numRows)
    , dimX(numCols)
    , minX(0)
    , maxX(numCols)
    , minY(0)
    , maxY(numRows)
    , valid(true)
    , mType(Bilinear)
  {
    if ( ((numCols==0) && (numRows!=0)) || ((numCols!=0) && (numRows==0)) )
      {
        dimY=0;
        dimX=0;
        mData.resize(dimY);
      }
    else
      {
        for (int i=0;i<dimY;++i)
          mData[i].resize(dimX);
        dimY=numRows;
        dimX=numCols;
      }
    for (int i=0;i<numRows;++i)
      {
        for (int j=0;j<numCols;++j)
          mData[i][j]=value;
      }
  }

  void resize(int numRows, int numCols, const double &value = double() )
  {
    mData.resize(numRows);
    for (size_t i=0;i<mData.size();++i)
      {
        mData[i].resize(numCols);
        for (size_t j=0;j<mData[i].size();++j)
          mData[i][j]=value;
      }
    dimY=numRows;
    dimX=numCols;
    valid = true;
  }


  void initializeFromFile(const std::string& fileName)
  {
    std::ifstream stream(fileName.c_str(), std::ios_base::binary);
    valid = stream.good();
    if (!valid) return;
    stream >> dimX >> dimY;
    stream >> minX >> maxX >> minY >> maxY;
    mData.resize(dimY);
    for (int y = 0; y < dimY; ++y)
      {
        mData[y].resize(dimX);
        for (int x = 0; x < dimX; ++x)
          {
            stream >> mData[y][x];
          }
      }
    stream.close();
  }

  void writeToFile(const std::string& fileName) const
  {
    std::ofstream data(fileName.c_str(), std::ios_base::binary);

    data << dimX << " ";
    data << dimY << std::endl;

    for (int y = 0; y < dimY; ++y)
      {
        for (int x = 0; x < dimX; ++x)
          {
            data << mData[y][x] << " ";
          }
        data << std::endl;
      }

    data.flush();
    data.close();
  }

  bool isValid() const { return valid; }

  int xIndexForDouble(double x) const
  {
    if (x < minX)
      x = minX;
    if (x > maxX)
      x = maxX;
    int index = floor(((x - minX)/maxX)*(double)dimX);
    if (index < 0) index = 0;
    if (index >= dimX) index = dimX -1;
    return index;
  }

  int yIndexForDouble(double y) const
  {
    if (y < minY)
      y = minY;
    if (y > maxY)
      y = maxY;
    int index = floor(((y - minY)/maxY)*(double)dimY);
    if (index < 0) index = 0;
    if (index >= dimY) index = dimY -1;
    return index;
  }

  double xDataForIndex(int index) const
  {
    double value =
        ((double)index/(double)dimX)*(maxX-minX)
        + minX;
    return value;
  }

  double yDataForIndex(int index) const
  {
    double value =
        ((double)index/(double)dimY)*(maxY-minY)
        + minY;
    return value;
  }

  double operator ()(int x, int y) const
  {
    if (x >= dimX) x=dimX -1;
    if (x < 0) x = 0;
    if (y >= dimY) y=dimY -1;
    if (y < 0) y = 0;

    return mData[y][x];
  }


  double &operator()(int x, int y)
  {
    if (x >= dimX) x=dimX -1;
    if (x < 0) x = 0;
    if (y >= dimY) y=dimY -1;
    if (y < 0) y = 0;

    return mData[y][x];
  }

  double get(double x, double y) const
  {
    int xpos = xIndexForDouble(x);
    int ypos = yIndexForDouble(y);

    switch (mType) {
      case Bilinear:
        {
          return bilinearAround(
                xpos, ypos,
                x,
                y);
        }
      default:
        return mData[ypos][xpos];
        break;
      }
  }

  int Rows() const
  {
    return dimY;
  }

  int Cols() const
  {
    return dimX;
  }

  double getMinX() const
  {
    return minX;
  }

  double getMaxX() const
  {
    return maxX;
  }

  double getMinY() const
  {
    return minY;
  }

  double getMaxY() const
  {
    return maxY;
  }

  InterpolationType getInterpolationType() const
  {
    return mType;
  }

  void setInterpolationType(InterpolationType type)
  {
    mType = type;
  }

  double bilinearAround(
                        double x,
                        double y) const
  {
    int xpos = xIndexForDouble(x);
    int ypos = yIndexForDouble(y);
    return bilinearAround(
          xpos, ypos,
          x,
          y);
  }

  double bilinearAround(int xIndex, int yIndex,
                        double x,
                        double y) const
  {
    double x1 = xDataForIndex(xIndex);
    double y1 = yDataForIndex(yIndex);
    double diffXX1 = x - x1;
    double diffYY1 = y - y1;
    int xNeighbor = xIndex +
        (std::signbit(diffXX1) ? -1 : 1);
    if (xNeighbor  >= dimX) xNeighbor = dimX -1;
    if (xNeighbor < 0) xNeighbor = 0;
    int yNeighbor = yIndex +
        (std::signbit(diffYY1) ? -1 : 1);
    if (yNeighbor  >= dimY) yNeighbor = dimY -1;
    if (yNeighbor < 0) yNeighbor = 0;

    double x2 = xDataForIndex(xNeighbor);
    double y2 = yDataForIndex(yNeighbor);


    //    cout << "(" << xIndex << " " << yIndex << ") "
    //         << "(" << xNeighbor << " " << yNeighbor << ") "
    //         << endl;

    double x1y1 = mData[yIndex][xIndex];
    if (x2 == x1 || y2 == y1) return  x1y1 ;
    double x2y1 = mData[yIndex][xNeighbor];
    double x1y2 = mData[yNeighbor][xIndex];
    double x2y2 = mData[yNeighbor][xNeighbor];


    return
          ((y2-y)/(y2-y1))
        * ((x2-x)/(x2-x1)*x1y1 + (x-x1)/(x2-x1)*x2y1)
        + ((y-y1)/(y2-y1))
        * ((x2-x)/(x2-x1)*x1y2 + (x-x1)/(x2-x1)*x2y2);

  }

  void setBounds(double xmin, double xmax, double ymin, double ymax)
  {
    minX = xmin;
    maxX = xmax;
    minY = ymin;
    maxY = ymax;
  }
};



template<typename GV, typename RF>
void getCoordinates(GV &gv, std::vector< std::pair<RF,RF> > &coordinates)
{
    coordinates.clear();
    for(auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
    {
        const auto& center = it->geometry().center();

        coordinates.push_back(std::make_pair((RF)center[0], (RF)center[1]));
    }
}

template<typename V, typename RF>
void  setNthElement (V &v, size_t n, size_t numComponents, const std::vector< std::pair<RF,RF> > &coordinates, const DataHelper* data)
{
    if (data == NULL) return;

    size_t numCoords = coordinates.size();
    size_t currentPos = 0;

    size_t pos= 0;
    for (auto it=v.begin(); it!=v.end();)
    {
        // skip first n components
        size_t skip = 0;
        while (skip < n)
        {
            ++it;
            ++skip;
            ++pos;
        }

        RF x = coordinates[currentPos].first;
        RF y = coordinates[currentPos].second;

        // change value position
        *it = (*data).get(x,y);

        ++it;
        ++skip;
        ++pos;


        // skip remaining components
        while(skip < numComponents)
        {
            ++it;
            ++skip;
            ++pos;
        }


        // advance position
        ++currentPos;

        if (currentPos == numCoords)
            currentPos = 0;
    }
}


template<typename V>
double calculateDifference (V &vnew, V &old)
{
  double result = 0;
  for (auto it=vnew.begin(), it2 = old.begin(); it!=vnew.end();++it, ++it2)
  {
    result += pow((*it) - (*it2), 2) ;
  }
  return sqrt(result);
}



class EventData
{
public:
    EventData(const Dune::ParameterTree& param)
     : mData(DataHelper::forFile(param.get<std::string>("file")))
     , mStartTime(param.get<double>("start", 0))
     , mAppliedEvent(false)
     , mVariableIndex(param.get<int>("target"))
     , mUniformValue(param.get<double>("uniform", 0))
    {

    }

    bool hasFired() const { return mAppliedEvent; }

    void setAppliedEvent(bool appliedEvent) { mAppliedEvent = appliedEvent; }

    int getVariableIndex() const { return mVariableIndex; }

    bool hasData() const { return mData != NULL; }

    const DataHelper* getData() const { return mData; }

    double getUniformValue() const { return mUniformValue; }

    bool shouldFire(double time) const { return !mAppliedEvent &&  time >= mStartTime;  }

protected:
    DataHelper* mData;
    double mStartTime;
    int mVariableIndex;
    bool mAppliedEvent;
    double mUniformValue;
};

#define EVENT_SUPPORT_INITIALIZE(gv)\
  std::vector< std::pair<RF, RF> > coordinates;\
  getCoordinates(gv, coordinates);\


#define EVENT_SUPPORT_HANDLE_EVENTS(time, numVariables, events, uold, unew )\
    {\
      std::vector<EventData>::iterator eventIt = events.begin();\
      while(eventIt != events.end())\
      {\
        if (eventIt->shouldFire(time))\
        {\
            setNthElement(unew, eventIt->getVariableIndex(), numVariables, coordinates, eventIt->getData());\
            uold = unew;\
            eventIt->setAppliedEvent(true);\
            appliedEvent = true;\
            std::cout << "applied event!" << std::endl;\
        }\
        ++eventIt;\
      }\
    }\

#define EVENT_SUPPORT_READ_DATA(target, parser, filename, skipVariable)\
  {\
    try\
    {\
        Dune::ParameterTree eventData;\
        parser.readINITree(filename, eventData);\
        \
        skipVariable = eventData.get<bool>("skipOutputUntilEvent", false);\
        \
        const auto& keys = eventData.getSubKeys();\
        for (auto it= keys.begin(); it != keys.end(); ++it )\
        {\
            target.push_back(EventData(eventData.sub(*it)));\
        }\
    }\
    catch(...)\
    {\
    }\
  }\

#define EVENT_SUPPORT_GLOBALS(mEvents, skipOutputUntilEvent, appliedEvent)\
    std::vector<EventData> mEvents;\
    bool skipOutputUntilEvent = false;\
    bool appliedEvent = false\


#endif // DATAHELPER_HH
