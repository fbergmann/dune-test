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


#endif // DATAHELPER_HH
