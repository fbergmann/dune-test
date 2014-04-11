#ifndef DATAHELPER_HH
#define DATAHELPER_HH

#include <vector>
#include <string>
#include <fstream>
#include <math.h>

class DataHelper
{

private:
    std::vector<std::vector<double>> mData;
    int maxX;
    int maxY;
    bool valid;

public:
    DataHelper()
        : mData(0)
        , maxX(0)
        , maxY(0)
        , valid(false)
    {
    }

    DataHelper(const std::string& fileName)
        : mData(0)
        , maxX(0)
        , maxY(0)
        , valid(false)
    {
        initializeFromFile(fileName);
    }

    void initializeFromFile(const std::string& fileName)
    {
        std::ifstream stream(fileName.c_str(), std::ios_base::binary);
        valid = stream.good();
        if (!valid) return;
        stream >> maxX >> maxY;
        mData.resize(maxX);
        for (int x = 0; x < maxX; ++x)
        {
            mData[x].resize(maxY);
            for (int y = 0; y < maxY; ++y)
            {
                stream >> mData[x][y];
            }
        }
        stream.close();
    }

    bool isValid() const { return valid; }

    double operator ()(double x, double y) const
    {
        int xpos = floor(x);
        if (xpos > maxX) xpos = maxX - 1;
        if (xpos < 0) xpos = 0;
        int ypos = floor(y);
        if (ypos > maxY) ypos = maxY - 1;
        if (ypos < 0) ypos = 0;

        return mData[xpos][ypos];
    }

};


#endif // DATAHELPER_HH
