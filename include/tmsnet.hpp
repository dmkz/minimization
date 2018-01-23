#pragma once

#include <string>

typedef long double Real;

template<typename T> struct Point
{
    uint32_t N;
    std::vector<T> coordinate;

    Point();

    Point(uint32_t _N);

    Point(uint32_t _N, const std::vector<T>& _coordinate);
};

template <typename T>
Point<T>::Point()
{
    N = 1;
	coordinate.resize(N);
}

template <typename T>
Point<T>::Point(uint32_t _N)
{
    N = _N;
	coordinate.resize(_N);
}

template <typename T>
Point<T>::Point(uint32_t _N, const std::vector<T>& _coordinate)
{
    N = _N;
    coordinate = _coordinate;
}


typedef Point<Real> PointReal;

struct NetPoint {
    uint32_t point;
    double value;
};

class TMSNet
{
public:
    virtual int Init() = 0;
    virtual int Init(uint32_t _N, uint32_t _D, std::string dir_file) = 0;
	virtual int Reset() = 0;
	
    // Генерация одной точки
    virtual PointReal GeneratePoint() = 0;
};
