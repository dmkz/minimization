#pragma once

typedef long double Real;

using namespace std;

template<typename T> struct Point
{
    long unsigned int N;
    vector<T> coordinate;

    Point();

    Point(long unsigned int _N);

    Point(long unsigned int _N, const vector<T>& _coordinate);
};

template <typename T>
Point<T>::Point()
{
    N = 1;
    coordinate = vector<T>(1, 0);
}

template <typename T>
Point<T>::Point(long unsigned int _N)
{
    N = _N;
    coordinate = vector<T>(_N);
}

template <typename T>
Point<T>::Point(long unsigned int _N, const vector<T>& _coordinate)
{
    N = _N;
    coordinate = _coordinate;
}


typedef Point<Real> PointReal;

struct NetPoint {
    unsigned long point;
    double value;
};

class TMSNet
{
public:
    // Генерация одной точки
    PointReal GeneratePoint()
    {
        return PointReal();
    };
};
