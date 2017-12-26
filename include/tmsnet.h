#pragma once

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
    // Генерация одной точки
    PointReal GeneratePoint()
    {
        return PointReal();
    };
};
