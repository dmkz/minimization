/*
 * Авторы: Бураханова А., Золкин А.
 */

#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cassert>

typedef long double Real;

template<typename T> struct Point
{
    uint32_t N;
    std::vector<T> coordinate;

    Point();

    Point(uint32_t N_);

    Point(uint32_t N_, const std::vector<T>& _coordinate);
};

template <typename T>
Point<T>::Point()
{
    N = 1;
	coordinate.resize(N);
}

template <typename T>
Point<T>::Point(uint32_t N_)
{
    N = N_;
	coordinate.resize(N_);
}

template <typename T>
Point<T>::Point(uint32_t N_, const std::vector<T>& _coordinate)
{
    N = N_;
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
    virtual int Init(uint32_t N_, uint32_t D_, std::string dir_file) = 0;
	virtual int Reset() = 0;
	
    // Генерация одной точки
    virtual PointReal GeneratePoint() = 0;
};
