#pragma once

/*
    Библиотека со всей математикой. 
    Автор: Дмитрий Козырев
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cassert>

#ifndef M_PI
    #define M_PI 3.14159265358979323846 
#endif

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifndef M_E
    #define M_E 2.71828182845904523536
#endif

// Объявление типов:
typedef long double Real;
typedef std::vector<Real> Vector;
typedef std::vector<Vector> Matrix;
typedef Real(*Function)(const Vector & x);

// Точность сравнения вещественных чисел:
const Real COMPARE_EPS = 0.0000000000000001L;

// Математические константы (M_E, M_PI, M_PI_2) в файле реализации

// Операторы вывода в поток:
std::ostream& operator<<(std::ostream& os, const Vector& v);
std::ostream& operator<<(std::ostream& os, const Matrix& m);

// Операция умножения вектора на число и различные ее вариации:
Vector& operator*=(Vector& v, const Real value);
Vector operator*(const Real value, Vector v);
Vector operator*(Vector v, const Real value);

// Унарный минус для вектора:
Vector operator-(Vector v);

// Сложение и вычитание векторов:
Vector& operator+=(Vector & v1, const Vector& v2);
Vector& operator-=(Vector & v1, const Vector& v2);
Vector operator+(Vector v1, const Vector& v2);
Vector operator-(Vector v1, const Vector& v2);

// Скалярное произведение векторов:
Real dot(const Vector& v1, const Vector& v2);

// Норма вектора:
Real norm(const Vector& v);

// Умножение матрицы на вектор:
Vector operator*(const Matrix& m, const Vector& v);

// Проверка равенства вектора нулю
bool is_zero(const Vector& v);

// Базисный вектор в пространстве R^n:
// на i-ом месте 1, на всех остальных нули
Vector id_vect(int size, int i);

// Явное взятие градиента от функции:
Vector grad(Function f, const Vector& point, const Real h = 1e-4);

// Явное вычисление матрицы Гессе (слишком затратно, годится только для проверки результатов)
// в точке x с шагом h и погрешностью O(h)
Matrix hess(Function f, const Vector& x, const Real h = 1e-4);


template<class T>
inline const T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

template<class T>
inline void SWAP(T &a, T &b)
{
	T dum = a; a = b; b = dum;
}

template<class T>
inline const T MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}

template<class T>
inline const T SQR(const T a) { return a * a; }