#pragma once

#include "math.hpp"

// Апроксимация умножения градиента функции f в точке х на вектор dx
// Погрешность O(||h*dx||^2), где h - выбранный шаг дифференцирования
ld grad_prod_vect(Function f, const Vector& x, const Vector& dx);

// Апроксимация умножения матрицы Гессе в точке x на вектор dx 
// Погрешность O(||h*dx||^2)
Vector hess_prod_vect(Function f, const Vector& x, const Vector& dx);

// Метод сопряженных градиентов
Vector conjugade_gradient(Matrix A, Vector b, ld c, Vector x);

// Медленный алгоритм Hessian Free (без модификации)
std::pair<Vector, int> slow_hessian_free(Function f, Vector x, int iter_limit = 100);

// Быстрый алгоритм Hessian Free (с модификацией)
std::pair<Vector, int> hessian_free(Function f, Vector x, int iter_limit = 100);