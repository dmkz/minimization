#pragma once

#include "math.hpp"

// Апроксимация умножения градиента функции f в точке х на вектор dx
// Погрешность O(||h||^2*||dx||^3), где h - выбранный шаг дифференцирования
// Авторы: Козырев Дмитрий, Бадави Полина
ld grad_prod_vect(Function f, const Vector& x, const Vector& dx);

// Апроксимация умножения матрицы Гессе в точке x на вектор dx 
// Погрешность O(||h||^2*||dx||^3)
// Авторы: Козырев Дмитрий, Бадави Полина
Vector hess_prod_vect(Function f, const Vector& x, const Vector& dx);

// Метод сопряженных градиентов
// Авторы: Козырев Дмитрий, Бадави Полина
Vector conjugade_gradient(Matrix A, Vector b, Vector x);

// Медленный алгоритм Hessian Free (без модификации)
// Авторы: Козырев Дмитрий, Бадави Полина
std::pair<Vector, int> slow_hessian_free(Function f, Vector x, int iter_limit = 100);

// Быстрый алгоритм Hessian Free (с модификацией)
// Авторы: Козырев Дмитрий, Бадави Полина
std::pair<Vector, int> hessian_free(Function f, Vector x, int iter_limit = 100);