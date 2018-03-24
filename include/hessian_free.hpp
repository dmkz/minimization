#pragma once

#include "math.hpp"
#include "iteration_object.hpp"

// Апроксимация умножения градиента функции f в точке х на вектор dx
// Погрешность O(||h||^2*||dx||^3), где h - выбранный шаг дифференцирования
// Авторы: Козырев Дмитрий, Бадави Полина
Real grad_prod_vect(Function f, const Vector& x, const Vector& dx);

// Апроксимация умножения матрицы Гессе в точке x на вектор dx 
// Погрешность O(||h||^2*||dx||^3)
// Авторы: Козырев Дмитрий, Бадави Полина
Vector hess_prod_vect(Function f, const Vector& x, const Vector& dx);

// Метод сопряженных градиентов
// Авторы: Козырев Дмитрий, Бадави Полина
Vector conjugade_gradient(Matrix A, Vector b, Vector x);

// Медленный алгоритм Hessian Free (без модификации)
// Авторы: Козырев Дмитрий, Бадави Полина
void slow_hessian_free(Function f, Vector x, BasicIterationObject* iter_object);
// f - указатель на целевую функцию
// x - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации

// Быстрый алгоритм Hessian Free (с модификацией)
// Авторы: Козырев Дмитрий, Бадави Полина
void hessian_free(Function f, Vector x, BasicIterationObject* iter_object);
// f - указатель на целевую функцию
// x - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации