#pragma once

#include "math.hpp"
#include "StopCondition.hpp"

// Быстрый (модифицированный) алгоритм Hessian Free
// Авторы: Козырев Дмитрий, Бадави Полина
IterationData hessian_free(Function f, Vector x, const StopCondition& stop_condition = default_stop_condition);
// f - указатель на целевую функцию
// x - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации

// Медленный алгоритм Hessian Free (без модификации)
// Авторы: Козырев Дмитрий, Бадави Полина
IterationData slow_hessian_free(Function f, Vector x, const StopCondition& stop_condition = default_stop_condition);
// f - указатель на целевую функцию
// x - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации

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