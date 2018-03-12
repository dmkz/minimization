#pragma once

#include "math.hpp"
#include "iteration_object.hpp"

// Авторы: Амиров Аскар, Китаев Станислав

Matrix out_pr(Vector& x, Vector& y);

Matrix hes_upd(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

Real search_alpha(Function f, const Vector& x, const Vector& p, int iter_limit);

void dfp(Function f, Vector start_point, BasicIterationObject* iter_object);
// f - указатель на целевую функцию
// start_point - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации