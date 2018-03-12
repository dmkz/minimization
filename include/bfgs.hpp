#pragma once

#include "math.hpp"
#include "iteration_object.hpp"

// Авторы: Михнев Денис (реализация), Грошева Екатерина (теория)

Matrix out_pr_bfgs(Vector& x, Vector& y);

Matrix hes_upd_bfgs(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

Real search_alpha_bfgs(Function f, Vector& x, Vector& p, int iter_limit);

void bfgs(Function f, Vector start_point, BasicIterationObject* iter_object);
// f - указатель на целевую функцию
// start_point - начальное приближение
// iter_object - объект итерации
// Результат работы метода будет лежать в объекте итерации