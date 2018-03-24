#pragma once

#include "math.hpp"
#include "StopCondition.hpp"

// Авторы: Михнев Денис (реализация), Грошева Екатерина (теория)

IterationData bfgs(Function f, Vector start_point, const StopCondition& stop_condition = default_stop_condition);
// f - указатель на целевую функцию
// start_point - начальное приближение
// stop_condition - критерий остановы
// Результат работы метода будет лежать в структуре данных о последней итерации

Matrix out_pr_bfgs(Vector& x, Vector& y);

Matrix hes_upd_bfgs(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

Real search_alpha_bfgs(Function f, Vector& x, Vector& p, int iter_limit);