#pragma once

#include "math.hpp"

// Авторы: Михнев Денис (реализация), Грошева Екатерина (теория)

Matrix out_pr_bfgs(Vector& x, Vector& y);

Matrix hes_upd_bfgs(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

ld search_alpha_bfgs(Function f, Vector& x, Vector& p, int iter_limit);

std::pair<Vector, int> bfgs(Function f, Vector start_point, int iter_limit = 100);