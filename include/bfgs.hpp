#pragma once

#include "math.hpp"

// Авторы: Михнев Денис, Грошева Екатерина
Matrix out_pr(Vector& x, Vector& y);

// Авторы: Михнев Денис, Грошева Екатерина
Matrix hes_upd(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

// Авторы: Михнев Денис, Грошева Екатерина
ld search_alpha(Function f, Vector& x, Vector& p, int iter_limit);

// Авторы: Михнев Денис, Грошева Екатерина
std::pair<Vector, int> bfgs(Function f, Vector start_point, int iter_limit = 100);