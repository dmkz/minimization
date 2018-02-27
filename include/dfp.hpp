#pragma once

#include "math.hpp"

// Авторы: Амиров Аскар, Китаев Станислав

Matrix out_pr(Vector& x, Vector& y);

Matrix hes_upd(Function f, Matrix& B, Vector& x_cur, Vector& x_prv);

ld search_alpha(Function f, Vector& x, Vector& p, int iter_limit);

std::pair<Vector, int> dfp(Function f, Vector start_point, int iter_limit = 100);
