#pragma once

// #include "math.hpp"
// #include <cmath>  
// #include <math.h>
// #include <limits>

// Авторы: Бабичев Денис (теория), Бессонов Трофим (тестирование, теория), Данилов Алексей (реализация), Киселев Николай (реализация)

#include <fstream>
#include <complex>
#include <iostream>
#include <cmath>
#include "math.hpp"


ld brent(const ld ax, const ld bx, const ld cx, ld f(const ld, int, Vector *, Vector *, Function),
	const ld tol, int ncom, Vector *pcom_p, Vector *xicom_p, Function);

ld f1dim(const ld x, int ncom, Vector *pcom_p, Vector *xicom_p, Function func);


std::pair<Vector, Vector> linmin(Vector &pInit, Vector &xiInit, Function func);

std::pair<Vector, Vector> mnbrak(ld axInit, ld bxInit, ld func(const ld, int, Vector *, Vector *, Function), int ncom, Vector *pcom_p, Vector *xicom_p, Function);
      
std::pair<Vector, int> powell(Function func, Vector p, int iter_limit);